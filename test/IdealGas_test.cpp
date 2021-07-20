/*
 * The MIT License
 *
 * Copyright (c) 2015-2017 The University of Utah
 *
  * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

/*
 * IdealGas_test.cpp
 *
 *  Created on: October 15, 2014
 *      Author: Nathan Yonkee
 */

#include <numeric>
#include <iostream>
#include "TestHelper.h"
#include "LinearMassFracs.h"
#include <pokitt/thermo/Density.h>
#include <pokitt/thermo/Pressure.h>
#include <pokitt/thermo/SpecificVol.h>
#include <pokitt/MixtureMolWeight.h>

#include <expression/ExprLib.h>

#include <spatialops/structured/Grid.h>
#include <spatialops/structured/FieldComparisons.h>
#include <spatialops/util/TimeLogger.h>

#include <boost/foreach.hpp>
#include <boost/program_options.hpp>

#include <cantera/thermo/ThermoPhase.h>

namespace SO = SpatialOps;
typedef SO::SVolField  CellField;
typedef SO::SpatFldPtr<CellField> CellFieldPtrT;

namespace po = boost::program_options;

using namespace pokitt;

//==============================================================================

enum GasQuantity{
  P,
  RHO,
  NU
};

std::string property_name( const GasQuantity q )
{
  std::string name;
  switch(q){
    case P  : name = "pressure";     break;
    case RHO: name = "density";      break;
    case NU : name = "SpecificVol"; break;
  }
  return name;
}

//==============================================================================

const std::vector< std::vector<double> >
extract_mass_fracs( const Expr::TagList yiTags, Expr::FieldManagerList& fml ){
  CellField& yi = fml.field_ref< CellField >(yiTags[0]);
  const size_t nPts = yi.window_with_ghost().glob_npts();
  const int nSpec = yiTags.size();

  std::vector< std::vector<double> > massFracs;
  for( size_t i=0; i<nPts; ++i ){
    std::vector<double> massFrac( nSpec, 0.0 );
    massFracs.push_back(massFrac);
  }
  for( size_t n=0; n<nSpec; ++n ){
    yi = fml.field_ref< CellField >( yiTags[n] );
#   ifdef ENABLE_CUDA
    yi.set_device_as_active( CPU_INDEX );
#   endif
    size_t i=0;
    for( CellField::const_iterator iY = yi.begin(); iY != yi.end(); ++iY, ++i ){
      massFracs[i][n] = *iY;
    }
  }
  return massFracs;
}

//==============================================================================

CellFieldPtrT
get_cantera_result( const bool timings,
                    const size_t canteraReps,
                    const GasQuantity gasQuantity,
                    IdealGasPtr gasMix,
                    Expr::FieldManagerList& fml,
                    const Expr::TagList& yiTags,
                    const Expr::Tag& tTag,
                    const Expr::Tag& refTag )
{
  const std::vector< std::vector<double> > massFracs = extract_mass_fracs( yiTags, fml );

  CellField& temp        = fml.field_ref< CellField >( tTag );
  CellField& refQuantity = fml.field_ref< CellField >( refTag );
# ifdef ENABLE_CUDA
  temp.set_device_as_active       ( CPU_INDEX );
  refQuantity.set_device_as_active( CPU_INDEX );
# endif

  CellFieldPtrT canteraResult = SO::SpatialFieldStore::get<CellField>(temp, CPU_INDEX);

  std::vector< std::vector<double> >::const_iterator iMass;
  CellField::const_iterator                          iTemp;
  CellField::const_iterator                          iRef;
  CellField::iterator                                iCant;
  SO::Timer gasTime;
  gasTime.start();
  for( size_t rep=0; rep < canteraReps; ++rep ){
    iTemp = temp.begin();
    iMass = massFracs.begin();
    iRef  = refQuantity.begin();
    for(iCant = canteraResult->begin(); iCant!=canteraResult->end(); ++iTemp, ++iRef, ++iMass, ++iCant){
      gasMix->setMassFractions_NoNorm( &(*iMass)[0] );
      switch( gasQuantity ){
      case P:
        gasMix->setState_TR( *iTemp, *iRef );
        *iCant=gasMix->pressure();
        break;
      case RHO:
        gasMix->setState_TP( *iTemp, *iRef );
        *iCant=gasMix->density();
        break;
      case NU:
        gasMix->setState_TP( *iTemp, *iRef );
        *iCant=1/gasMix->density();
        break;
      }
    }
  }
  gasTime.stop();
  if( timings ) std::cout << "Cantera " + property_name(gasQuantity) + " time " << gasTime.elapsed_time()/canteraReps << std::endl;
  return canteraResult;
}

//==============================================================================

bool driver( const bool timings,
             const size_t pokittReps,
             const size_t canteraReps,
             const GasQuantity gasQuantity )
{
  TestHelper status( !timings );
  IdealGasPtr gasMix = CanteraObjects::get_thermo();
  const int nSpec = gasMix->nSpecies();

  const Expr::Tag xTag( "XCoord",      Expr::STATE_NONE );
  const Expr::Tag tTag( "Temperature", Expr::STATE_NONE );
  Expr::TagList yiTags;
  for( size_t n=0; n<nSpec; ++n ){
    yiTags.push_back( Expr::Tag( "yi_" + boost::lexical_cast<std::string>(n), Expr::STATE_NONE ) );
  }
  const Expr::Tag mmwTag ("Mixture Molecular Weight", Expr::STATE_NONE);
  Expr::Tag refTag;
  switch( gasQuantity ){
  case P:   refTag = Expr::Tag( "Density",  Expr::STATE_NONE ); break;
  case RHO: refTag = Expr::Tag( "Pressure", Expr::STATE_NONE ); break;
  case NU:  refTag = Expr::Tag( "Pressure", Expr::STATE_NONE ); break;
  }
  const Expr::Tag gasTag (property_name(gasQuantity), Expr::STATE_NONE);

  // we use an initialization tree to avoid recalculations when timing the execution
  Expr::ExpressionFactory initFactory;
  std::set<Expr::ExpressionID> initIDs;
  Expr::ExpressionID tID; // temperature
  Expr::ExpressionID rID; // reference quantity
  Expr::ExpressionID mID; // mixture mol weight
  {
    typedef Expr::PlaceHolder     <CellField>::Builder XCoord;
    typedef Expr::LinearFunction  <CellField>::Builder Temperature;
    typedef Expr::ConstantExpr    <CellField>::Builder RefQuantity;
    typedef       LinearMassFracs <CellField>::Builder MassFracs;
    typedef       MixtureMolWeight<CellField>::Builder MixtureMolWeight;

    initFactory.register_expression(       new XCoord          ( xTag                    ) );
    initFactory.register_expression(       new MassFracs       ( yiTags, xTag            ) );
    tID = initFactory.register_expression( new Temperature     ( tTag,   xTag, 1000, 500 ) );
    mID = initFactory.register_expression( new MixtureMolWeight( mmwTag, yiTags, MASS    ) );
    switch( gasQuantity ){
    case P:   rID = initFactory.register_expression( new RefQuantity( refTag, gasMix->density()  ) ); break;
    case RHO: rID = initFactory.register_expression( new RefQuantity( refTag, gasMix->pressure() ) ); break;
    case NU:  rID = initFactory.register_expression( new RefQuantity( refTag, gasMix->pressure() ) ); break;
    }
    initIDs.insert( tID );
    initIDs.insert( rID );
    initIDs.insert( mID );
  }

  Expr::ExpressionFactory execFactory;
  Expr::ExpressionID execID;
  {
    typedef Expr::PlaceHolder<CellField>::Builder MixtureMolWeight;
    typedef Expr::PlaceHolder<CellField>::Builder RefQuantity;
    typedef Expr::PlaceHolder<CellField>::Builder Temperature;
    typedef       Pressure   <CellField>::Builder Pressure;
    typedef       Density    <CellField>::Builder Density;
    typedef       SpecificVol<CellField>::Builder SpecVol;

    execFactory.register_expression( new MixtureMolWeight ( mmwTag ) );
    execFactory.register_expression( new RefQuantity      ( refTag ) );
    execFactory.register_expression( new Temperature      ( tTag   ) );
    switch( gasQuantity ){
    case P:
      execID = execFactory.register_expression( new Pressure ( gasTag, tTag, refTag, mmwTag) ); break;
    case RHO:
      execID = execFactory.register_expression( new Density  ( gasTag, tTag, refTag, mmwTag) ); break;
    case NU:
      execID = execFactory.register_expression( new SpecVol  ( gasTag, tTag, refTag, mmwTag) ); break;
    }
  }

  Expr::ExpressionTree initTree( initIDs, initFactory, 0 );
  Expr::ExpressionTree execTree( execID , execFactory, 0 );
  {
    std::ofstream gasGraph( "IdealGas_init.dot" );
    initTree.write_tree( gasGraph );
  }
  {
    std::ofstream gasGraph( "IdealGas.dot" );
    execTree.write_tree( gasGraph );
  }

  Expr::FieldManagerList fml;
  initTree.register_fields( fml );
  initTree.bind_fields( fml );
  initTree.lock_fields( fml );

  execTree.register_fields( fml );
  execTree.bind_fields( fml );
  execTree.lock_fields( fml );

  std::vector<SO::IntVec> sizeVec;
  if( timings ){
    sizeVec.push_back( SO::IntVec(2046, 1022, 1) );
    sizeVec.push_back( SO::IntVec(1022, 1022, 1) );
    sizeVec.push_back( SO::IntVec(510,  510,  1) );
    sizeVec.push_back( SO::IntVec(254,  254,  1) );
    sizeVec.push_back( SO::IntVec(126,  126,  1) );
    sizeVec.push_back( SO::IntVec(62,   62,   1) );
    sizeVec.push_back( SO::IntVec(30,   30,   1) );
  }
  else{
    sizeVec.push_back( SO::IntVec( 20,  1,  1) );
  }

  for( std::vector<SO::IntVec>::iterator iSize = sizeVec.begin(); iSize!= sizeVec.end(); ++iSize){

    SO::IntVec gridSize = *iSize;
    fml.allocate_fields( Expr::FieldAllocInfo( gridSize, 0, 0, false, false, false ) );
    SO::Grid grid( gridSize, SO::DoubleVec(1,1,1) );

    CellField& xcoord = fml.field_ref< CellField >( xTag );
    grid.set_coord<SO::XDIR>( xcoord );
    const int nPoints = xcoord.window_with_ghost().glob_npts();
#   ifdef ENABLE_CUDA
    xcoord.set_device_as_active( GPU_INDEX );
#   endif
    initTree.execute_tree();

    if( timings ){
      std::cout << std::endl << property_name(gasQuantity) + " test - " << nPoints << std::endl;
      execTree.execute_tree();// sets memory high-water mark
    }

    SO::Timer timer;
    std::vector< double > times;
    for( size_t rep = 0; rep < pokittReps; ++rep ){
      timer.reset();
      execTree.execute_tree();
      times.push_back( timer.stop() );
    }

    if( timings ){
      std::sort( times.begin(), times.end() );
      const int chop = floor(pokittReps/4);
      const double avgTime = std::accumulate( times.begin() + chop, times.end()-chop, 0.0 )/(pokittReps-2*chop);
      std::cout << "PoKiTT  " + property_name(gasQuantity) + " time " << avgTime << std::endl;
    }

    CellFieldPtrT canteraResult = get_cantera_result( timings, canteraReps, gasQuantity, gasMix, fml, yiTags, tTag, refTag );
    CellField& gasField = fml.field_ref< CellField >( gasTag );
#   ifdef ENABLE_CUDA
    gasField.set_device_as_active( CPU_INDEX );
#   endif

    status( field_equal(gasField, *canteraResult, 1e-12), gasTag.name() );

    fml.deallocate_fields();

  }
  return status.ok();
}

//==============================================================================

int main( int iarg, char* carg[] )
{
  std::string inputFileName;
  std::string inpGroup;
  bool timings = false;
  size_t pokittReps = 1;
  size_t canteraReps = 1;

  // parse the command line options input describing the problem
  po::options_description desc("Supported Options");
  desc.add_options()
               ( "help", "print help message" )
               ( "yaml-input-file", po::value<std::string>(&inputFileName)->default_value("h2o2.yaml"), "Cantera yaml input file name" )
               ( "phase", po::value<std::string>(&inpGroup), "name of phase in Cantera yaml input file" )
               ( "timings", "Generate comparison timings between Cantera and PoKiTT across several problem sizes" )
               ( "pokitt-reps", po::value<size_t>(&pokittReps), "Repeat the PoKiTT tests and report the average execution time")
               ( "cantera-reps", po::value<size_t>(&canteraReps), "Repeat the Cantera tests and report the average execution time");

  try {
    po::variables_map args;
    po::store( po::parse_command_line(iarg,carg,desc), args );
    po::notify(args);

    timings = args.count("timings") || args.count("pokitt-reps") || args.count("cantera-reps");

    if (args.count("help")) {
      std::cout << desc << "\n";
      return -1;
    }
  }

  catch( std::exception& err ){
    std::cout << "Error parsing input arguments\n" << err.what() << std::endl << desc << std::endl;
    return -2;
  }

  try{
    const CanteraObjects::Setup setup( "Mix", inputFileName, inpGroup );
    CanteraObjects::setup_cantera( setup );

    TestHelper status( false );
    status( driver( timings, pokittReps, canteraReps, RHO ) );
    status( driver( timings, pokittReps, canteraReps, P   ) );
    status( driver( timings, pokittReps, canteraReps, NU  ) );

    if( status.ok() ){
      std::cout << "\nPASS\n";
      return 0;
    }
  }
  catch( Cantera::CanteraError& err ){
    std::cout << err.what() << std::endl;
  }
  catch( std::exception& err ){
    std::cout << err.what() << std::endl;
  }

  std::cout << "\nFAIL\n";
  return -1;
}

//==============================================================================
