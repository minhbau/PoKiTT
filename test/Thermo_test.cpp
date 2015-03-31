/*
 * The MIT License
 *
 * Copyright (c) 2015 The University of Utah
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
 * Thermo_test.cpp
 *
 *  Created on: August 28, 2014
 *      Author: Nathan Yonkee
 */

#include <numeric>
#include <iostream>
#include "TestHelper.h"
#include "LinearMassFracs.h"
#include <pokitt/thermo/HeatCapacity_Cp.h>
#include <pokitt/thermo/HeatCapacity_Cv.h>
#include <pokitt/thermo/Enthalpy.h>
#include <pokitt/thermo/InternalEnergy.h>

#include <expression/ExprLib.h>

#include <spatialops/structured/Grid.h>
#include <spatialops/structured/FieldComparisons.h>
#include <spatialops/util/TimeLogger.h>

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>

#include <cantera/IdealGasMix.h>

namespace SO = SpatialOps;
typedef SO::SVolField  CellField;
typedef SO::SpatFldPtr<CellField> CellFieldPtrT;

namespace po = boost::program_options;

using namespace pokitt;

//==============================================================================

enum ThermoQuantity{
  CP,
  CV,
  ENTH,
  E
};

std::string thermo_name( const ThermoQuantity q )
{
  std::string name;
  switch(q){
    case CP  : name = "cp"; break;
    case CV  : name = "cv"; break;
    case ENTH: name = "h";  break;
    case E   : name = "e";  break;
  }
  return name;
}

//==============================================================================

std::set< Expr::ExpressionID >
register_thermo_id( const bool mix,
                    const ThermoQuantity thermoQuantity,
                    Expr::ExpressionFactory& exprFactory,
                    const Expr::TagList& thermoTags,
                    const Expr::Tag& tTag,
                    const Expr::TagList& yiTags )
{
  const int nSpec = thermoTags.size();
  std::set< Expr::ExpressionID > thermoID;

  if( mix ){
    typedef HeatCapacity_Cp<CellField>::Builder Cp;
    typedef HeatCapacity_Cv<CellField>::Builder Cv;
    typedef Enthalpy       <CellField>::Builder Enthalpy;
    typedef InternalEnergy <CellField>::Builder IntEnergy;
    Expr::ExpressionBuilder* builder = NULL;
    switch( thermoQuantity ){
      case CP  : builder = new        Cp( thermoTags[0], tTag, yiTags ); break;
      case CV  : builder = new        Cv( thermoTags[0], tTag, yiTags ); break;
      case ENTH: builder = new  Enthalpy( thermoTags[0], tTag, yiTags ); break;
      case E   : builder = new IntEnergy( thermoTags[0], tTag, yiTags ); break;
    } // switch( thermoQuantity )
    thermoID.insert( exprFactory.register_expression(builder) );
  }
  else{ // species
    typedef SpeciesHeatCapacity_Cp<CellField>::Builder SpecCp;
    typedef SpeciesHeatCapacity_Cv<CellField>::Builder SpecCv;
    typedef SpeciesEnthalpy       <CellField>::Builder SpecEnth;
    typedef SpeciesInternalEnergy <CellField>::Builder SpecE;
    for( size_t n=0; n<nSpec; ++n ){
      switch( thermoQuantity ){
        case CP  : thermoID.insert(exprFactory.register_expression( new SpecCp  (thermoTags[n], tTag, n) )); break;
        case CV  : thermoID.insert(exprFactory.register_expression( new SpecCv  (thermoTags[n], tTag, n) )); break;
        case ENTH: thermoID.insert(exprFactory.register_expression( new SpecEnth(thermoTags[n], tTag, n) )); break;
        case E   : thermoID.insert(exprFactory.register_expression( new SpecE   (thermoTags[n], tTag, n) )); break;
      } // switch(thermoQuantity)
    } // species loop
  }
  return thermoID;
}

//==============================================================================

void write_trees( const bool mix, const ThermoQuantity tq, const Expr::ExpressionTree& init, const Expr::ExpressionTree& exec)
{
  if( mix ){
    std::ofstream initOut( (thermo_name(tq) + "Mixture_init.dot").c_str() );
    init.write_tree(initOut);
    std::ofstream execOut( (thermo_name(tq) + "Mixture.dot").c_str() );
    exec.write_tree(execOut);
  }
  else{
    std::ofstream initOut( (thermo_name(tq) + "Species_init.dot").c_str() );
    init.write_tree(initOut);
    std::ofstream execOut( (thermo_name(tq) + "Species.dot").c_str() );
    exec.write_tree(execOut);
  }
}

//==============================================================================

const std::vector< std::vector<double> >
extract_mass_fracs( const Expr::TagList yiTags, Expr::FieldManagerList& fml ){
  CellField& yi0 = fml.field_ref< CellField >(yiTags[0]);
  const size_t nPts = yi0.window_with_ghost().glob_npts();
  const int nSpec = yiTags.size();

  std::vector< std::vector<double> > massFracs;
  for( size_t i=0; i<nPts; ++i ){
    std::vector<double> massFrac( nSpec, 0.0 );
    massFracs.push_back(massFrac);
  }
  for( size_t n=0; n<nSpec; ++n ){
    CellField& yi = fml.field_ref< CellField >( yiTags[n] );
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

std::vector< CellFieldPtrT >
get_cantera_results( const bool mix,
                     const bool timings,
                     const size_t canteraReps,
                     const ThermoQuantity thermoQuantity,
                     Cantera_CXX::IdealGasMix& gasMix,
                     Expr::FieldManagerList& fml,
                     const Expr::Tag& tTag,
                     const Expr::TagList& yiTags )
{
  const std::vector< std::vector<double> > massFracs = extract_mass_fracs( yiTags, fml );
  const double press = gasMix.pressure();
  const int nSpec = gasMix.nSpecies();
  const std::vector<double>& molecularWeights = gasMix.molecularWeights();

  CellField& temp   = fml.field_ref< CellField >( tTag );
# ifdef ENABLE_CUDA
  temp.set_device_as_active  ( CPU_INDEX );
# endif
  std::vector< CellFieldPtrT > canteraResults;
  if( mix ){
    canteraResults.push_back( SO::SpatialFieldStore::get<CellField>(temp) );
  }
  else{
    for( size_t n=0; n < nSpec; ++n){
      canteraResults.push_back( SO::SpatialFieldStore::get<CellField>(temp) );
    }
  }

  std::vector< std::vector<double> >::const_iterator iMass;
  CellField::const_iterator                          iTemp;
  CellField::iterator                                iCant;

  SpatialOps::Timer thermoTimer;
  if( mix ){
    const CellField::iterator iCantEnd = canteraResults[0]->end();
    thermoTimer.start();
    for( size_t rep=0; rep < canteraReps; ++rep ){
      iTemp = temp.begin();
      iMass = massFracs.begin();
      for( iCant = canteraResults[0]->begin(); iCant!=iCantEnd; ++iTemp, ++iMass, ++iCant ){
        gasMix.setMassFractions_NoNorm( &(*iMass)[0] );
        gasMix.setState_TP( *iTemp, press );
        switch(thermoQuantity){
        case CP  : *iCant=gasMix.cp_mass();        break;
        case CV  : *iCant=gasMix.cv_mass();        break;
        case ENTH: *iCant=gasMix.enthalpy_mass();  break;
        case E   : *iCant=gasMix.intEnergy_mass(); break;
        } // switch(thermoQuantity)
      }
    }
    thermoTimer.stop();
  }

  else{ // species
    std::vector<double> thermoResult(nSpec,0.0);
    const CellField::const_iterator iTempEnd = temp.end();
    thermoTimer.start();
    for( size_t rep=0; rep < canteraReps; ++rep ){
      size_t i = 0;
      iMass = massFracs.begin();
      for( iTemp = temp.begin(); iTemp != iTempEnd; ++iTemp, ++iMass, ++i){
        gasMix.setMassFractions_NoNorm( &(*iMass)[0] );
        gasMix.setState_TP( *iTemp, press );
        switch( thermoQuantity ){
        case CP  : gasMix.getPartialMolarCp(&thermoResult[0]);          break;
        case CV  : gasMix.getPartialMolarCp(&thermoResult[0]);          break;
        case ENTH: gasMix.getPartialMolarEnthalpies(&thermoResult[0]);  break;
        case E   : gasMix.getPartialMolarIntEnergies(&thermoResult[0]); break;
        }
        for( size_t n=0; n<nSpec; ++n ){
          switch( thermoQuantity ){
          case CP  : (*canteraResults[n])[i] =   thermoResult[n] / molecularWeights[n];                           break; // convert to mass basis for field comparison
          case CV  : (*canteraResults[n])[i] = ( thermoResult[n] - Cantera::GasConstant ) / molecularWeights[n];  break; // convert from molar cp to mass cv for field comparison
          case ENTH: (*canteraResults[n])[i] =   thermoResult[n] / molecularWeights[n];                           break; // convert to mass basis for field comparison
          case E   : (*canteraResults[n])[i] =   thermoResult[n] / molecularWeights[n];                           break; // convert to mass basis for field comparison
          }
        }
      }
    }
    thermoTimer.stop();
  }

  if( timings ) std::cout << "Cantera " + thermo_name(thermoQuantity) + " time " << thermoTimer.elapsed_time()/canteraReps << std::endl;
  return canteraResults;
}

//==============================================================================

bool driver( const bool timings,
             const size_t pokittReps,
             const size_t canteraReps,
             const bool mix,
             const ThermoQuantity thermoQuantity )
{
  TestHelper status( !timings );
  Cantera_CXX::IdealGasMix* const gasMix = CanteraObjects::get_gasmix();
  const int nSpec=gasMix->nSpecies();

  const Expr::Tag xTag( "XCoord", Expr::STATE_NONE );
  Expr::TagList yiTags;
  for( size_t n=0; n<nSpec; ++n ){
    yiTags.push_back( Expr::Tag( "yi_" + boost::lexical_cast<std::string>(n), Expr::STATE_NONE ) );
  }
  const Expr::Tag tTag  ( "Temperature", Expr::STATE_NONE );
  Expr::TagList thermoTags;
  if( mix )
    thermoTags.push_back(Expr::Tag( thermo_name(thermoQuantity) + " mix", Expr::STATE_NONE));
  else{
    for( size_t n=0; n<nSpec; ++n )
      thermoTags.push_back( Expr::Tag( thermo_name(thermoQuantity) + boost::lexical_cast<std::string>(n), Expr::STATE_NONE ) );
  }

  // we use an initialization tree to avoid recalculations when timing the execution
  Expr::ExpressionFactory initFactory;
  std::set<Expr::ExpressionID> initIDs;
  Expr::ExpressionID tID; // perturbed temperature
  Expr::ExpressionID yID; // mass fractions
  {
    typedef Expr::PlaceHolder     <CellField>::Builder XCoord;
    typedef       LinearMassFracs <CellField>::Builder MassFracs;
    typedef Expr::LinearFunction  <CellField>::Builder Temperature;

    initFactory.register_expression(       new XCoord       ( xTag )                  );
    yID = initFactory.register_expression( new MassFracs    ( yiTags, xTag )          );
    tID = initFactory.register_expression( new Temperature  ( tTag ,xTag, 1000, 500 ) );
    initIDs.insert( tID );
    initIDs.insert( yID );
  }

  Expr::ExpressionFactory execFactory;
  std::set<Expr::ExpressionID> execIDs;
  {
    typedef Expr::PlaceHolder      <CellField>::Builder MassFracs;
    typedef Expr::PlaceHolder      <CellField>::Builder Temperature;

    execFactory.register_expression( new MassFracs   ( yiTags ) );
    execFactory.register_expression( new Temperature ( tTag   ) );

    execIDs = register_thermo_id( mix, thermoQuantity, execFactory,
                                  thermoTags, tTag, yiTags );
  }

  Expr::ExpressionTree initTree( initIDs, initFactory, 0 );
  Expr::ExpressionTree execTree( execIDs, execFactory, 0 );

  write_trees( mix, thermoQuantity, initTree, execTree );

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
      std::cout << std::endl << thermo_name(thermoQuantity) << " test - " << nPoints << std::endl;
      execTree.execute_tree(); // sets memory high-water mark
    }

    SpatialOps::Timer timer;
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
      std::cout << "PoKiTT  " + thermo_name(thermoQuantity) + " time " << avgTime << std::endl;
    }

    const std::vector< CellFieldPtrT > canteraResults = get_cantera_results( mix,
                                                                             timings,
                                                                             canteraReps,
                                                                             thermoQuantity,
                                                                             *gasMix,
                                                                             fml,
                                                                             tTag,
                                                                             yiTags );
#   ifdef ENABLE_CUDA
    BOOST_FOREACH( Expr::Tag thermoTag, thermoTags){
      CellField& thermo = fml.field_ref< CellField >( thermoTag );
      thermo.set_device_as_active(CPU_INDEX);
    }
#   endif

    std::vector< CellFieldPtrT >::const_iterator iCantera = canteraResults.begin();
    BOOST_FOREACH( const Expr::Tag& thermoTag, thermoTags ){
      CellField& thermo = fml.field_ref< CellField >( thermoTag );
      switch( thermoQuantity ){
        case CP  :
        case CV  : {
          status( field_equal( thermo, **iCantera, 1e-14 ) || field_equal_abs( thermo, **iCantera, 1e-11 ), thermoTag.name() );
          break;
        }
        case E   :
        case ENTH: {
          status( field_equal( thermo, **iCantera, 1e-11 ) || field_equal_abs( thermo, **iCantera, 1e-8 ), thermoTag.name() );
          break;
        }
      }
      ++iCantera;
    }

    fml.deallocate_fields();
  } // number of points

  return status.ok();
}

//==============================================================================

int main( int iarg, char* carg[] )
{
  std::string inputFileName;
  std::string inpGroup;
  bool mix     = false;
  bool timings = false;
  size_t pokittReps = 1;
  size_t canteraReps = 1;

  // parse the command line options input describing the problem
  try {
    po::options_description desc("Supported Options");
    desc.add_options()
           ( "help", "print help message" )
           ( "xml-input-file", po::value<std::string>(&inputFileName)->default_value("h2o2.xml"), "Cantera xml input file name" )
           ( "phase", po::value<std::string>(&inpGroup), "name of phase in Cantera xml input file" )
           ( "mix", "Triggers mixture heat capacity test.  Otherwise, species heat capacities are tested." )
           ( "timings", "Generate comparison timings between Cantera and PoKiTT across several problem sizes" )
           ( "pokitt-reps", po::value<size_t>(&pokittReps), "Repeat the PoKiTT tests and report the average execution time")
           ( "cantera-reps", po::value<size_t>(&canteraReps), "Repeat the Cantera tests and report the average execution time");

    po::variables_map args;
    po::store( po::parse_command_line(iarg,carg,desc), args );
    po::notify(args);

    mix = args.count("mix");
    timings = args.count("timings") || args.count("pokitt-reps") || args.count("cantera-reps");

    if (args.count("help")) {
      std::cout << desc << "\n";
      return -1;
    }
  }

  catch( std::exception& err ){
    std::cout << "Error parsing input arguments\n" << err.what() << std::endl;
    return -2;
  }

  try{
    const CanteraObjects::Setup setup( "Mix", inputFileName, inpGroup );
    CanteraObjects::setup_cantera( setup );

    TestHelper status( !timings );
    status( driver(  timings, pokittReps, canteraReps, mix, CP   ), thermo_name(CP  ) );
    status( driver(  timings, pokittReps, canteraReps, mix, CV   ), thermo_name(CV  ) );
    status( driver(  timings, pokittReps, canteraReps, mix, ENTH ), thermo_name(ENTH) );
    status( driver(  timings, pokittReps, canteraReps, mix, E    ), thermo_name(E   ) );

    if( status.ok() ){
      std::cout << "\nPASS\n";
      return 0;
    }
  }
  catch( Cantera::CanteraError& ){
    Cantera::showErrors();
  }
  catch( std::exception& err ){
    std::cout << err.what() << std::endl;
  }

  std::cout << "\nFAIL\n";
  return -1;
}

//==============================================================================
