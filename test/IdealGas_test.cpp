/*
 * IdealGas_test.cpp
 *
 *  Created on: October 15, 2014
 *      Author: Nathan Yonkee
 */

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

#include <cantera/kernel/ct_defs.h> // contains value of Cantera::GasConstant
#include <cantera/IdealGasMix.h>

namespace So = SpatialOps;
typedef So::SVolField  CellField;
typedef So::SpatFldPtr<CellField> CellFieldPtrT;

namespace Cantera_CXX{ class IdealGasMix; }

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
                    Cantera_CXX::IdealGasMix& gasMix,
                    Expr::FieldManagerList& fml,
                    const Expr::TagList& yiTags,
                    const Expr::Tag& tTag,
                    const Expr::Tag& refTag )
{
  CellField& temp = fml.field_ref< CellField >( tTag );
  CellField& refQuantity = fml.field_ref< CellField >( refTag );
  const std::vector< std::vector<double> > massFracs = extract_mass_fracs( yiTags, fml );
# ifdef ENABLE_CUDA
  temp.set_device_as_active       ( CPU_INDEX );
  refQuantity.set_device_as_active( CPU_INDEX );
# endif
  const int nSpec=gasMix.nSpecies();

  CellFieldPtrT canteraResult = So::SpatialFieldStore::get<CellField>(temp, CPU_INDEX);

  std::vector< std::vector<double> >::const_iterator iMass = massFracs.begin();
  CellField::const_iterator iTemp;
  CellField::const_iterator iRef;
  Timer gasTime;
  gasTime.start();
  for( size_t rep=0; rep < canteraReps; ++rep ){
    iTemp = temp.begin();
    iMass = massFracs.begin();
    iRef  = refQuantity.begin();
    for(CellField::iterator iCant = canteraResult->begin(); iCant!=canteraResult->end(); ++iTemp, ++iRef, ++iMass, ++iCant){
      switch( gasQuantity ){
      case P:
        gasMix.setState_TRY( *iTemp, *iRef, &(*iMass)[0]);
        *iCant=gasMix.pressure();
        break;
      case RHO:
        gasMix.setState_TPY( *iTemp, *iRef, &(*iMass)[0]);
        *iCant=gasMix.density();
        break;
      case NU:
        gasMix.setState_TPY( *iTemp, *iRef, &(*iMass)[0]);
        *iCant=1/gasMix.density();
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
  Cantera_CXX::IdealGasMix* const gasMix = CanteraObjects::get_gasmix();

  const int nSpec = gasMix->nSpecies();

  typedef Expr::PlaceHolder     <CellField>::Builder XCoord;
  typedef Expr::LinearFunction  <CellField>::Builder Temperature;
  typedef Expr::ConstantExpr    <CellField>::Builder RefQuantity;
  typedef       LinearMassFracs <CellField>::Builder MassFracs;
  typedef       MixtureMolWeight<CellField>::Builder MixtureMolWeight;
  typedef       Pressure        <CellField>::Builder Pressure;
  typedef       Density         <CellField>::Builder Density;
  typedef       SpecificVol     <CellField>::Builder SpecVol;

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
  case NU:  refTag = Expr::Tag( "SpecVol",  Expr::STATE_NONE ); break;
  }
  const Expr::Tag gasTag (property_name(gasQuantity), Expr::STATE_NONE);

  Expr::ExpressionFactory exprFactory;

  exprFactory.register_expression( new XCoord          ( xTag           ) );
  exprFactory.register_expression( new MassFracs       ( yiTags, xTag   ) );
  exprFactory.register_expression( new Temperature     ( tTag,   xTag, 1000, 500 ) );
  exprFactory.register_expression( new MixtureMolWeight( mmwTag, yiTags ) );

  Expr::ExpressionID gas_id;
  switch( gasQuantity ){
  case P:

    gas_id = exprFactory.register_expression( new Pressure   ( gasTag, tTag, refTag, mmwTag) );
    exprFactory.register_expression(          new RefQuantity( refTag, gasMix->density()   ) );
    break;
  case RHO:
    gas_id = exprFactory.register_expression( new Density    ( gasTag, tTag, refTag, mmwTag) );
    exprFactory.register_expression(          new RefQuantity( refTag, gasMix->pressure() ) );
    break;
  case NU:
    gas_id = exprFactory.register_expression( new SpecVol    ( gasTag, tTag, refTag, mmwTag) );
    exprFactory.register_expression(          new RefQuantity( refTag, gasMix->pressure() ) );
    break;
  }

  Expr::ExpressionTree gasTree( gas_id , exprFactory, 0 );
  {
    std::ofstream gasGraph( "IdealGas.dot" );
    gasTree.write_tree( gasGraph );
  }
  Expr::FieldManagerList fml;
  gasTree.register_fields( fml );
  gasTree.bind_fields( fml );
  gasTree.lock_fields( fml );

  std::vector<So::IntVec> sizeVec;
  if( timings ){
    sizeVec.push_back( So::IntVec(126,126,126) );
    sizeVec.push_back( So::IntVec( 62, 62, 62) );
    sizeVec.push_back( So::IntVec( 30, 30, 30) );
    sizeVec.push_back( So::IntVec( 14, 14, 14) );
    sizeVec.push_back( So::IntVec(  6,  6,  6) );
  }
  else{
    sizeVec.push_back( So::IntVec( 20,  1,  1) );
  }

  for( std::vector<So::IntVec>::iterator iSize = sizeVec.begin(); iSize!= sizeVec.end(); ++iSize){

    So::IntVec gridSize = *iSize;
    fml.allocate_fields( Expr::FieldAllocInfo( gridSize, 0, 0, false, false, false ) );
    So::Grid grid( gridSize, So::DoubleVec(1,1,1) );

    CellField& xcoord = fml.field_ref< CellField >( xTag );
    grid.set_coord<So::XDIR>( xcoord );
    const int nPoints = xcoord.window_with_ghost().glob_npts();
#   ifdef ENABLE_CUDA
    xcoord.set_device_as_active( GPU_INDEX );
#   endif

    if( timings ) std::cout << std::endl << property_name(gasQuantity) + " test - " << nPoints << std::endl;

    Timer gasTimer;
    gasTimer.start();
    for( size_t rep = 0; rep < pokittReps; ++rep ){
      gasTree.execute_tree();
    }
    gasTimer.stop();

    if( timings ) std::cout << "PoKiTT  " + property_name(gasQuantity) + " time " << gasTimer.elapsed_time()/pokittReps << std::endl;

    CellFieldPtrT canteraResult = get_cantera_result( timings, canteraReps, gasQuantity, *gasMix, fml, yiTags, tTag, refTag );
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
               ( "xml-input-file", po::value<std::string>(&inputFileName), "Cantera xml input file name" )
               ( "phase", po::value<std::string>(&inpGroup), "name of phase in Cantera xml input file" )
               ( "timings", "Generate comparison timings between Cantera and PoKiTT across several problem sizes" )
               ( "pokitt-reps", po::value<size_t>(&pokittReps), "Repeat the PoKiTT tests and report the average execution time")
               ( "cantera-reps", po::value<size_t>(&canteraReps), "Repeat the Cantera tests and report the average execution time");

  try {
    po::variables_map args;
    po::store( po::parse_command_line(iarg,carg,desc), args );
    po::notify(args);

    timings = args.count("timings") > 0;

    if (!args.count("xml-input-file")){
      std::cout << "You must enter an xml input file for Cantera" << std::endl;
      std::cout << desc << std::endl;
      return 1;
    }

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
