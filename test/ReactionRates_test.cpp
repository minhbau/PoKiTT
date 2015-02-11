/*
 * ReactionRate_test.cpp
 *
 *  Created on: Mar 31, 2014
 *      Author: Nathan Yonkee
 */

#include <iostream>
#include "TestHelper.h"
#include "LinearMassFracs.h"
#include <pokitt/MixtureMolWeight.h>
#include <pokitt/kinetics/ReactionRates.h>

#include <expression/ExprLib.h>

#include <spatialops/structured/Grid.h>
#include <spatialops/structured/FieldComparisons.h>
#include <spatialops/util/TimeLogger.h>

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>

#include <cantera/IdealGasMix.h>

namespace Cantera_CXX{ class IdealGasMix; }

namespace SO = SpatialOps;
typedef SO::SVolField  CellField;
typedef SO::SpatFldPtr<CellField> CellFieldPtrT;

namespace po = boost::program_options;

using namespace pokitt;

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

const std::vector< CellFieldPtrT >
get_cantera_results( const bool timings,
                     const size_t canteraReps,
                     Cantera_CXX::IdealGasMix& gasMix,
                     Expr::FieldManagerList& fml,
                     const Expr::Tag& tTag,
                     const Expr::TagList& yiTags,
                     const Expr::Tag& pTag){

  const std::vector<double>& molecularWeights = gasMix.molecularWeights();
  const int nSpec = gasMix.nSpecies();

  const std::vector< std::vector<double> > massFracs = extract_mass_fracs( yiTags, fml );

  CellField& temp   = fml.field_ref< CellField >( tTag );
  CellField& press  = fml.field_ref< CellField >( pTag );
# ifdef ENABLE_CUDA
  temp.set_device_as_active  ( CPU_INDEX );
  press.set_device_as_active ( CPU_INDEX );
# endif

  std::vector< CellFieldPtrT > canteraResults;
  for( size_t n=0; n < nSpec; ++n){
    canteraResults.push_back( SO::SpatialFieldStore::get<CellField>(temp) );
  }

  std::vector< std::vector<double> >::const_iterator iMass;
  CellField::const_iterator                          iTemp;
  CellField::const_iterator                          iPress;
  CellField::iterator                                iCant;

  std::vector<double> rResult(nSpec,0.0);
  Timer cTimer;
  cTimer.start();
  for( size_t rep=0; rep < canteraReps; ++rep ){
    iPress = press.begin();
    iMass  = massFracs.begin();
    size_t i = 0;
    for( iTemp = temp.begin(); iTemp != temp.end(); ++iTemp, ++iMass, ++i){
      gasMix.setMassFractions_NoNorm( &(*iMass)[0] );
      gasMix.setState_TP( *iTemp, *iPress );
      gasMix.getNetProductionRates(&rResult[0]);
      for( size_t n=0; n<nSpec; ++n){
        (*canteraResults[n])[i] = rResult[n] * molecularWeights[n];
      }
    }
  }
  cTimer.stop();
  if( timings ) std::cout << "Cantera reaction rate time " << cTimer.elapsed_time()/canteraReps << std::endl;

  return canteraResults;
}

bool driver( const bool timings,
             const size_t pokittReps,
             const size_t canteraReps )
{
  TestHelper status( !timings );
  Cantera_CXX::IdealGasMix* const gasMix = CanteraObjects::get_gasmix();
  const int nSpec=gasMix->nSpecies();

  const Expr::Tag xTag  ( "XCoord",      Expr::STATE_NONE );
  const Expr::Tag tTag  ( "Temperature", Expr::STATE_NONE);
  const Expr::Tag pTag  ( "Pressure",    Expr::STATE_NONE);
  const Expr::Tag mmwTag( "mmw",         Expr::STATE_NONE);
  Expr::TagList yiTags;
  Expr::TagList rTags;
  for( size_t n=0; n<nSpec; ++n ){
    yiTags.push_back( Expr::Tag( "yi_" + boost::lexical_cast<std::string>(n), Expr::STATE_NONE ) );
    rTags.push_back(  Expr::Tag( "ri" + boost::lexical_cast<std::string>(n), Expr::STATE_NONE ) );
  }

  // we use an initialization tree to avoid recalculations when timing the execution
  Expr::ExpressionFactory initFactory;
  std::set<Expr::ExpressionID> initIDs;
  Expr::ExpressionID tID; // temperature
  Expr::ExpressionID pID; // pressure
  Expr::ExpressionID mID; // mixture molecular weight
  Expr::ExpressionID yID; // mass fractions
  {
    typedef Expr::PlaceHolder     <CellField>::Builder XCoord;
    typedef       LinearMassFracs <CellField>::Builder MassFracs;
    typedef Expr::ConstantExpr    <CellField>::Builder Pressure;
    typedef       MixtureMolWeight<CellField>::Builder MixMolWeight;
    typedef Expr::LinearFunction  <CellField>::Builder Temperature;

    initFactory.register_expression(       new XCoord       ( xTag )                     );
    yID = initFactory.register_expression( new MassFracs    ( yiTags, xTag )             );
    pID = initFactory.register_expression( new Pressure     ( pTag, gasMix->pressure() ) );
    mID = initFactory.register_expression( new MixMolWeight ( mmwTag, yiTags )           );
    tID = initFactory.register_expression( new Temperature  ( tTag ,xTag, 1000, 500 )    );
    initIDs.insert( tID );
    initIDs.insert( pID );
    initIDs.insert( mID );
    initIDs.insert( yID );
  }

  Expr::ExpressionFactory execFactory;
  Expr::ExpressionID execID;
  {
    typedef Expr::PlaceHolder   <CellField>::Builder MassFracs;
    typedef Expr::PlaceHolder   <CellField>::Builder Temperature;
    typedef Expr::PlaceHolder   <CellField>::Builder Pressure;
    typedef Expr::PlaceHolder   <CellField>::Builder MixMolWeight;
    typedef       ReactionRates <CellField>::Builder ReactionRates;

    execFactory.register_expression( new MassFracs   ( yiTags ) );
    execFactory.register_expression( new Temperature ( tTag   ) );
    execFactory.register_expression( new Pressure    ( pTag   ) );
    execFactory.register_expression( new Pressure    ( mmwTag ) );
    execID = execFactory.register_expression( new ReactionRates(rTags, tTag, pTag, yiTags, mmwTag) );
  }

  Expr::ExpressionTree initTree( initIDs, initFactory, 0 );
  Expr::ExpressionTree execTree( execID , execFactory, 0 );
  {
    std::ofstream fout( "ReactionRate_init.dot" );
    initTree.write_tree(fout);
  }
  {
    std::ofstream fout( "ReactionRate.dot" );
    execTree.write_tree(fout);
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
    sizeVec.push_back( SO::IntVec(126,126,126) );
    sizeVec.push_back( SO::IntVec( 62, 62, 62) );
    sizeVec.push_back( SO::IntVec( 30, 30, 30) );
    sizeVec.push_back( SO::IntVec( 14, 14, 14) );
    sizeVec.push_back( SO::IntVec(  6,  6,  6) );
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

    if( timings ) std::cout << std::endl << "Reaction rates test - " << nPoints << std::endl;

    initTree.execute_tree();
    Timer rxnTimer;
    rxnTimer.start();
    for( size_t rep = 0; rep < pokittReps; ++rep ){
      execTree.execute_tree();
    }
    rxnTimer.stop();
    if( timings ) std::cout << "PoKiTT  reaction rate time " << rxnTimer.elapsed_time()/pokittReps << std::endl;

    const std::vector< CellFieldPtrT > canteraResults = get_cantera_results( timings,
                                                                             canteraReps,
                                                                             *gasMix,
                                                                             fml,
                                                                             tTag,
                                                                             yiTags,
                                                                             pTag );
#   ifdef ENABLE_CUDA
    BOOST_FOREACH( Expr::Tag rTag, rTags){
      CellField& r = fml.field_ref< CellField >( rTag );
      r.add_device(CPU_INDEX);
    }
#   endif

    std::vector< CellFieldPtrT >::const_iterator iCantera = canteraResults.begin();
    BOOST_FOREACH( const Expr::Tag& rTag, rTags ){
      CellField& r = fml.field_ref< CellField >( rTag );
      status( field_equal( r, **iCantera, 1e-8 ) || field_equal_abs( r, **iCantera, 1e-10 ), rTag.name() );
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
  bool timings = false;
  size_t pokittReps = 1;
  size_t canteraReps = 1;

  // parse the command line options input describing the problem
  try {
    po::options_description desc("Supported Options");
    desc.add_options()
           ( "help", "print help message" )
           ( "xml-input-file", po::value<std::string>(&inputFileName), "Cantera xml input file name" )
           ( "phase", po::value<std::string>(&inpGroup), "name of phase in Cantera xml input file" )
           ( "timings", "Generate comparison timings between Cantera and PoKiTT across several problem sizes" )
           ( "pokitt-reps", po::value<size_t>(&pokittReps), "Repeat the PoKiTT tests and report the average execution time")
           ( "cantera-reps", po::value<size_t>(&canteraReps), "Repeat the Cantera tests and report the average execution time");

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
    std::cout << "Error parsing input arguments\n" << err.what() << std::endl;
    return -2;
  }

  try{
    const CanteraObjects::Setup setup( "Mix", inputFileName, inpGroup );
    CanteraObjects::setup_cantera( setup );

    TestHelper status( !timings );
    status( driver( timings, pokittReps, canteraReps ), "Reaction Rates" );

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
