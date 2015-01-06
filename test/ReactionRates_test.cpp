/*
 * ReactionRate_test.cpp
 *
 *  Created on: Mar 31, 2014
 *      Author: Nathan Yonkee
 */

#include <iostream>
#include <stdio.h>
#include <fstream>
#include "TestHelper.h"

#include <expression/ExprLib.h>

#include <pokitt/MixtureMolWeight.h>
#include <pokitt/kinetics/ReactionRates.h>

#include <spatialops/structured/Grid.h>
#include <spatialops/structured/FieldComparisons.h>
#include <spatialops/util/TimeLogger.h>

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>

#include <cantera/IdealGasMix.h>

namespace Cantera_CXX{ class IdealGasMix; }

namespace So = SpatialOps;
typedef   So::SVolField   CellField;
typedef So::SpatFldPtr<CellField> CellFieldPtrT;

namespace po = boost::program_options;

using namespace pokitt;

//==============================================================================

const std::vector< std::vector<double> >
calculate_mass_fracs(const int nSpec, CellField& xcoord, const Expr::TagList yiTags, Expr::FieldManagerList& fml){
  Expr::FieldMgrSelector<CellField>::type& cellFM = fml.field_manager<CellField>();
  CellField& yi = cellFM.field_ref(yiTags[0]);
  CellFieldPtrT sum = So::SpatialFieldStore::get<CellField>(yi);
  *sum<<=0.0;
  for( size_t n=0; n<nSpec; ++n ){
    CellField& yi = cellFM.field_ref(yiTags[n]);
    yi <<= n + 1 + xcoord;
    *sum <<= *sum + yi;
  }
  BOOST_FOREACH( Expr::Tag yiTag, yiTags){
    CellField& yi = cellFM.field_ref(yiTag);
    yi <<= yi / *sum;
  }
  std::vector< std::vector<double> > massFracs;
  size_t nPts = xcoord.window_with_ghost().glob_npts();
  for( size_t i=0; i<nPts; ++i ){
    std::vector<double> massFrac( nSpec, 0.0 );
    massFracs.push_back(massFrac);
  }
  for( size_t n=0; n<nSpec; ++n ){
    CellField& yi = cellFM.field_ref(yiTags[n]);
#   ifdef ENABLE_CUDA
    yi.set_device_as_active( CPU_INDEX );
#   endif
    size_t i=0;
    for( CellField::iterator iY = yi.begin(); iY != yi.end(); ++iY, ++i ){
      massFracs[i][n] = *iY;
    }
#   ifdef ENABLE_CUDA
    yi.set_device_as_active( GPU_INDEX );
#   endif
  }
  return massFracs;
}

//==============================================================================

const std::vector< CellFieldPtrT >
get_cantera_results( const bool timings,
                     const size_t canteraReps,
                     Cantera_CXX::IdealGasMix& gasMix,
                     const std::vector< std::vector<double> >& massFracs,
                     CellField& temp){
  using namespace SpatialOps;

  const double refPressure=gasMix.pressure();
  const std::vector<double>& molecularWeights = gasMix.molecularWeights();
  const int nSpec = gasMix.nSpecies();

  std::vector< CellFieldPtrT > canteraResults;
  for( size_t n=0; n < nSpec; ++n){
    canteraResults.push_back( SpatialFieldStore::get<CellField>(temp) );
  }

  CellField::const_iterator                          iTemp    = temp.begin();
  const CellField::const_iterator                    iTempEnd = temp.end();
  std::vector< std::vector<double> >::const_iterator iMass = massFracs.begin();
  std::vector<double> rResult(nSpec,0.0);

  Timer cTimer;
  cTimer.start();
  for( size_t rep=0; rep < canteraReps; ++rep ){
    iMass = massFracs.begin();
    size_t i = 0;
    for( iTemp = temp.begin(); iTemp != iTempEnd; ++iTemp, ++iMass, ++i){
      gasMix.setState_TPY( *iTemp, refPressure, &(*iMass)[0]);
      gasMix.getNetProductionRates(&rResult[0]);
      for( size_t n=0; n<nSpec; ++n){
        (*canteraResults[n])[i] = rResult[n];
      }
    }
  }
  cTimer.stop();
  if( timings ) std::cout << "Cantera reaction rate time " << cTimer.elapsed_time()/canteraReps << std::endl;

  for( size_t n=0; n<nSpec; ++n){
    *canteraResults[n] <<= *canteraResults[n] * molecularWeights[n]; // convert to mass basis for field comparison
  }

  return canteraResults;
}

bool driver( const bool timings,
             const size_t pokittReps,
             const size_t canteraReps )
{
  TestHelper status( !timings );
  Cantera_CXX::IdealGasMix* const gasMix = CanteraObjects::get_gasmix();

  const int nSpec=gasMix->nSpecies();
  const double refPressure=gasMix->pressure();

  typedef Expr::PlaceHolder <CellField> Temp;
  typedef Expr::PlaceHolder <CellField> Pressure;
  typedef Expr::PlaceHolder <CellField> MassFracs;
  typedef MixtureMolWeight  <CellField> MixtureMolWeight;
  typedef ReactionRates     <CellField> ReactionRates;

  const Expr::Tag tTag  ( "Temperature", Expr::STATE_NONE);
  const Expr::Tag pTag  ( "Pressure",    Expr::STATE_NONE);
  const Expr::Tag mmwTag( "mmw",         Expr::STATE_NONE);

  Expr::TagList yiTags;
  for( size_t n=0; n<nSpec; ++n ){
    yiTags.push_back( Expr::Tag( "yi_" + boost::lexical_cast<std::string>(n), Expr::STATE_NONE ) );
  }
  Expr::TagList rTags;
  for( size_t n=0; n<nSpec; ++n )
    rTags.push_back( Expr::Tag( "ri" + boost::lexical_cast<std::string>(n), Expr::STATE_NONE ) );

  Expr::ExpressionFactory exprFactory;

  BOOST_FOREACH( Expr::Tag yiTag, yiTags){
    exprFactory.register_expression( new MassFracs::Builder( yiTag ) );
  }
  exprFactory.register_expression( new Temp             ::Builder( tTag  ) );
  exprFactory.register_expression( new Pressure         ::Builder( pTag  ) );
  exprFactory.register_expression( new MixtureMolWeight ::Builder( mmwTag, yiTags ));

  Expr::ExpressionID rRate_id = exprFactory.register_expression( new ReactionRates::Builder(rTags, tTag, pTag, yiTags, mmwTag) );

  std::vector<So::IntVec> ptvec;
  if( timings ){
    ptvec.push_back( So::IntVec(126,126,126) );
    ptvec.push_back( So::IntVec( 62, 62, 62) );
    ptvec.push_back( So::IntVec( 30, 30, 30) );
    ptvec.push_back( So::IntVec( 14, 14, 14) );
    ptvec.push_back( So::IntVec(  6,  6,  6) );
  }
  else{
    ptvec.push_back( So::IntVec(20,1,1) );
  }

  for( std::vector<So::IntVec>::iterator iPts = ptvec.begin(); iPts!= ptvec.end(); ++iPts){

    Expr::ExpressionTree tree( rRate_id, exprFactory, 0 );
    {
      std::ofstream fout( "ReactionRate.dot" );
      tree.write_tree(fout);
    }

    So::IntVec nPts = *iPts;
    const So::BoundaryCellInfo cellBCInfo = So::BoundaryCellInfo::build<CellField>(false,false,false);
    const So::GhostData cellGhosts(1);
    const So::MemoryWindow vwindow( So::get_window_with_ghost(nPts,cellGhosts,cellBCInfo) );
    CellField xcoord( vwindow, cellBCInfo, cellGhosts, NULL );

    std::vector<double> length(3,1.0);
    So::Grid grid( nPts, length );
    grid.set_coord<SpatialOps::XDIR>( xcoord );
#   ifdef ENABLE_CUDA
    xcoord.add_device( GPU_INDEX );
#   endif

    Expr::FieldManagerList fml;
    tree.register_fields( fml );
    fml.allocate_fields( Expr::FieldAllocInfo( nPts, 0, 0, false, false, false ) );
    tree.bind_fields( fml );

    using namespace SpatialOps;
    Expr::FieldMgrSelector<CellField>::type& cellFM = fml.field_manager<CellField>();

    CellField& temp = cellFM.field_ref(tTag);
    temp <<= 500.0 + 1000 * xcoord;

    CellField& p = cellFM.field_ref(pTag);
    p <<= refPressure;

    const std::vector< std::vector<double> > massFracs = calculate_mass_fracs( nSpec, xcoord, yiTags, fml );

    tree.lock_fields(fml);

    if( timings ) std::cout << std::endl << "Reaction rates test - " << vwindow.glob_npts() << std::endl;

    Timer rxnTimer;
    rxnTimer.start();
    for( size_t rep = 0; rep < pokittReps; ++rep ){
      tree.execute_tree();
    }
    rxnTimer.stop();
    if( timings ) std::cout << "PoKiTT  reaction rate time " << rxnTimer.elapsed_time()/pokittReps << std::endl;

#   ifdef ENABLE_CUDA
    BOOST_FOREACH( Expr::Tag rTag, rTags){
      CellField& r = fml.field_manager<CellField>().field_ref(rTag);
      r.add_device(CPU_INDEX);
    }
    temp.set_device_as_active( CPU_INDEX );
#   endif

    const std::vector< CellFieldPtrT > canteraResults = get_cantera_results( timings,
                                                                             canteraReps,
                                                                             *gasMix,
                                                                             massFracs,
                                                                             temp );

    std::vector< CellFieldPtrT >::const_iterator iCantera = canteraResults.begin();
    BOOST_FOREACH( const Expr::Tag& rTag, rTags ){
      CellField& r = cellFM.field_ref(rTag);
      status( field_equal( r, **iCantera, 1e-8 ) || field_equal_abs( r, **iCantera, 1e-10 ), rTag.name() );
      ++iCantera;
    }

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
