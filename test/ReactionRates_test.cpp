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

#include <test/TemperaturePowers.h>
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

void calculate_mass_fracs(const int nSpec, CellField& xcoord, const Expr::TagList yiTags, Expr::FieldManagerList& fml){
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
}

//==============================================================================

const std::vector< std::vector<double> >
mass_fracs(const int nPts, const int nSpec){
  std::vector< std::vector<double> > massFracs;
  double sum;
  for( size_t i=0; i < nPts+2; ++i){
    std::vector<double> massFrac;
    sum = 0.0;
    for( size_t n=0; n < nSpec; ++n){
      massFrac.push_back(1 + n + (i-0.5)/ nPts);
      sum += massFrac[n];
    }
    for( size_t n=0; n < nSpec; ++n)
      massFrac[n] = massFrac[n]/sum;
    massFracs.push_back(massFrac);
  }
  return massFracs;
}

//==============================================================================

const std::vector< CellFieldPtrT >
get_cantera_results( const bool timings,
                     const size_t repeats,
                     Cantera_CXX::IdealGasMix& gasMix,
                     const int nPts,
                     CellField& temp){
  using namespace SpatialOps;

  const double refPressure=gasMix.pressure();
  const std::vector<double>& molecularWeights = gasMix.molecularWeights();
  const int nSpec = gasMix.nSpecies();

  const std::vector< std::vector<double> > massFracs = mass_fracs( nPts, nSpec);

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
  for( size_t rep=0; rep < repeats; ++rep ){
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
  if( timings ) std::cout << "Cantera reaction rate time " << cTimer.elapsed_time()/repeats << std::endl;

  for( size_t n=0; n<nSpec; ++n){
    *canteraResults[n] <<= *canteraResults[n] * molecularWeights[n]; // convert to mass basis for field comparison
  }

  return canteraResults;
}

bool driver( const bool timings, const size_t repeats )
{
  TestHelper status( !timings );
  Cantera_CXX::IdealGasMix* const gasMix = CanteraObjects::get_gasmix();

  const int nSpec=gasMix->nSpecies();
  const double refPressure=gasMix->pressure();

  typedef Expr::PlaceHolder <CellField> Temp;
  typedef TemperaturePowers <CellField> TemperaturePowers;
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
  exprFactory.register_expression( new TemperaturePowers::Builder( tTag  ) );
  exprFactory.register_expression( new Pressure         ::Builder( pTag  ) );
  exprFactory.register_expression( new MixtureMolWeight ::Builder( mmwTag, yiTags ));

  Expr::ExpressionID rRate_id = exprFactory.register_expression( new ReactionRates::Builder(rTags, tTag, pTag, yiTags, mmwTag) );

  std::vector<int> ptvec;
  if( timings ){
    ptvec.push_back(8*8*8);
    ptvec.push_back(16*16*16);
    ptvec.push_back(32*32*32);
    ptvec.push_back(64*64*64);
    ptvec.push_back(128*128*128);
  }
  else{
    ptvec.push_back(10);
  }

  for( std::vector<int>::iterator iPts = ptvec.begin(); iPts!= ptvec.end(); ++iPts){

    Expr::ExpressionTree tree( rRate_id, exprFactory, 0 );
    {
      std::ofstream fout( "ReactionRate.dot" );
      tree.write_tree(fout);
    }

    So::IntVec npts(*iPts,1,1);
    const So::BoundaryCellInfo cellBCInfo = So::BoundaryCellInfo::build<CellField>(false,false,false);
    const So::GhostData cellGhosts(1);
    const So::MemoryWindow vwindow( So::get_window_with_ghost(npts,cellGhosts,cellBCInfo) );
    CellField xcoord( vwindow, cellBCInfo, cellGhosts, NULL );

    std::vector<double> length(3,1.0);
    So::Grid grid( npts, length );
    grid.set_coord<SpatialOps::XDIR>( xcoord );
#   ifdef ENABLE_CUDA
    xcoord.add_device( GPU_INDEX );
#   endif

    Expr::FieldManagerList fml;
    tree.register_fields( fml );
    fml.allocate_fields( Expr::FieldAllocInfo( npts, 0, 0, false, false, false ) );
    tree.bind_fields( fml );

    using namespace SpatialOps;
    Expr::FieldMgrSelector<CellField>::type& cellFM = fml.field_manager<CellField>();

    CellField& temp = cellFM.field_ref(tTag);
    temp <<= 500.0 + 1000 * xcoord;

    CellField& p = cellFM.field_ref(pTag);
    p <<= refPressure;

    calculate_mass_fracs( nSpec, xcoord, yiTags, fml );

    tree.lock_fields(fml);

    if( timings ) std::cout << std::endl << "Reaction rates test - " << *iPts << std::endl;

    Timer rxnTimer;
    rxnTimer.start();
    for( size_t rep = 0; rep < repeats; ++rep ){
      tree.execute_tree();
    }
    rxnTimer.stop();
    if( timings ) std::cout << "PoKiTT  reaction rate time " << rxnTimer.elapsed_time()/repeats << std::endl;

#   ifdef ENABLE_CUDA
    BOOST_FOREACH( Expr::Tag rTag, rTags){
      CellField& r = fml.field_manager<CellField>().field_ref(rTag);
      r.add_device(CPU_INDEX);
    }
    temp.set_device_as_active( CPU_INDEX );
#   endif

    const std::vector< CellFieldPtrT > canteraResults = get_cantera_results( timings,
                                                                             repeats,
                                                                             *gasMix,
                                                                             *iPts,
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
  size_t repeats = 1;

  // parse the command line options input describing the problem
  try {
    po::options_description desc("Supported Options");
    desc.add_options()
           ( "help", "print help message" )
           ( "xml-input-file", po::value<std::string>(&inputFileName), "Cantera xml input file name" )
           ( "phase", po::value<std::string>(&inpGroup), "name of phase in Cantera xml input file" )
           ( "timings", "Generate comparison timings between Cantera and PoKiTT across several problem sizes" )
           ( "repeats", po::value<size_t>(&repeats), "Repeat the tests and report the average execution time");

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
    status( driver( timings, repeats ), "Reaction Rates" );

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
