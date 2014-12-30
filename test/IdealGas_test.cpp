/*
 * IdealGas_test.cpp
 *
 *  Created on: October 15, 2014
 *      Author: Nathan Yonkee
 */

#include <iostream>
#include <stdio.h>
#include <fstream>
#include "TestHelper.h"

#include <pokitt/thermo/Density.h>
#include <pokitt/thermo/Pressure.h>
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
typedef So::SVolField   CellField;
typedef So::SpatFldPtr<CellField> CellFieldPtrT;

namespace Cantera_CXX{ class IdealGasMix; }

namespace po = boost::program_options;

using namespace pokitt;

//==============================================================================

enum GasQuantity{
  P,
  RHO
};

std::string property_name( const GasQuantity q )
{
  switch(q){
    case P  : return "pressure";
    case RHO: return "density";
  }
}

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

CellFieldPtrT
get_cantera_result( const bool timings,
                    const size_t canteraReps,
                    const GasQuantity gasQuantity,
                    Cantera_CXX::IdealGasMix& gasMix,
                    const std::vector< std::vector<double> >& massFracs,
                    CellField& temp,
                    CellField& refQuantity)
{
# ifdef ENABLE_CUDA
  temp.set_device_as_active       ( CPU_INDEX );
  refQuantity.set_device_as_active( CPU_INDEX );
# endif
  const int nSpec=gasMix.nSpecies();

  CellFieldPtrT canteraResult = So::SpatialFieldStore::get<CellField>(temp, CPU_INDEX);

  std::vector< std::vector<double> >::const_iterator iMass = massFracs.begin();
  CellField::const_iterator iTemp = temp.begin();
  CellField::const_iterator iRef  = refQuantity.begin();
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

  typedef Expr::PlaceHolder <CellField> Temperature;
  typedef Expr::PlaceHolder <CellField> RefQuantity;
  typedef Expr::PlaceHolder <CellField> MassFracs;
  typedef MixtureMolWeight  <CellField> MixtureMolWeight;

  const Expr::Tag tTag(   "Temperature",        Expr::STATE_NONE );
  Expr::Tag refTag;
  switch( gasQuantity ){
  case P:
    refTag = Expr::Tag( "Density",  Expr::STATE_NONE ); break;
  case RHO:
    refTag = Expr::Tag( "Pressure", Expr::STATE_NONE ); break;
  }
  Expr::TagList yiTags;
  for( size_t n=0; n<nSpec; ++n ){
    yiTags.push_back( Expr::Tag( "yi_" + boost::lexical_cast<std::string>(n), Expr::STATE_NONE ) );
  }
  const Expr::Tag mmwTag ("Mixture Molecular Weight", Expr::STATE_NONE);
  const Expr::Tag gasTag (property_name(gasQuantity), Expr::STATE_NONE);

  Expr::ExpressionFactory exprFactory;

  BOOST_FOREACH( const Expr::Tag& yiTag, yiTags){
    exprFactory.register_expression( new MassFracs::Builder( yiTag ) );
  }
  exprFactory.register_expression( new Temperature::Builder     ( tTag           ) );
  exprFactory.register_expression( new RefQuantity::Builder     ( refTag         ) );
  exprFactory.register_expression( new MixtureMolWeight::Builder( mmwTag, yiTags ) );

  Expr::ExpressionID gas_id;
  switch( gasQuantity ){
  case P:
    gas_id = exprFactory.register_expression( new Pressure<CellField>::Builder(gasTag, tTag, refTag, mmwTag) );
    break;
  case RHO:
    gas_id = exprFactory.register_expression( new Density<CellField>::Builder (gasTag, tTag, refTag, mmwTag) );
    break;
  }

  std::vector<So::IntVec> ptvec;
  if( timings ){
    ptvec.push_back( So::IntVec(  6,  6,  6) );
    ptvec.push_back( So::IntVec( 14, 14, 14) );
    ptvec.push_back( So::IntVec( 30, 30, 30) );
    ptvec.push_back( So::IntVec( 62, 62, 62) );
    ptvec.push_back( So::IntVec(126,126,126) );
  }
  else{
    ptvec.push_back( So::IntVec(20,1,1) );
  }

  for( std::vector<So::IntVec>::iterator iPts = ptvec.begin(); iPts!= ptvec.end(); ++iPts){

    Expr::ExpressionTree gasTree( gas_id , exprFactory, 0 );
    {
      std::ofstream gasGraph( "IdealGas.dot" );
      gasTree.write_tree( gasGraph );
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

    gasTree.register_fields( fml );
    fml.allocate_fields( Expr::FieldAllocInfo( nPts, 0, 0, false, false, false ) );
    gasTree.bind_fields( fml );

    using namespace SpatialOps;
    Expr::FieldMgrSelector<CellField>::type& cellFM = fml.field_manager< CellField>();

    CellField& temp        = cellFM.field_ref( tTag   );
    CellField& refQuantity = cellFM.field_ref( refTag );
    CellField& mixMW       = cellFM.field_ref( mmwTag );
    CellField& gasField    = cellFM.field_ref( gasTag );

    temp <<= 500.0 + 1000 * xcoord;

    switch( gasQuantity ){
      case P:   refQuantity <<= gasMix->density();  break;
      case RHO: refQuantity <<= gasMix->pressure(); break;
    }

    const std::vector< std::vector<double> > massFracs = calculate_mass_fracs( nSpec, xcoord, yiTags, fml );

    gasTree.lock_fields( fml );  // prevent fields from being deallocated so that we can get them after graph execution.

    if( timings ) std::cout << std::endl << property_name(gasQuantity) + " test - " << vwindow.glob_npts() << std::endl;

    Timer gasTimer;
    gasTimer.start();
    for( size_t rep = 0; rep < pokittReps; ++rep ){
      gasTree.execute_tree();
    }
    gasTimer.stop();

    if( timings ) std::cout << "PoKiTT  " + property_name(gasQuantity) + " time " << gasTimer.elapsed_time()/pokittReps << std::endl;

#   ifdef ENABLE_CUDA
    gasField.set_device_as_active( CPU_INDEX );
#   endif

    CellFieldPtrT canteraResult = get_cantera_result( timings, canteraReps, gasQuantity, *gasMix, massFracs, temp, refQuantity);

    status( field_equal(gasField, *canteraResult, 1e-12), gasTag.name() );

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
