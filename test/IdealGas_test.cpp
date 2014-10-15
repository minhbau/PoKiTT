/*
 * Density_test.cpp
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

#include <boost/timer.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>

#include <cantera/kernel/ct_defs.h> // contains value of Cantera::GasConstant
#include <cantera/IdealGasMix.h>

namespace So = SpatialOps;
typedef So::SVolField   CellField;

namespace Cantera_CXX{ class IdealGasMix; }

namespace po = boost::program_options;

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

const std::vector< std::vector<double> > mass_fracs(const int nPts, const int nSpec){
  std::vector< std::vector<double> > massfracs;
  for( size_t i=0; i<nPts+2; ++i){
    std::vector<double> massfrac;
    double sum = 0.0;
    for( size_t n=0; n<nSpec; ++n){
      massfrac.push_back(1 + n + (i-0.5)/ nPts);
      sum+=massfrac[n];
    }
    for( size_t n=0; n<nSpec; ++n)
      massfrac[n] = massfrac[n]/sum;
    massfracs.push_back(massfrac);
  }
  return massfracs;
}

//==============================================================================

So::SpatFldPtr<CellField>
get_cantera_result( const bool timings,
                    const GasQuantity gasQuantity,
                    Cantera_CXX::IdealGasMix& gasMix,
                    const int nPts,
                    CellField& temp,
                    CellField& xcoord)
{
# ifdef ENABLE_CUDA
  temp.set_device_as_active  ( CPU_INDEX );
  xcoord.set_device_as_active( CPU_INDEX );
# endif
  const int nSpec=gasMix.nSpecies();
  double refQuantity;
  switch( gasQuantity ){
  case P: refQuantity = gasMix.density(); break;
  case RHO: refQuantity = gasMix.pressure(); break;
  }

  const std::vector< std::vector<double> > massFracs = mass_fracs( nPts, nSpec);
  So::SpatFldPtr<CellField> canteraResult = So::SpatialFieldStore::get<CellField>(temp, CPU_INDEX);

  std::vector< std::vector<double> >::const_iterator iMass = massFracs.begin();
  CellField::const_iterator iTemp = temp.begin();
  boost::timer gasTime;
  for(CellField::iterator iCant = canteraResult->begin(); iCant!=canteraResult->end(); ++iTemp, ++iMass, ++iCant){
    switch( gasQuantity ){
    case P:
      gasMix.setState_TRY( *iTemp, refQuantity, &(*iMass)[0]);
      *iCant=gasMix.pressure();
      break;
    case RHO:
      gasMix.setState_TPY( *iTemp, refQuantity, &(*iMass)[0]);
      *iCant=gasMix.density();
      break;
    }
  }
  if( timings ) std::cout << "Cantera " + property_name(gasQuantity) + " time " << gasTime.elapsed() << std::endl;
  return canteraResult;
}

//==============================================================================

bool driver( bool timings,
             const GasQuantity gasQuantity )
{
  TestHelper status( !timings );
  Cantera_CXX::IdealGasMix* const gasMix = CanteraObjects::get_gasmix();

  const int nSpec = gasMix->nSpecies();
  const std::vector<double>& molecularWeights = gasMix->molecularWeights();
  size_t n;

  double refQuantity;
  switch( gasQuantity ){
  case P: refQuantity = gasMix->density(); break;
  case RHO: refQuantity = gasMix->pressure(); break;
  }

  typedef Expr::PlaceHolder <CellField> Temperature;
  typedef Expr::PlaceHolder <CellField> RefQuantity;
  typedef Expr::PlaceHolder <CellField> MassFracs;
  typedef MixtureMolWeight  <CellField> MixtureMolWeight;

  const Expr::Tag tTag ( "Temperature", Expr::STATE_NONE);
  const Expr::Tag refTag  ( "Reference Quantity", Expr::STATE_NONE);
  const Expr::Tag yiTag ( "yi",          Expr::STATE_NONE );
  Expr::TagList yiTags;
  for( n=0; n<nSpec; ++n ){
    std::ostringstream name;
    name << yiTag.name() << "_" << n;
    yiTags.push_back( Expr::Tag(name.str(),yiTag.context()) );
  }
  const Expr::Tag mmwTag ("Mixture Molecular Weight", Expr::STATE_NONE);
  const Expr::Tag gasTag (property_name(gasQuantity), Expr::STATE_NONE);

  Expr::ExpressionFactory exprFactory;

  exprFactory.register_expression( new Temperature::Builder (tTag) );
  exprFactory.register_expression( new RefQuantity::Builder (refTag) );
  BOOST_FOREACH( Expr::Tag yiTag, yiTags){
    exprFactory.register_expression( new MassFracs      ::Builder( yiTag ) );
  }
  exprFactory.register_expression( new MixtureMolWeight ::Builder( mmwTag, yiTag, molecularWeights));

  Expr::ExpressionID gas_id;
  switch( gasQuantity ){
  case P:
    gas_id = exprFactory.register_expression( new Pressure<CellField>::Builder (gasTag, tTag, refTag, mmwTag) );
    break;
  case RHO:
    gas_id = exprFactory.register_expression( new Density<CellField>::Builder (gasTag, tTag, refTag, mmwTag) );
    break;
  }

  std::vector<int> ptvec;
  if( timings ){
    ptvec.push_back(8*8*8);
    ptvec.push_back(16*16*16);
    ptvec.push_back(32*32*32);
    ptvec.push_back(64*64*64);
#   ifdef ENABLE_CUDA
    ptvec.push_back(128*128*128);
#   endif
  }
  else ptvec.push_back(10);

  for( std::vector<int>::iterator iPts = ptvec.begin(); iPts!= ptvec.end(); ++iPts){

    Expr::ExpressionTree gasTree( gas_id , exprFactory, 0 );
    {
      std::ofstream gasGraph( "IdealGas.dot" );
      gasTree.write_tree( gasGraph );
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

    gasTree.register_fields( fml );
    fml.allocate_fields( Expr::FieldAllocInfo( npts, 0, 0, false, false, false ) );
    gasTree.bind_fields( fml );

    using namespace SpatialOps;
    Expr::FieldMgrSelector<CellField>::type& cellFM = fml.field_manager< CellField>();

    CellField& temp        = cellFM.field_ref( tTag );
    CellField& refField = cellFM.field_ref( refTag );
    CellField& mixMW       = cellFM.field_ref( mmwTag );
    CellField& gasField = cellFM.field_ref( gasTag );

    temp <<= 500.0 + 1000 * xcoord;
    refField <<= refQuantity;

    SpatFldPtr<CellField> sum  = SpatialFieldStore::get<CellField>(temp);
    *sum<<=0.0;
    for( n=0; n<nSpec; ++n ){ //normalize the mass fractions
      CellField& yi = fml.field_manager<CellField>().field_ref(yiTags[n]);
      yi <<= n + 1 + xcoord;
      *sum <<= *sum + yi;
    }
    BOOST_FOREACH( Expr::Tag yiTag, yiTags){
      CellField& yi = cellFM.field_ref(yiTag);
      yi <<= yi / *sum;
    }

    gasTree.lock_fields( fml );  // prevent fields from being deallocated so that we can get them after graph execution.

    if( timings ) std::cout << std::endl << property_name(gasQuantity) + " test - " << *iPts << std::endl;

    boost::timer gasTimer;
    gasTree.execute_tree();

    if( timings ) std::cout << "PoKiTT  " + property_name(gasQuantity) + " time " << gasTimer.elapsed() << std::endl;

#ifdef ENABLE_CUDA
    gasField.set_device_as_active( CPU_INDEX );
#endif

    SpatFldPtr<CellField> canteraResult = get_cantera_result( timings, gasQuantity, *gasMix, *iPts, temp, xcoord );

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

  // parse the command line options input describing the problem
  try {
    po::options_description desc("Supported Options");
    desc.add_options()
           ( "help", "print help message" )
           ( "xml-input-file", po::value<std::string>(&inputFileName), "Cantera xml input file name" )
           ( "phase", po::value<std::string>(&inpGroup), "name of phase in Cantera xml input file" )
           ( "timings", "Generate comparison timings between Cantera and PoKiTT across several problem sizes" );

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

    TestHelper status( false );
    status( driver( timings, RHO ) );
    status( driver( timings, P   ) );

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
