/*
 * Transport_test.cpp
 *
 *  Created on: September 9, 2014
 *      Author: Nathan Yonkee
 */

#include <iostream>
#include <stdio.h>
#include <fstream>
#include "TestHelper.h"

#include <pokitt/transport/ThermalCondMix.h>
#include <pokitt/transport/DiffusionCoeffMix.h>
#include <pokitt/transport/ViscosityMix.h>
#include <pokitt/MixtureMolWeight.h>

#include <expression/ExprLib.h>

#include <spatialops/structured/Grid.h>
#include <spatialops/structured/FieldComparisons.h>

#include <boost/lexical_cast.hpp>
#include <boost/timer.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>

#include <cantera/IdealGasMix.h>

namespace So = SpatialOps;
typedef   So::SVolField CellField;

namespace Cantera_CXX{ class IdealGasMix; } //location of polynomial

namespace po = boost::program_options;

//==============================================================================

enum TransportQuantity{
  DIFF_MASS,
  DIFF_MOL,
  TCOND,
  VISC
};

std::string transport_name( const TransportQuantity q )
{
  switch(q){
    case DIFF_MASS : return "Diffusion Coefficient Mass";
    case DIFF_MOL : return "Diffusion Coefficient Mol";
    case TCOND: return "Thermal Conductivity";
    case VISC : return "Viscosity";
  }
}

//==============================================================================

const Expr::TagList make_trans_tags( TransportQuantity transportQuantity, const int nSpec ){
  Expr::TagList transTags;
  switch( transportQuantity ){
  case DIFF_MASS:
    for( size_t n=0; n<nSpec; ++n )
      transTags.push_back( Expr::Tag( "DiMass" + boost::lexical_cast<std::string>(n), Expr::STATE_NONE ) );
    break;
  case DIFF_MOL:
    for( size_t n=0; n<nSpec; ++n )
      transTags.push_back( Expr::Tag( "DiMol" + boost::lexical_cast<std::string>(n), Expr::STATE_NONE ) );
    break;
  case TCOND:
    transTags.push_back(   Expr::Tag( "Thermal Conductivity Mix",                 Expr::STATE_NONE ) );
    break;
  case VISC:
    transTags.push_back(   Expr::Tag( "Viscosity Mix",                            Expr::STATE_NONE ) );
  }
  return transTags;
}

//==============================================================================

std::vector< So::SpatFldPtr<CellField> >
get_cantera_results( const bool timings,
                     const TransportQuantity transportQuantity,
                     Cantera::MixTransport& mixTrans,
                     const int nPts,
                     const CellField& temp )
{
  using namespace SpatialOps;

  Cantera::ThermoPhase& canteraThermo = mixTrans.thermo();
  const double refPressure = canteraThermo.pressure();
  const std::vector<double>& molecularWeights = canteraThermo.molecularWeights();
  const int nSpec = canteraThermo.nSpecies();

  std::vector< std::vector<double> > massFracs;
  for( size_t i=0; i<nPts+2; ++i){
    std::vector<double> massFrac;
    double sum = 0.0;
    for( size_t n=0; n<nSpec; ++n){
      massFrac.push_back(1 + n + (i-0.5)/ nPts);
      sum+=massFrac[n];
    }
    for( size_t n=0; n<nSpec; ++n)
      massFrac[n] = massFrac[n]/sum;
    massFracs.push_back(massFrac);
  }

  std::vector< SpatFldPtr<CellField> > canteraResults;
  switch( transportQuantity ){
  case DIFF_MASS:
    for( size_t n=0; n < nSpec; ++n){
      canteraResults.push_back(SpatialFieldStore::get<CellField>(temp));
    }
    break;
  case DIFF_MOL:
    for( size_t n=0; n < nSpec; ++n){
      canteraResults.push_back(SpatialFieldStore::get<CellField>(temp));
    }
    break;
  default:
    canteraResults.push_back(SpatialFieldStore::get<CellField>(temp));
  }

  CellField::const_iterator iTemp = temp.begin();
  std::vector< std::vector<double> >::iterator iMass = massFracs.begin();
  boost::timer transportTimer;

  if( transportQuantity == DIFF_MASS || transportQuantity == DIFF_MOL){
    std::vector<double> d_result(nSpec,0.0);
    for( size_t i=0; i<nPts+2; ++iTemp, ++iMass, ++i){
      canteraThermo.setState_TPY( *iTemp, refPressure, &(*iMass)[0]);
      switch(transportQuantity ){
      case DIFF_MASS:
      mixTrans.getMixDiffCoeffsMass(&d_result[0]);
      for( size_t n=0; n<nSpec; ++n){
        (*canteraResults[n])[i] = d_result[n];
      }
      break;
      case DIFF_MOL:
      mixTrans.getMixDiffCoeffs(&d_result[0]);
      for( size_t n=0; n<nSpec; ++n){
        (*canteraResults[n])[i] = d_result[n];
      }
      break;
      }
    }
  }
  else{
    CellField::iterator iCantEnd = canteraResults[0]->end();
    for(CellField::iterator iCant = canteraResults[0]->begin(); iCant!=iCantEnd; ++iTemp, ++iMass, ++iCant){
      canteraThermo.setState_TPY( *iTemp, refPressure, &(*iMass)[0]);
      switch( transportQuantity ){
      case TCOND:
        *iCant = mixTrans.thermalConductivity();
        break;
      case VISC:
        *iCant = mixTrans.viscosity();
        break;
      }
    }
  }

  if( timings ) std::cout << "Cantera " + transport_name(transportQuantity) + " time " << transportTimer.elapsed() << std::endl;
  return canteraResults;
}

//==============================================================================

bool driver( const bool timings,
             const TransportQuantity transportQuantity )
{
  TestHelper status( !timings );

  Cantera::Transport* transport = CanteraObjects::get_transport();
  Cantera::MixTransport* mixTrans;

  if( transport->model() ==210 || transport->model()==211)
    mixTrans = dynamic_cast<Cantera::MixTransport*>( transport );
  else {
    std::cout<<"error, transport not mixture\ntransport model is " << transport->model() << std::endl;
    return -1;
  }

  const int nSpec=mixTrans->thermo().nSpecies();
  size_t n;
  const double refPressure=mixTrans->thermo().pressure();
  const std::vector<double>& molecularWeights = mixTrans->thermo().molecularWeights();

  typedef Expr::PlaceHolder <CellField> Temp;
  typedef Expr::PlaceHolder <CellField> Pressure;
  typedef Expr::PlaceHolder <CellField> MassFracs;
  typedef MixtureMolWeight  <CellField> MixtureMolWeight;

  const Expr::Tag tTag ( "Temperature", Expr::STATE_NONE );
  const Expr::Tag pTag ( "Pressure"   , Expr::STATE_NONE);
  const Expr::Tag yiTag ( "yi", Expr::STATE_NONE );
  Expr::TagList yiTags;
  for( n=0; n<nSpec; ++n ){
    std::ostringstream name;
    name << "yi_" << n;
    yiTags.push_back( Expr::Tag( name.str(), Expr::STATE_NONE ) );
  }
  const Expr::Tag mmwTag( "mmw", Expr::STATE_NONE);
  const Expr::TagList transportTags = make_trans_tags( transportQuantity, nSpec );

  Expr::ExpressionFactory exprFactory;

  exprFactory.register_expression( new Temp::Builder(tTag) );
  exprFactory.register_expression( new Pressure ::Builder (pTag                 ) );
  BOOST_FOREACH( Expr::Tag yiTag, yiTags){
    exprFactory.register_expression( new MassFracs::Builder (yiTag) );
  }
  exprFactory.register_expression( new MixtureMolWeight::Builder( mmwTag, yiTag, molecularWeights));
  Expr::ExpressionID transportID;
  switch( transportQuantity ){
  case DIFF_MASS:
    typedef DiffusionCoeff <CellField> DiffusionCoeffMix;
    transportID = exprFactory.register_expression( new DiffusionCoeffMix::Builder (transportTags, tTag, pTag ,yiTag, mmwTag) );
    break;
  case DIFF_MOL:
    typedef DiffusionCoeffMol <CellField> DiffusionCoeffMixMol;
    transportID = exprFactory.register_expression( new DiffusionCoeffMixMol::Builder (transportTags, tTag, pTag ,yiTag, mmwTag) );
    break;
  case TCOND:
    typedef ThermalConductivity <CellField> ThermalConductivityMix;
    transportID = exprFactory.register_expression( new ThermalConductivityMix::Builder (transportTags[0] ,tTag ,yiTag, mmwTag) );
    break;
  case VISC:
    typedef Viscosity <CellField> ViscosityMix;
    transportID = exprFactory.register_expression( new ViscosityMix::Builder (transportTags[0] ,tTag ,yiTag) );
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
  else{
    ptvec.push_back(10);
  }

  for( std::vector<int>::iterator iPts = ptvec.begin(); iPts!= ptvec.end(); ++iPts){

    Expr::ExpressionTree transportTree( transportID, exprFactory, 0 );

    {
      std::ofstream fout( (transport_name(transportQuantity) + ".dot").c_str() );
      transportTree.write_tree(fout);
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

    transportTree.register_fields( fml );
    fml.allocate_fields( Expr::FieldAllocInfo( npts, 0, 0, false, false, false ) );
    transportTree.bind_fields( fml );

    using namespace SpatialOps;
    Expr::FieldMgrSelector<CellField>::type& cellFM = fml.field_manager<CellField>();

    CellField& temp = cellFM.field_ref(tTag);
    temp <<= 500.0 + 1000.0 * xcoord;

    if( transportQuantity == DIFF_MASS || transportQuantity == DIFF_MOL ){
      CellField& p = cellFM.field_ref(pTag);
      p <<= refPressure;
    }

    SpatFldPtr<CellField> sum  = SpatialFieldStore::get<CellField>(temp);
    *sum<<=0.0;
    for( n=0; n<nSpec; ++n ){
      CellField& yi = cellFM.field_ref(yiTags[n]);
      yi <<= n + 1 + xcoord;
      *sum <<= *sum + yi;
    }
    BOOST_FOREACH( Expr::Tag yiTag, yiTags){
      CellField& yi = cellFM.field_ref(yiTag);
      yi <<= yi / *sum;
    }

    transportTree.lock_fields( fml );  // prevent fields from being deallocated so that we can get them after graph execution.

    if( timings ) std::cout << std::endl << transport_name(transportQuantity) << " test - " << *iPts << std::endl;

    boost::timer transportTimer;
    transportTree.execute_tree();

    if( timings ) std::cout << "PoKiTT  " + transport_name(transportQuantity) + " time " << transportTimer.elapsed() << std::endl;

#   ifdef ENABLE_CUDA
    BOOST_FOREACH( Expr::Tag transportTag, transportTags){
      CellField& trans = fml.field_manager<CellField>().field_ref(transportTag);
      trans.add_device(CPU_INDEX);
    }
    temp.set_device_as_active( CPU_INDEX );
#   endif

    const std::vector< SpatFldPtr<CellField> > canteraResults = get_cantera_results( timings,
                                                                                     transportQuantity,
                                                                                     *mixTrans,
                                                                                     *iPts,
                                                                                     temp );

    std::vector< SpatFldPtr<CellField> >::const_iterator iCantera = canteraResults.begin();
    BOOST_FOREACH( const Expr::Tag& transportTag, transportTags ){
      CellField& transport = cellFM.field_ref(transportTag);
      status( field_equal( transport, **iCantera, 1e-12 ), transportTag.name() );
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

    TestHelper status( !timings );
    status( driver( timings, DIFF_MASS  ), transport_name(DIFF_MASS ) );
    status( driver( timings, DIFF_MOL  ), transport_name(DIFF_MOL ) );
    status( driver( timings, TCOND ), transport_name(TCOND) );
    status( driver( timings, VISC  ), transport_name(VISC ) );

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
