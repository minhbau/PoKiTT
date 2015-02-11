/*
 * Transport_test.cpp
 *
 *  Created on: September 9, 2014
 *      Author: Nathan Yonkee
 */

#include <iostream>
#include "TestHelper.h"
#include "LinearMassFracs.h"
#include <pokitt/transport/ThermalCondMix.h>
#include <pokitt/transport/DiffusionCoeffMix.h>
#include <pokitt/transport/ViscosityMix.h>
#include <pokitt/MixtureMolWeight.h>

#include <expression/ExprLib.h>

#include <spatialops/structured/Grid.h>
#include <spatialops/structured/FieldComparisons.h>
#include <spatialops/util/TimeLogger.h>

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>

#include <cantera/IdealGasMix.h>

namespace SO = SpatialOps;
typedef   SO::SVolField CellField;
typedef SO::SpatFldPtr<CellField> CellFieldPtrT;

namespace Cantera_CXX{ class IdealGasMix; } //location of polynomial

namespace po = boost::program_options;

using namespace pokitt;

//==============================================================================

enum TransportQuantity{
  DIFF_MASS,
  DIFF_MOL,
  TCOND,
  VISC
};

std::string transport_name( const TransportQuantity q )
{
  std::string name;
  switch(q){
    case DIFF_MASS: name = "Diffusion Coefficient Mass"; break;
    case DIFF_MOL : name = "Diffusion Coefficient Mol";  break;
    case TCOND    : name = "Thermal Conductivity";       break;
    case VISC     : name = "Viscosity";                  break;
  }
  return name;
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
    transTags.push_back(   Expr::Tag( "Thermal Conductivity Mix", Expr::STATE_NONE ) );
    break;
  case VISC:
    transTags.push_back(   Expr::Tag( "Viscosity Mix",            Expr::STATE_NONE ) );
  }
  return transTags;
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

const std::vector< CellFieldPtrT >
get_cantera_results( const bool timings,
                     const size_t canteraReps,
                     const TransportQuantity transportQuantity,
                     Cantera::MixTransport& mixTrans,
                     Expr::FieldManagerList& fml,
                     const Expr::Tag& tTag,
                     const Expr::TagList& yiTags,
                     const Expr::Tag& pTag){
  Cantera::ThermoPhase& canteraThermo = mixTrans.thermo();
  const std::vector<double>& molecularWeights = canteraThermo.molecularWeights();
  const int nSpec = canteraThermo.nSpecies();

  const std::vector< std::vector<double> > massFracs = extract_mass_fracs( yiTags, fml );

    CellField& temp   = fml.field_ref< CellField >( tTag );
    CellField& press  = fml.field_ref< CellField >( pTag );
  # ifdef ENABLE_CUDA
    temp.set_device_as_active  ( CPU_INDEX );
    press.set_device_as_active ( CPU_INDEX );
  # endif

  std::vector< CellFieldPtrT > canteraResults;
  switch( transportQuantity ){
  case DIFF_MASS:
    for( size_t n=0; n < nSpec; ++n){
      canteraResults.push_back(SO::SpatialFieldStore::get<CellField>(temp));
    }
    break;
  case DIFF_MOL:
    for( size_t n=0; n < nSpec; ++n){
      canteraResults.push_back(SO::SpatialFieldStore::get<CellField>(temp));
    }
    break;
  default:
    canteraResults.push_back(SO::SpatialFieldStore::get<CellField>(temp));
  }

  std::vector< std::vector<double> >::const_iterator iMass;
  CellField::const_iterator                          iTemp;
  CellField::const_iterator                          iPress;
  CellField::iterator                                iCant;

  Timer transportTimer;
  if( transportQuantity == DIFF_MASS || transportQuantity == DIFF_MOL){
    std::vector<double> d_result(nSpec,0.0);
    transportTimer.start();
    for( size_t rep=0; rep < canteraReps; ++rep ){
      iPress = press.begin();
      iMass  = massFracs.begin();
      size_t i = 0;
      for( iTemp = temp.begin(); iTemp != temp.end(); ++iTemp, ++iPress, ++iMass, ++i){
        canteraThermo.setMassFractions_NoNorm( &(*iMass)[0] );
        canteraThermo.setState_TP( *iTemp, *iPress );
        switch(transportQuantity ){
        case DIFF_MASS:
          mixTrans.getMixDiffCoeffsMass(&d_result[0]); break;
        case DIFF_MOL:
          mixTrans.getMixDiffCoeffs    (&d_result[0]); break;
        }
        for( size_t n=0; n<nSpec; ++n){
          (*canteraResults[n])[i] = d_result[n];
        }
      }
    }
  }
  else{
    transportTimer.start();
    for( size_t rep=0; rep < canteraReps; ++rep ){
      iPress = press.begin();
      iTemp  = temp.begin();
      iMass  = massFracs.begin();
      CellField::iterator iCantEnd = canteraResults[0]->end();
      for(CellField::iterator iCant = canteraResults[0]->begin(); iCant!=iCantEnd; ++iPress, ++iTemp, ++iMass, ++iCant){
        canteraThermo.setMassFractions_NoNorm( &(*iMass)[0] );
        canteraThermo.setState_TP( *iTemp, *iPress );
        switch( transportQuantity ){
        case TCOND:
          *iCant = mixTrans.thermalConductivity(); break;
        case VISC:
          *iCant = mixTrans.viscosity(); break;
        }
      }
    }
  }
  transportTimer.stop();

  if( timings ) std::cout << "Cantera " + transport_name(transportQuantity) + " time " << transportTimer.elapsed_time()/canteraReps << std::endl;
  return canteraResults;
}

//==============================================================================

bool driver( const bool timings,
             const size_t pokittReps,
             const size_t canteraReps,
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
  Cantera::ThermoPhase& cThermo = mixTrans->thermo();
  const int nSpec=cThermo.nSpecies();

  const Expr::Tag xTag  ( "XCoord",      Expr::STATE_NONE );
  const Expr::Tag tTag  ( "Temperature", Expr::STATE_NONE);
  const Expr::Tag pTag  ( "Pressure",    Expr::STATE_NONE);
  const Expr::Tag mmwTag( "mmw",         Expr::STATE_NONE);
  Expr::TagList yiTags;
  for( size_t n=0; n<nSpec; ++n ){
    yiTags.push_back( Expr::Tag( "yi_" + boost::lexical_cast<std::string>(n), Expr::STATE_NONE ) );
  }
  const Expr::TagList transportTags = make_trans_tags( transportQuantity, nSpec );

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
    pID = initFactory.register_expression( new Pressure     ( pTag, cThermo.pressure() ));
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
    typedef DiffusionCoeff      <CellField>::Builder DiffusionCoeffMix;
    typedef DiffusionCoeffMol   <CellField>::Builder DiffusionCoeffMixMol;
    typedef ThermalConductivity <CellField>::Builder TConductivityMix;
    typedef Viscosity           <CellField>::Builder ViscosityMix;

    execFactory.register_expression( new MassFracs   ( yiTags ) );
    execFactory.register_expression( new Temperature ( tTag   ) );
    execFactory.register_expression( new Pressure    ( pTag   ) );
    execFactory.register_expression( new MixMolWeight( mmwTag ) );
    switch( transportQuantity ){
    case DIFF_MASS:
      execID = execFactory.register_expression( new DiffusionCoeffMix   ( transportTags, tTag, pTag, yiTags, mmwTag) );
      break;
    case DIFF_MOL:
      execID = execFactory.register_expression( new DiffusionCoeffMixMol( transportTags, tTag, pTag, yiTags, mmwTag) );
      break;
    case TCOND:
      execID = execFactory.register_expression( new TConductivityMix    ( transportTags[0], tTag, yiTags, mmwTag) );
      break;
    case VISC:
      execID = execFactory.register_expression( new ViscosityMix        ( transportTags[0], tTag, yiTags) );
      break;
    }
  }

  Expr::ExpressionTree initTree( initIDs, initFactory, 0 );
  Expr::ExpressionTree execTree( execID , execFactory, 0 );
  {
    std::ofstream fout( (transport_name(transportQuantity) + "_init.dot").c_str() );
    initTree.write_tree(fout);
  }
  {
    std::ofstream fout( (transport_name(transportQuantity) + ".dot").c_str() );
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

    if( timings ) std::cout << std::endl << transport_name(transportQuantity) << " test - " << nPoints << std::endl;

    initTree.execute_tree();
    Timer transportTimer;
    transportTimer.start();
    for( size_t rep = 0; rep < pokittReps; ++rep ){
      execTree.execute_tree();
    }
    transportTimer.stop();

    if( timings ) std::cout << "PoKiTT  " + transport_name(transportQuantity) + " time " << transportTimer.elapsed_time()/pokittReps << std::endl;

    const std::vector< CellFieldPtrT > canteraResults = get_cantera_results( timings,
                                                                             canteraReps,
                                                                             transportQuantity,
                                                                             *mixTrans,
                                                                             fml,
                                                                             tTag,
                                                                             yiTags,
                                                                             pTag );
#   ifdef ENABLE_CUDA
    BOOST_FOREACH( const Expr::Tag& transportTag, transportTags){
      CellField& tran = fml.field_ref< CellField >( transportTag );
      tran.add_device(CPU_INDEX);
    }
#   endif

    std::vector< CellFieldPtrT >::const_iterator iCantera = canteraResults.begin();
    BOOST_FOREACH( const Expr::Tag& transportTag, transportTags ){
      CellField& tran = fml.field_ref< CellField >( transportTag );
      status( field_equal( tran, **iCantera, 1e-12 ), transportTag.name() );
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
    status( driver( timings, pokittReps, canteraReps, DIFF_MASS ), transport_name(DIFF_MASS ) );
    status( driver( timings, pokittReps, canteraReps, DIFF_MOL  ), transport_name(DIFF_MOL  ) );
    status( driver( timings, pokittReps, canteraReps, TCOND     ), transport_name(TCOND     ) );
    status( driver( timings, pokittReps, canteraReps, VISC      ), transport_name(VISC      ) );

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
