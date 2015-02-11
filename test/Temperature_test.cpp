/*
 * Temperature_test.cpp
 *
 *  Created on: August 21, 2014
 *      Author: Nathan Yonkee
 */

#include <iostream>
#include "TestHelper.h"
#include "LinearMassFracs.h"
#include <pokitt/thermo/Temperature.h>
#include <pokitt/thermo/SpecificVol.h>
#include <pokitt/thermo/Enthalpy.h>
#include <pokitt/thermo/InternalEnergy.h>
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
typedef SO::SVolField  CellField;
typedef SO::SpatFldPtr<CellField> CellFieldPtrT;

namespace Cantera_CXX{ class IdealGasMix; } //location of polynomial

namespace po = boost::program_options;

using namespace pokitt;

//==============================================================================

enum EnergyType{
  H,
  E0
};

std::string energy_name( const EnergyType e )
{
  std::string name;
  switch(e){
  case H  : name = "Enthalpy"; break;
  case E0 : name = "E0";       break;
  }
  return name;
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

CellFieldPtrT
get_cantera_result( const bool timings,
                    const size_t canteraReps,
                    const EnergyType energyType,
                    Cantera_CXX::IdealGasMix& gasMix,
                    Expr::FieldManagerList& fml,
                    const Expr::Tag& tTag,
                    const Expr::TagList& yiTags,
                    const Expr::Tag& energyTag,
                    const Expr::Tag& nuTag,
                    const Expr::Tag& pTag)
{
  const std::vector< std::vector<double> > massFracs = extract_mass_fracs( yiTags, fml );

  CellField& temp   = fml.field_ref< CellField >( tTag );
  CellField& energy = fml.field_ref< CellField >( energyTag );
  CellField& press  = fml.field_ref< CellField >( pTag );
  CellField& volume = fml.field_ref< CellField >( nuTag );
# ifdef ENABLE_CUDA
  temp.set_device_as_active  ( CPU_INDEX );
  energy.set_device_as_active( CPU_INDEX );
  press.set_device_as_active ( CPU_INDEX );
  volume.set_device_as_active( CPU_INDEX );
# endif
  CellFieldPtrT canteraResult = SO::SpatialFieldStore::get<CellField>(temp, CPU_INDEX);

  std::vector< std::vector<double> >::const_iterator iMass;
  CellField::const_iterator                          iTemp;
  CellField::const_iterator                          iEnergy;
  CellField::const_iterator                          iVolume;
  CellField::const_iterator                          iPress;
  CellField::iterator                                iCant;
  Timer tTime;
  tTime.start();
  for( size_t rep=0; rep < canteraReps; ++rep ){
    iEnergy = energy.begin();
    iVolume = volume.begin();
    iPress  = press.begin();
    iTemp   = temp.begin();
    iMass   = massFracs.begin();
    for( iCant = canteraResult->begin(); iCant != canteraResult->end(); ++iVolume, ++iPress, ++iEnergy, ++iTemp, ++iMass, ++iCant){
      gasMix.setMassFractions_NoNorm( &(*iMass)[0] );
      switch( energyType ){
      case H:  gasMix.setState_HP( *iEnergy, *iPress  ); break;
      case E0: gasMix.setState_UV( *iEnergy, *iVolume ); break;
      }
      *iCant=gasMix.temperature();
    }
  }
  tTime.stop();
  if( timings ) std::cout << "Cantera T from " + energy_name(energyType) + " time " << tTime.elapsed_time()/canteraReps << std::endl;
  return canteraResult;
}

//==============================================================================

bool driver( const bool timings,
             const size_t pokittReps,
             const size_t canteraReps,
             const EnergyType energyType)
{
  TestHelper status( !timings );
  Cantera_CXX::IdealGasMix* const gasMix = CanteraObjects::get_gasmix();
  const int nSpec=gasMix->nSpecies();

  const Expr::Tag xTag( "XCoord", Expr::STATE_NONE );
  Expr::TagList yiTags;
  for( size_t n=0; n<nSpec; ++n ){
    yiTags.push_back( Expr::Tag( "yi_" + boost::lexical_cast<std::string>(n), Expr::STATE_NONE ) );
  }
  const Expr::Tag tInputTag( "Input Temperature",        Expr::STATE_NONE);
  const Expr::Tag pTag     ( "Pressure",                 Expr::STATE_NONE);
  const Expr::Tag nuTag    ( "Specific Volume",          Expr::STATE_NONE);
  const Expr::Tag mmwTag   ( "Mixture Molecular Weight", Expr::STATE_NONE);
  const Expr::Tag energyTag( energy_name(energyType) ,   Expr::STATE_NONE);
  const Expr::Tag keTag    ( "kinetic energy",           Expr::STATE_NONE);
  const Expr::Tag tTag     ( "Temperature",              Expr::STATE_NONE);

  // we use an initialization tree to avoid recalculations when timing the execution
  Expr::ExpressionFactory initFactory;
  std::set<Expr::ExpressionID> initIDs;
  Expr::ExpressionID tID; // perturbed temperature
  Expr::ExpressionID eID; // input energy for T solver
  Expr::ExpressionID vID; // specific volume
  Expr::ExpressionID kID; // kinetic energy
  {
    typedef Expr::PlaceHolder     <CellField>::Builder XCoord;
    typedef       LinearMassFracs <CellField>::Builder MassFracs;
    typedef Expr::ConstantExpr    <CellField>::Builder Pressure;
    typedef       SpecificVol     <CellField>::Builder SpecVol;
    typedef       MixtureMolWeight<CellField>::Builder MixMolWeight;
    typedef       Enthalpy        <CellField>::Builder Enthalpy;
    typedef       InternalEnergy  <CellField>::Builder IntEnergy;
    typedef Expr::ConstantExpr    <CellField>::Builder KinEnergy;
    typedef Expr::LinearFunction  <CellField>::Builder Temperature;
    typedef Expr::LinearFunction  <CellField>::Builder TPerturbation;

    initFactory.register_expression(        new XCoord       ( xTag )                           );
    initFactory.register_expression(        new MassFracs    ( yiTags, xTag )                   );
    initFactory.register_expression(        new Pressure     ( pTag, gasMix->pressure() )       );
    initFactory.register_expression(        new MixMolWeight ( mmwTag, yiTags )                 );
    initFactory.register_expression(        new Temperature  ( tInputTag ,xTag, 1000, 500 )     );
    vID = initFactory.register_expression(  new SpecVol      ( nuTag, tInputTag, pTag, mmwTag ) );
    tID = initFactory.register_expression(  new TPerturbation( tTag ,tInputTag, 1.1, -50 )      );
    switch( energyType ){
    case H:
      eID = initFactory.register_expression( new Enthalpy ( energyTag, tInputTag, yiTags ) );
      break;
    case E0:
      eID = initFactory.register_expression( new IntEnergy( energyTag, tInputTag, yiTags ) );
      kID = initFactory.register_expression( new KinEnergy( keTag, 0.0                   ) );
      initIDs.insert( kID );
    }
    initIDs.insert( tID );
    initIDs.insert( eID );
    initIDs.insert( vID );
  }

  Expr::ExpressionFactory execFactory;
  Expr::ExpressionID execID;
  {
    typedef Expr::PlaceHolder      <CellField>::Builder MassFracs;
    typedef Expr::PlaceHolder      <CellField>::Builder Energy;
    typedef Expr::PlaceHolder      <CellField>::Builder KineticEnergy;
    typedef       TemperatureFromE0<CellField>::Builder TemperatureE0;
    typedef       Temperature      <CellField>::Builder Temperature;

    execFactory.register_expression( new MassFracs ( yiTags    ) );
    execFactory.register_expression( new Energy    ( energyTag ) );
    switch( energyType ){
    case H:
      execID = execFactory.register_expression( new Temperature ( tTag, yiTags, energyTag ) );
      break;
    case E0:
      execFactory.register_expression(          new KineticEnergy (keTag )                         );
      execID = execFactory.register_expression( new TemperatureE0 (tTag, yiTags, energyTag, keTag) );
      break;
    }
  }

  Expr::ExpressionTree initTree( initIDs, initFactory, 0 );
  Expr::ExpressionTree execTree( execID , execFactory, 0 );
  {
    std::ofstream tGraph( ("TemperatureFrom"+energy_name(energyType)+"_init.dot").c_str() );
    initTree.write_tree( tGraph );
  }
  {
    std::ofstream tGraph( ("TemperatureFrom"+energy_name(energyType)+".dot").c_str() );
    execTree.write_tree( tGraph );
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

    if( timings ) std::cout << std::endl << "T from " << energy_name(energyType) << " test - " << nPoints << std::endl;

    Timer tTimer;
    for( size_t rep = 0; rep < pokittReps; ++rep ){
      initTree.execute_tree(); // set initial guess
      tTimer.start();
      execTree.execute_tree();
      tTimer.stop();
    }

    if( timings ) std::cout << "PoKiTT  T from " + energy_name(energyType) + " time " << tTimer.elapsed_time()/pokittReps << std::endl;

    initTree.execute_tree(); // set initial guess for Cantera
    CellFieldPtrT canteraResult = get_cantera_result( timings, canteraReps, energyType, *gasMix, fml, tTag, yiTags, energyTag, nuTag, pTag );
    execTree.execute_tree(); // resolve to compare with Cantera
    CellField& temp = fml.field_ref< CellField >( tTag );
    status( field_equal(temp, *canteraResult, 5e-4), tTag.name() );

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
    status( driver( timings, pokittReps, canteraReps, H  ), "T from " + energy_name(H ) );
    std::cout << std::endl;
    status( driver( timings, pokittReps, canteraReps, E0 ), "T from " + energy_name(E0) );

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
