/*
 * Temperature_test.cpp
 *
 *  Created on: August 21, 2014
 *      Author: Nathan Yonkee
 */

#include <iostream>
#include <stdio.h>
#include <fstream>
#include "TestHelper.h"

#include <pokitt/thermo/Temperature.h>
#include <test/TemperaturePowers.h>
#include <pokitt/MixtureMolWeight.h>

#include <expression/ExprLib.h>

#include <spatialops/structured/Grid.h>
#include <spatialops/structured/FieldComparisons.h>
#include <spatialops/util/TimeLogger.h>

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>

#include <cantera/kernel/ct_defs.h> // contains value of Cantera::GasConstant
#include <cantera/IdealGasMix.h>

namespace So = SpatialOps;
typedef So::SVolField   CellField;

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
  switch(e){
  case H  : return "enthalpy";
  case E0 : return "E0";
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

void calculate_energy(CellField& energy, CellField& temp, Cantera_CXX::IdealGasMix& gasMix, const EnergyType energyType, const int nPts){
# ifdef ENABLE_CUDA
  energy.set_device_as_active( CPU_INDEX );
  temp.set_device_as_active( CPU_INDEX );
# endif

  const int nSpec=gasMix.nSpecies();
  const double refPressure=gasMix.pressure();
  const std::vector< std::vector<double> > massFracs = mass_fracs( nPts, nSpec);

  CellField::const_iterator iTemp = temp.begin();
  std::vector< std::vector<double> >::const_iterator iMass = massFracs.begin();
  CellField::iterator iEnergyEnd = energy.end();
  for( CellField::iterator iEnergy = energy.begin(); iEnergy!=iEnergyEnd; ++iTemp, ++iMass, ++iEnergy ){
    gasMix.setState_TPY( *iTemp, refPressure, &(*iMass)[0]);
    switch(energyType){
    case H  : *iEnergy=gasMix.enthalpy_mass();  break;
    case E0 : *iEnergy=gasMix.intEnergy_mass(); break;
    } // switch(energyType)
  }
}

//==============================================================================

So::SpatFldPtr<CellField>
get_cantera_result( const bool timings,
                    const EnergyType energyType,
                    Cantera_CXX::IdealGasMix& gasMix,
                    const int nPts,
                    CellField& temp,
                    CellField& mixMW,
                    CellField& xcoord,
                    const CellField& energy)
{
# ifdef ENABLE_CUDA
  temp.set_device_as_active  ( CPU_INDEX );
  xcoord.set_device_as_active( CPU_INDEX );
  mixMW.set_device_as_active ( CPU_INDEX );
# endif
  const int nSpec=gasMix.nSpecies();
  const double refPressure=gasMix.pressure();

  const std::vector< std::vector<double> > massFracs = mass_fracs( nPts, nSpec);
  std::vector<double> tVec;
  for( size_t i=0; i<nPts+2; ++i)
    tVec.push_back( 525.0 + 950.0 * (i-0.5)/ nPts);
  So::SpatFldPtr<CellField> canteraResult = So::SpatialFieldStore::get<CellField>(temp, CPU_INDEX);
  So::SpatFldPtr<CellField> canteraVolume = So::SpatialFieldStore::get<CellField>(temp, CPU_INDEX);
  *canteraVolume <<= Cantera::GasConstant * (500.0 + 1000 * xcoord) * mixMW / refPressure;

  std::vector< std::vector<double> >::const_iterator iMass = massFracs.begin();
  std::vector<double>::iterator iTemp = tVec.begin();
  CellField::const_iterator iEnergy = energy.begin();
  CellField::const_iterator iVolume = canteraVolume->begin();
  Timer tTime;
  tTime.start();
  for(CellField::iterator iCant = canteraResult->begin(); iCant!=canteraResult->end();++iVolume, ++iEnergy, ++iTemp, ++iMass, ++iCant){
    gasMix.setState_TPY( *iTemp, refPressure, &(*iMass)[0]);
    switch( energyType ){
    case H:  gasMix.setState_HP( *iEnergy, refPressure ); break;
    case E0: gasMix.setState_UV( *iEnergy, *iVolume    ); break;
    }
    *iCant=gasMix.temperature();
  }
  tTime.stop();
  if( timings ) std::cout << "Cantera T from " + energy_name(energyType) + " time " << tTime.elapsed_time() << std::endl;
  return canteraResult;
}

//==============================================================================

bool TPowers_equal( TestHelper& status, So::SpatFldPtr<CellField> canteraResult, Expr::FieldMgrSelector<CellField>::type& cellFM, const Expr::TagList& tPowerTags ){

  So::SpatFldPtr<CellField> canteraTPower = So::SpatialFieldStore::get<CellField>(*canteraResult, CPU_INDEX);
# ifdef ENABLE_CUDA
  BOOST_FOREACH( const Expr::Tag& tPowerTag, tPowerTags){
    CellField& tpow = cellFM.field_ref(tPowerTag);
    tpow.set_device_as_active( CPU_INDEX );
  }
# endif

  *canteraTPower <<= *canteraResult * *canteraResult;
  CellField& t2 = cellFM.field_ref(tPowerTags[0]);
  status( field_equal(t2, *canteraTPower , 5e-4), tPowerTags[0].name() );

  *canteraTPower <<= *canteraResult * *canteraResult * *canteraResult;
  CellField& t3 = cellFM.field_ref(tPowerTags[1]);
  status( field_equal(t3, *canteraTPower , 5e-4), tPowerTags[1].name() );

  *canteraTPower <<= *canteraResult * *canteraResult * *canteraResult * *canteraResult;
  CellField& t4 = cellFM.field_ref(tPowerTags[2]);
  status( field_equal(t4, *canteraTPower , 5e-4), tPowerTags[2].name() );

  *canteraTPower <<= *canteraResult * *canteraResult * *canteraResult * *canteraResult * *canteraResult;
  CellField& t5 = cellFM.field_ref(tPowerTags[3]);
  status( field_equal(t5, *canteraTPower , 5e-4), tPowerTags[3].name() );

  *canteraTPower <<= 1 / *canteraResult;
  CellField& recipT = cellFM.field_ref(tPowerTags[4]);
  status( field_equal(recipT, *canteraTPower , 5e-4), tPowerTags[4].name() );

  *canteraTPower <<=  1 / *canteraResult / *canteraResult;
  CellField& recipRecipT = cellFM.field_ref(tPowerTags[5]);
  status( field_equal(recipRecipT, *canteraTPower , 5e-4), tPowerTags[5].name() );

  return status.ok();
}

//==============================================================================

bool driver( bool timings, EnergyType energyType)
{
  TestHelper status( !timings );
  Cantera_CXX::IdealGasMix* const gasMix = CanteraObjects::get_gasmix();

  const int nSpec=gasMix->nSpecies();
  const std::vector<double>& molecularWeights = gasMix->molecularWeights();

  typedef Expr::PlaceHolder <CellField> MassFracs;
  typedef Expr::PlaceHolder <CellField> Energy;
  typedef Expr::PlaceHolder <CellField> KineticEnergy;

  Expr::TagList yiTags;
  for( size_t n=0; n<nSpec; ++n ){
    std::ostringstream name;
    name << "yi_" << n;
    yiTags.push_back( Expr::Tag(name.str(),Expr::STATE_NONE) );
  }
  const Expr::Tag energyTag( energy_name(energyType) , Expr::STATE_NONE);
  const Expr::Tag keTag( "kinetic energy", Expr::STATE_NONE);
  Expr::TagList tPowerTags;
  const Expr::Tag tTag( "Temperature", Expr::STATE_NONE);

  Expr::ExpressionFactory exprFactory;

  BOOST_FOREACH( const Expr::Tag& yiTag, yiTags){
    exprFactory.register_expression( new MassFracs::Builder(yiTag) );
  }
  exprFactory.register_expression( new Energy::Builder (energyTag) );

  Expr::ExpressionID temp_id;
  switch( energyType ){
  case H:
    typedef Temperature <CellField> Temperature;
    temp_id = exprFactory.register_expression( new Temperature ::Builder (tTag, yiTags, energyTag) );
    tPowerTags = Temperature::temperature_powers_tags();
    break;
  case E0:
    typedef TemperatureFromE0 <CellField> TemperatureE0;
    exprFactory.register_expression( new KineticEnergy::Builder (keTag) );
    temp_id = exprFactory.register_expression( new TemperatureE0 ::Builder (tTag, yiTags, energyTag, keTag) );
    tPowerTags = TemperatureE0::temperature_powers_tags();
    break;
  }

  std::vector<int> ptvec;
  if( timings ){
    ptvec.push_back(8*8*8);
    ptvec.push_back(16*16*16);
    ptvec.push_back(32*32*32);
    ptvec.push_back(64*64*64);
    ptvec.push_back(128*128*128);
  }
  else ptvec.push_back(10);

  for( std::vector<int>::iterator iPts = ptvec.begin(); iPts!= ptvec.end(); ++iPts){

    Expr::ExpressionTree tTree( temp_id , exprFactory, 0 );
    {
      std::ofstream tGraph( ("TemperatureFrom"+energy_name(energyType)+".dot").c_str() );
      tTree.write_tree( tGraph );
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

    tTree.register_fields( fml );
    fml.allocate_fields( Expr::FieldAllocInfo( npts, 0, 0, false, false, false ) );
    tTree.bind_fields( fml );

    using namespace SpatialOps;
    Expr::FieldMgrSelector<CellField>::type& cellFM = fml.field_manager< CellField>();

    CellField& temp = cellFM.field_ref(tTag);
    temp <<= 500.0 + 1000 * xcoord;

    SpatFldPtr<CellField> sum  = SpatialFieldStore::get<CellField>(temp);
    *sum<<=0.0;
    for( size_t n=0; n<nSpec; ++n ){
      CellField& yi = cellFM.field_ref(yiTags[n]);
      yi <<= n + 1 + xcoord;
      *sum <<= *sum + yi;
    }
    BOOST_FOREACH( const Expr::Tag& yiTag, yiTags){
      CellField& yi = cellFM.field_ref(yiTag);
      yi <<= yi / *sum;
    }

    SpatFldPtr<CellField> mixMW = SpatialFieldStore::get<CellField>(temp);
    for( size_t n=0; n<nSpec; ++n ){
      CellField& yi = cellFM.field_ref(yiTags[n]);
      *mixMW <<= *mixMW + yi / molecularWeights[n];
    }
    *mixMW <<= 1.0 / *mixMW;

    CellField& energy = cellFM.field_ref(energyTag);
    calculate_energy( energy, temp, *gasMix, energyType, *iPts);

    if(energyType == E0){
      CellField& ke = cellFM.field_ref(keTag);
      ke <<= 0.0;
    }

    temp <<= temp + 25 - 50 * xcoord;

    tTree.lock_fields( fml );  // prevent fields from being deallocated so that we can get them after graph execution.

    if( timings ) std::cout << std::endl << energy_name(energyType) << " test - " << *iPts << std::endl;

    Timer tTimer;
    tTimer.start();
    tTree.execute_tree();
    tTimer.stop();
    if( timings ) std::cout << "PoKiTT  T from " + energy_name(energyType) + " time " << tTimer.elapsed_time() << std::endl;

    SpatFldPtr<CellField> canteraResult = get_cantera_result( timings, energyType, *gasMix, *iPts, temp, *mixMW, xcoord, energy );

    status( field_equal(temp, *canteraResult, 5e-4), tTag.name() );

    status( TPowers_equal( status, canteraResult, cellFM, tPowerTags ), "TPowers" );

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

    TestHelper status( !timings );
    status( driver( timings, H  ), "T from " + energy_name(H ) );
    std::cout << std::endl;
    status( driver( timings, E0 ), "T from " + energy_name(E0) );

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
