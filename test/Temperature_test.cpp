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

#include <boost/lexical_cast.hpp>
#include <boost/timer.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>

#include <cantera/kernel/ct_defs.h> // contains value of Cantera::GasConstant
#include <cantera/kernel/speciesThermoTypes.h> // contains definitions for which polynomial is being used
#include <cantera/IdealGasMix.h>

namespace So = SpatialOps;
typedef So::SVolField   CellField;

namespace Cantera_CXX{ class IdealGasMix; } //location of polynomial

namespace po = boost::program_options;

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

bool driver( bool timings, EnergyType energyType)
{
  TestHelper status( !timings );
  Cantera_CXX::IdealGasMix* const gasMix = CanteraObjects::get_gasmix();

  const int nSpec=gasMix->nSpecies();
  const double refPressure=gasMix->pressure();
  const std::vector<double>& molecularWeights = gasMix->molecularWeights();
  size_t n;

  typedef Expr::PlaceHolder <CellField> MassFracs;
  typedef MixtureMolWeight  <CellField> MixtureMolWeight;
  typedef Expr::PlaceHolder <CellField> Energy;
  typedef Expr::PlaceHolder <CellField> KineticEnergy;

  const Expr::Tag yiTag ( "yi", Expr::STATE_NONE );
  Expr::TagList yiTags;
  for( size_t n=0; n<nSpec; ++n ){
    std::ostringstream name;
    name << yiTag.name() << "_" << n;
    yiTags.push_back( Expr::Tag(name.str(),yiTag.context()) );
  }
  const Expr::Tag mmwTag ( "Mixture Mol Weight", Expr::STATE_NONE );
  const Expr::Tag energyTag  ( energy_name(energyType) , Expr::STATE_NONE);
  const Expr::Tag keTag  ( "kinetic energy", Expr::STATE_NONE);

  Expr::ExpressionFactory exprFactory;

  BOOST_FOREACH( Expr::Tag yiTag, yiTags){
    exprFactory.register_expression( new MassFracs::Builder (yiTag) );
  }
  exprFactory.register_expression( new Energy::Builder (energyTag) );
  exprFactory.register_expression( new MixtureMolWeight::Builder ( mmwTag, yiTag, molecularWeights));
  Expr::TagList tPowerTags;

  const Expr::Tag tTag ( "Temperature", Expr::STATE_NONE);
  Expr::ExpressionID temp_id;

  switch( energyType ){
  case H:
    typedef Temperature       <CellField> Temperature;
    temp_id = exprFactory.register_expression( new Temperature ::Builder (tTag, yiTag, energyTag) );
    tPowerTags = Temperature::temperature_powers_tags();
    break;
  case E0:
    typedef TemperatureFromE0 <CellField> TemperatureE0;
    exprFactory.register_expression( new KineticEnergy::Builder (keTag) );
    temp_id = exprFactory.register_expression( new TemperatureE0 ::Builder (tTag, yiTag, energyTag, keTag) );
    tPowerTags = Temperature::temperature_powers_tags();
    break;
  }

    std::vector<int> ptvec;
#   ifdef TIMINGS
    ptvec.push_back(8*8*8);
    ptvec.push_back(16*16*16);
    ptvec.push_back(32*32*32);
    ptvec.push_back(64*64*64);
#   ifdef ENABLE_CUDA
    ptvec.push_back(128*128*128);
#   endif
#   else
    ptvec.push_back(10);
#   endif

    for( std::vector<int>::iterator iPts = ptvec.begin(); iPts!= ptvec.end(); ++iPts){
      size_t i;

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
#     ifdef ENABLE_CUDA
      xcoord.add_device( GPU_INDEX );
#     endif

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
      for( n=0; n<nSpec; ++n ){
        CellField& yi = cellFM.field_ref(yiTags[n]);
        yi <<= n + 1 + xcoord;
        *sum <<= *sum + yi;
      }
      BOOST_FOREACH( Expr::Tag yiTag, yiTags){
        CellField& yi = cellFM.field_ref(yiTag);
        yi <<= yi / *sum;
      }

      SpatFldPtr<CellField> mixMW = SpatialFieldStore::get<CellField>(temp);
      for( size_t n=0; n<nSpec; ++n ){
        CellField& yi = cellFM.field_ref(yiTags[n]);
        *mixMW <<= *mixMW + yi / molecularWeights[n];
      }
      *mixMW <<= 1.0 / *mixMW;

      std::vector< std::vector<double> > massfracs;
      for( i=0; i<*iPts+2; ++i){
        std::vector<double> massfrac;
        double sum = 0.0;
        for( n=0; n<nSpec; ++n){
          massfrac.push_back(1 + n + (i-0.5)/ *iPts);
          sum+=massfrac[n];
        }
        for( n=0; n<nSpec; ++n)
          massfrac[n] = massfrac[n]/sum;
        massfracs.push_back(massfrac);
      }

      CellField& energy = cellFM.field_ref( energyTag );
      CellField::const_iterator iTemp = temp.begin();
      std::vector< std::vector<double> >::iterator iMass = massfracs.begin();
      CellField::iterator iEnergyEnd = energy.end();
      CellField::iterator iEnergy;
      for( iEnergy = energy.begin(); iEnergy!=iEnergyEnd; ++iTemp, ++iMass, ++iEnergy ){
        gasMix->setState_TPY( *iTemp, refPressure, &(*iMass)[0]);
        switch(energyType){
        case H  : *iEnergy=gasMix->enthalpy_mass();   break;
        case E0 : *iEnergy=gasMix->intEnergy_mass(); break;
        } // switch(energyType)
      }

      if(energyType == E0){
        CellField& ke = cellFM.field_ref(keTag);
        ke <<= 0.0;
      }

      temp <<= temp + 25 - 50 * xcoord;

//      calculate_internal_energy( e0Tag, tTag, yiTags, cellFM, gasMix);

//      std::vector<double> e0_massVec;
//      i=0;
//      for( itemp = tVec.begin(); itemp < itend; ++itemp, ++i){
//        gasMix->setState_TPY(*itemp,refPressure, &massfracs[i][0]);
//        e0_massVec.push_back( gasMix -> intEnergy_mass());
//      }
//
//      CellField& ke = cellFM.field_ref(keTag);
//      ke <<= 0.0;



      std::vector<double> tVecDiff;
      for( i=0; i<*iPts+2; ++i)
        tVecDiff.push_back( 525.0 + 950.0 * (i-0.5)/ *iPts);

      tTree.lock_fields( fml );  // prevent fields from being deallocated so that we can get them after graph execution.
//      te0Tree.lock_fields( fml );

#     ifdef TIMINGS
      std::cout << std::endl << setup.inputFile << " - " << *iPts << std::endl;
#     endif

      boost::timer tTimer;
      tTree.execute_tree();
#     ifdef TIMINGS
      std::cout << "PoKiTT temperature time  " << tTimer.elapsed() << std::endl;
#     endif

#     ifdef TIMINGS
      std::cout << "PoKiTT temperature from e0 time  " << te0Timer.elapsed() << std::endl;
#     endif

#     ifdef ENABLE_CUDA
      temp.add_device(CPU_INDEX);
#     endif

      std::vector< std::vector<double> >::iterator imass = massfracs.begin();
      std::vector<double>::iterator itempd = tVecDiff.begin();
      SpatFldPtr<CellField> canteraResult  = SpatialFieldStore::get<CellField>(temp);
      SpatFldPtr<CellField> canteraVolume = SpatialFieldStore::get<CellField>(temp);
      *canteraVolume <<= Cantera::GasConstant * (500.0 + 1000 * xcoord) * *mixMW / refPressure;
      i=0;
      iEnergy = energy.begin();
      CellField::iterator iVolume = canteraVolume->begin();
      for(CellField::iterator icant = canteraResult->begin(); icant!=canteraResult->end();++iVolume, ++iEnergy, ++itempd, ++imass, ++icant, ++i){
        gasMix->setState_TPY( *itempd, refPressure, &(*imass)[0]);
        switch( energyType ){
        case H: gasMix->setState_HP( *iEnergy, refPressure); break;
        case E0: gasMix->setState_UV( *iEnergy, *iVolume);
        }

        *icant=gasMix->temperature();
      }
      status( field_equal(temp, *canteraResult, 1e-6), "temperature from h");

      SpatFldPtr<CellField> canteraTPower = SpatialFieldStore::get<CellField>(*canteraResult);
      *canteraTPower <<= *canteraResult * *canteraResult;
      CellField& t2 = cellFM.field_ref(tPowerTags[1]);
      status( field_equal(t2, *canteraTPower , 1e-6), "temperature^2");

      *canteraTPower <<= *canteraResult * *canteraResult * *canteraResult;
      CellField& t3 = cellFM.field_ref(tPowerTags[2]);
      status( field_equal(t3, *canteraTPower , 1e-6), "temperature^3");

      *canteraTPower <<= *canteraResult * *canteraResult * *canteraResult * *canteraResult;
      CellField& t4 = cellFM.field_ref(tPowerTags[3]);
      status( field_equal(t4, *canteraTPower , 1e-6), "temperature^4");

      *canteraTPower <<= *canteraResult * *canteraResult * *canteraResult * *canteraResult * *canteraResult;
      CellField& t5 = cellFM.field_ref(tPowerTags[4]);
      status( field_equal(t5, *canteraTPower , 1e-6), "temperature^5");

      *canteraTPower <<= 1 / *canteraResult;
      CellField& recipT = cellFM.field_ref(tPowerTags[5]);
      status( field_equal(recipT, *canteraTPower , 1e-6), "temperature^-1");

      *canteraTPower <<=  1 / *canteraResult / *canteraResult;
      CellField& recipRecipT = cellFM.field_ref(tPowerTags[6]);
      status( field_equal(recipRecipT, *canteraTPower , 1e-6), "temperature^-2");


//      imass = massfracs.begin();
//      itempd = tVecDiff.begin();
//      i=0;
//      for(CellField::iterator icant = canteraResult->begin(); icant!=canteraResult->end(); ++itempd, ++imass, ++icant, ++i){
//        gasMix->setState_TPY( *itempd, refPressure, &(*imass)[0]);
//        double meanMW = gasMix->meanMolecularWeight();
//        gasMix->setState_UV( e0_massVec[i], Cantera::GasConstant*tVec[i]*meanMW/refPressure);
//        *icant=gasMix->temperature();
//      }
//      CellField& tempe0 = cellFM.field_ref(te0Tag);
//      status( field_equal(tempe0, *canteraResult, 1e-6), "temperature from e0");
    }
    if( status.ok() ){
      std::cout << "PASS\n";
      return true;
    }
  std::cout << "FAIL\n";
  return false;
}

int main( int iarg, char* carg[] )
{
  std::string inputFileName;
  std::string inpGroup;
  bool timings = false;
  EnergyType energyType;

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
    status( driver(  timings, H   ), "T from " + energy_name(H) );
    status( driver(  timings, E0   ), "T from " + energy_name(E0) );

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
