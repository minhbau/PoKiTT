/*
 * Temperature_test.cpp
 *
 *  Created on: August 21, 2014
 *      Author: Nathan Yonkee
 */

//#define TIMINGS

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

#include <cantera/kernel/ct_defs.h> // contains value of Cantera::GasConstant
#include <cantera/kernel/speciesThermoTypes.h> // contains definitions for which polynomial is being used
#include <cantera/IdealGasMix.h>

namespace So = SpatialOps;
typedef So::SVolField   CellField;

namespace Cantera_CXX{ class IdealGasMix; } //location of polynomial

void calculate_enthalpy( Expr::Tag hTag,
                         Expr::Tag tTag,
                         Expr::TagList yiTags,
                         Expr::FieldMgrSelector<CellField>::type& cellFM,
                         Cantera_CXX::IdealGasMix* const gasMix);

void calculate_internal_energy( Expr::Tag e0Tag,
                                Expr::Tag tTag,
                                Expr::TagList yiTags,
                                Expr::FieldMgrSelector<CellField>::type& cellFM,
                                Cantera_CXX::IdealGasMix* const gasMix);

int main()
{
  try {
    TestHelper status( true );
    const CanteraObjects::Setup setup( "Mix", "thermo_tester.xml", "const_cp"  );
    //const CanteraObjects::Setup setup( "Mix", "thermo_tester.xml", "shomate");
    //const CanteraObjects::Setup setup( "Mix", "h2o2.xml",          "ohmech"    );
    //const CanteraObjects::Setup setup( "Mix", "gri30.xml",         "gri30_mix" );
    //const CanteraObjects::Setup setup( "Mix", "ethanol_mech.xml",  "gas"       );

    CanteraObjects::setup_cantera(setup);
    Cantera_CXX::IdealGasMix* const gasMix = CanteraObjects::get_gasmix();

    const int nSpec=gasMix->nSpecies();
    size_t n;
    const double refPressure=gasMix->pressure();
    const std::vector<double>& molecularWeights = gasMix->molecularWeights();

    typedef Expr::PlaceHolder <CellField> Enthalpy;
    typedef Expr::PlaceHolder <CellField> IntEnergy;
    typedef Expr::PlaceHolder <CellField> KineticEnergy;
    typedef Expr::PlaceHolder <CellField> MassFracs;
    typedef MixtureMolWeight  <CellField> MixtureMolWeight;
    typedef Temperature       <CellField> Temperature;
    typedef TemperatureFromE0 <CellField> TemperatureE0;

    const Expr::Tag hTag ( "Enthalpy"   , Expr::STATE_NONE);
    const Expr::Tag e0Tag ( "Internal Energy"   , Expr::STATE_NONE);
    const Expr::Tag keTag ( "Kinetic Energy"   , Expr::STATE_NONE);
    const Expr::Tag yiTag ( "yi", Expr::STATE_NONE );
    Expr::TagList yiTags;
    for( size_t n=0; n<nSpec; ++n ){
      std::ostringstream name;
      name << yiTag.name() << "_" << n;
      yiTags.push_back( Expr::Tag(name.str(),yiTag.context()) );
    }
    const Expr::TagList tPowerTags = Temperature::temperature_powers_tags();
    const Expr::Tag mmwTag ( "Mixture Mol Weight", Expr::STATE_NONE );
    const Expr::Tag tTag ( "Temperature", Expr::STATE_NONE);
    const Expr::Tag te0Tag ( "Temperature E0", Expr::STATE_NONE);

    Expr::ExpressionFactory exprFactory;

    exprFactory.register_expression( new Enthalpy::Builder (hTag) );
    exprFactory.register_expression( new IntEnergy::Builder (e0Tag) );
    exprFactory.register_expression( new KineticEnergy::Builder (keTag) );
    BOOST_FOREACH( Expr::Tag yiTag, yiTags){
      exprFactory.register_expression( new MassFracs::Builder (yiTag) );
    }
    exprFactory.register_expression( new MixtureMolWeight::Builder ( mmwTag, yiTag, molecularWeights));
    const Expr::ExpressionID temp_id = exprFactory.register_expression( new Temperature ::Builder (tTag, yiTag, hTag) );
    const Expr::ExpressionID tempe0_id = exprFactory.register_expression( new TemperatureE0 ::Builder (te0Tag, yiTag, e0Tag, keTag) );

    Expr::ExpressionTree tTree( temp_id , exprFactory, 0 );
    Expr::ExpressionTree te0Tree( tempe0_id, exprFactory, 0 );

    {
      std::ofstream tGraph( "Temperature.dot" );
      std::ofstream te0Graph( "TemperatureFromE0.dot" );

      tTree.write_tree( tGraph );
      te0Tree.write_tree( te0Graph ) ;
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

    for( std::vector<int>::iterator ptit = ptvec.begin(); ptit!= ptvec.end(); ++ptit){
      size_t i;

      So::IntVec npts(*ptit,1,1);
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
//      te0Tree.register_fields( fml );

      fml.allocate_fields( Expr::FieldAllocInfo( npts, 0, 0, false, false, false ) );

      tTree.bind_fields( fml );
//      te0Tree.bind_fields( fml );

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

      std::vector<double> tVec;
      for( i=0; i<*ptit+2; ++i)
        tVec.push_back( 500.0 + 1000.0 * (i-0.5)/ *ptit);

      std::vector< std::vector<double> > massfracs;
      for( i=0; i<*ptit+2; ++i){
        std::vector<double> massfrac;
        double sum = 0.0;
        for( n=0; n<nSpec; ++n){
          massfrac.push_back(1 + n + (i-0.5)/ *ptit);
          sum+=massfrac[n];
        }
        for( n=0; n<nSpec; ++n)
          massfrac[n] = massfrac[n]/sum;
        massfracs.push_back(massfrac);
      }

      calculate_enthalpy( hTag, tTag, yiTags, cellFM, gasMix);

      std::vector<double>::const_iterator itemp;
      std::vector<double>::const_iterator itend = tVec.end();

      std::vector<double> enthalpy_massVec;
      i=0;
      for( itemp = tVec.begin(); itemp < itend; ++itemp, ++i){
        gasMix->setState_TPY(*itemp, refPressure, &massfracs[i][0]);
        enthalpy_massVec.push_back( gasMix -> enthalpy_mass());
      }

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

      temp <<= temp + 25 - 50 * xcoord;

      std::vector<double> tVecDiff;
      for( i=0; i<*ptit+2; ++i)
        tVecDiff.push_back( 525.0 + 950.0 * (i-0.5)/ *ptit);

      tTree.lock_fields( fml );  // prevent fields from being deallocated so that we can get them after graph execution.
//      te0Tree.lock_fields( fml );

#     ifdef TIMINGS
      std::cout << std::endl << setup.inputFile << " - " << *ptit << std::endl;
#     endif

      boost::timer tTimer;
      tTree.execute_tree();
#     ifdef TIMINGS
      std::cout << "PoKiTT temperature time  " << tTimer.elapsed() << std::endl;
#     endif

      boost::timer te0Timer;
//      te0Tree.execute_tree();
#     ifdef TIMINGS
      std::cout << "PoKiTT temperature from e0 time  " << te0Timer.elapsed() << std::endl;
#     endif

#     ifdef ENABLE_CUDA
      temp.add_device(CPU_INDEX);
#     endif

      std::vector< std::vector<double> >::iterator imass = massfracs.begin();
      std::vector<double>::iterator itempd = tVecDiff.begin();
      SpatFldPtr<CellField> canteraResult  = SpatialFieldStore::get<CellField>(temp);
      i=0;
      for(CellField::iterator icant = canteraResult->begin(); icant!=canteraResult->end(); ++itempd, ++imass, ++icant, ++i){
        gasMix->setState_TPY( *itempd, refPressure, &(*imass)[0]);
        gasMix->setState_HP( enthalpy_massVec[i], refPressure);
        *icant=gasMix->temperature();
      }
      status( field_equal(temp, *canteraResult, 1e-6), "temperature from h");
      *canteraResult <<= *canteraResult * *canteraResult;
      CellField& t2 = cellFM.field_ref(tPowerTags[0]);
      status( field_equal(t2, *canteraResult , 1e-6), "temperature^2");
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
      std::cout<<"here" <<std::endl;
    }
    std::cout<<"or here" <<std::endl;
    if( status.ok() ){
      std::cout << "PASS\n";
      return 0;
    }
  }
  catch( Cantera::CanteraError& ){
    Cantera::showErrors();
  }
  catch( std::exception& err ){
    std::cout << err.what() << std::endl;
  }

  std::cout << "FAIL\n";
  return -1;
}

void calculate_enthalpy( Expr::Tag hTag,
                         Expr::Tag tTag,
                         Expr::TagList yiTags,
                         Expr::FieldMgrSelector<CellField>::type& cellFM,
                         Cantera_CXX::IdealGasMix* const gasMix )
{
  const Cantera::SpeciesThermo* spThermo = &gasMix->speciesThermo();
  const int nSpec=gasMix->nSpecies();
  const std::vector<double>& molecularWeights = gasMix->molecularWeights();
  using namespace SpatialOps;
  CellField& enthalpy = cellFM.field_ref(hTag);
  CellField& temp = cellFM.field_ref(tTag);
  enthalpy <<= 0.0;
  for( size_t n=0; n<nSpec; ++n){
    int type; //type of correlation used
    std::vector<double> c(15,0); //vector of Cantera's coefficients, at most 15 for NASA7
    double minTemp;
    double maxTemp;
    double refPressure;

    spThermo->reportParams(n, type, &c[0], minTemp, maxTemp, refPressure);

    CellField& yi = cellFM.field_ref(yiTags[n]);
    double minTempScaled = minTemp/1000;
    double maxTempScaled = maxTemp/1000;
    if( type == SIMPLE )
      enthalpy <<= enthalpy + yi*(c[1] + c[3]*(temp-c[0]))
      / molecularWeights[n];
    else if (type == NASA2)
      enthalpy <<= enthalpy + Cantera::GasConstant * yi
            * cond( temp <= c[0] && temp >= minTemp, c[ 6] + c[1] * temp + c[2]/2 * temp*temp + c[ 3]/3 * temp*temp*temp + c[ 4]/4 * temp*temp*temp*temp + c[ 5]/5 * temp*temp*temp*temp*temp) // if low temp
                  ( temp >  c[0] && temp <= maxTemp, c[13] + c[8] * temp + c[9]/2 * temp*temp + c[10]/3 * temp*temp*temp + c[11]/4 * temp*temp*temp*temp+ c[12]/5 * temp*temp*temp*temp*temp)  // else if high range
                  ( temp < minTemp, c[1]*temp + c[2]*minTemp*(temp-minTemp/2) + c[ 3]*minTemp*minTemp*(temp-2*minTemp/3) + c[ 4]*pow(minTemp,3)*(temp-3*minTemp/4) + c[ 5]*pow(minTemp,4)*(temp-4*minTemp/5) + c[ 6]) // else if out of bounds - low
                  (                 c[8]*temp + c[9]*maxTemp*(temp-maxTemp/2) + c[10]*maxTemp*maxTemp*(temp-2*maxTemp/3) + c[11]*pow(maxTemp,3)*(temp-3*maxTemp/4) + c[12]*pow(maxTemp,4)*(temp-4*maxTemp/5) + c[13]) // else out of bounds - high
            / molecularWeights[n];
    else if( type == SHOMATE2 ){
      enthalpy <<= enthalpy + yi * 1e6
                * cond( temp <= c[0] && temp >= minTemp, c[ 6] + c[1]*temp*1e-3 + c[2]/2 * temp*temp*1e-6 + c[ 3]/3 * temp*temp*temp*1e-9 + c[ 4]/4 * temp*temp*temp*temp*1e-12 - c[ 5] *1e3/temp ) // if low temp
                      ( temp >  c[0] && temp <= maxTemp, c[13] + c[8]*temp*1e-3 + c[9]/2 * temp*temp*1e-6 + c[10]/3 * temp*temp*temp*1e-9 + c[11]/4 * temp*temp*temp*temp*1e-12 - c[12] * 1e3/temp )  // else if high range
                      ( temp < minTemp, c[1]*temp*1e-3 + c[2] *minTempScaled* ( temp*1e-3 - minTempScaled/2 ) + c[ 3] * minTempScaled*minTempScaled * ( temp*1e-3 - 2*minTempScaled/3 ) + c[ 4] * pow(minTempScaled,3) * ( temp*1e-3 - 3*minTempScaled/4 ) - c[ 5] * pow(minTempScaled,-1) * (-temp*1e-3/minTempScaled + 2 ) + c[ 6] )
                      (                 c[8]*maxTempScaled + c[9]/2 * maxTempScaled * maxTempScaled + c[10]/3 * pow(maxTempScaled,3) + c[11]/4 * pow(maxTempScaled,4) - c[12] * pow(maxTempScaled,-1) + c[13] + (temp*1e-3 - maxTempScaled)*(c[8] + c[9] * maxTempScaled + c[10] * maxTempScaled * maxTempScaled + c[11] * pow(maxTempScaled,3) + c[12] * pow(maxTempScaled,-2)))
                / molecularWeights[n];
    }
    else{
      std::ostringstream msg;
      msg << "Thermo model not supported\nModel = " << type << std::endl;
      throw std::invalid_argument( msg.str() );
    }
  }
}

void calculate_internal_energy( Expr::Tag e0Tag,
                                Expr::Tag tTag,
                                Expr::TagList yiTags,
                                Expr::FieldMgrSelector<CellField>::type& cellFM,
                                Cantera_CXX::IdealGasMix* const gasMix )
{
  const Cantera::SpeciesThermo* spThermo = &gasMix->speciesThermo();
  const int nSpec=gasMix->nSpecies();
  const std::vector<double>& molecularWeights = gasMix->molecularWeights();
  using namespace SpatialOps;
  CellField& e0 = cellFM.field_ref(e0Tag);
  CellField& temp = cellFM.field_ref(tTag);
  e0 <<= 0.0;
  for( size_t n=0; n<nSpec; ++n ){
    int type; //type of correlation used
    std::vector<double> c(15,0); //vector of Cantera's coefficients, at most 15 for NASA7
    double minTemp;
    double maxTemp;
    double refPressure;

    spThermo->reportParams(n, type, &c.at(0), minTemp, maxTemp, refPressure);
    double minTempScaled = minTemp/1000;
    double maxTempScaled = maxTemp/1000;
    CellField& yi = cellFM.field_ref(yiTags[n]);
    if( type == SIMPLE )
      e0<<= e0 + yi*(c[1] + c[3]*(temp-c[0]) - Cantera::GasConstant*temp)
      / molecularWeights[n];
    else if( type == NASA2 ){
      c[1]-=1;
      c[8]-=1;

      e0<<= e0 + Cantera::GasConstant * yi
            * cond( temp <= c[0] && temp >= minTemp, c[ 6] + c[1] * temp + c[2]/2 * temp*temp + c[ 3]/3 * temp*temp*temp + c[ 4]/4 * temp*temp*temp*temp + c[ 5]/5 * temp*temp*temp*temp*temp) // if low temp
                  ( temp >  c[0] && temp <= maxTemp, c[13] + c[8] * temp + c[9]/2 * temp*temp + c[10]/3 * temp*temp*temp + c[11]/4 * temp*temp*temp*temp+ c[12]/5 * temp*temp*temp*temp*temp)  // else if high range
                  ( temp < minTemp, c[1]*temp + c[2]*minTemp*(temp-minTemp/2) + c[ 3]*minTemp*minTemp*(temp-2*minTemp/3) + c[ 4]*pow(minTemp,3)*(temp-3*minTemp/4) + c[ 5]*pow(minTemp,4)*(temp-4*minTemp/5) + c[ 6]) // else if out of bounds - low
                  (                 c[8]*temp + c[9]*maxTemp*(temp-maxTemp/2) + c[10]*maxTemp*maxTemp*(temp-2*maxTemp/3) + c[11]*pow(maxTemp,3)*(temp-3*maxTemp/4) + c[12]*pow(maxTemp,4)*(temp-4*maxTemp/5) + c[13])
            / molecularWeights[n];
    }// else out of bounds - high
    else if( type == SHOMATE2 ){
      e0<<= e0 + yi *(-Cantera::GasConstant*temp + 1e6
          * cond( temp <= c[0] && temp >= minTemp, c[ 6] + c[1]*temp*1e-3 + c[2]/2 * temp*temp*1e-6 + c[ 3]/3 * temp*temp*temp*1e-9 + c[ 4]/4 * temp*temp*temp*temp*1e-12 - c[ 5] *1e3/temp ) // if low temp
                ( temp >  c[0] && temp <= maxTemp, c[13] + c[8]*temp*1e-3 + c[9]/2 * temp*temp*1e-6 + c[10]/3 * temp*temp*temp*1e-9 + c[11]/4 * temp*temp*temp*temp*1e-12 - c[12] * 1e3/temp )  // else if high range
                ( temp < minTemp, c[1]*temp*1e-3 + c[2] *minTempScaled* ( temp*1e-3 - minTempScaled/2 ) + c[ 3] * minTempScaled*minTempScaled * ( temp*1e-3 - 2*minTempScaled/3 ) + c[ 4] * pow(minTempScaled,3) * ( temp*1e-3 - 3*minTempScaled/4 ) - c[ 5] * pow(minTempScaled,-1) * (-temp*1e-3/minTempScaled + 2 ) + c[ 6] )
                (                 c[8]*maxTempScaled + c[9]/2 * maxTempScaled * maxTempScaled + c[10]/3 * pow(maxTempScaled,3) + c[11]/4 * pow(maxTempScaled,4) - c[12] * pow(maxTempScaled,-1) + c[13] + (temp*1e-3 - maxTempScaled)*(c[8] + c[9] * maxTempScaled + c[10] * maxTempScaled * maxTempScaled + c[11] * pow(maxTempScaled,3) + c[12] * pow(maxTempScaled,-2)))
          ) / molecularWeights[n];
    }
    else{
      std::ostringstream msg;
      msg << "Thermo model not supported\nModel = " << type << std::endl;
      throw std::invalid_argument( msg.str() );
    }
  }
}
