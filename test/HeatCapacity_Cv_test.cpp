/*
 * HeatCapacity_Cv_test.cpp
 *
 *  Created on: Mar 31, 2014
 *      Author: Nathan Yonkee
 */

//#define TIMINGS

#include <iostream>
#include <stdio.h>
#include <fstream>
#include "TestHelper.h"

#include <pokitt/thermo/TemperaturePowers.h>
#include <pokitt/thermo/HeatCapacity_Cv.h>

#include <expression/ExprLib.h>

#include <spatialops/structured/Grid.h>
#include <spatialops/structured/FieldComparisons.h>

#include <boost/lexical_cast.hpp>
#include <boost/timer.hpp>
#include <boost/foreach.hpp>

#include <cantera/kernel/ct_defs.h> // contains value of Cantera::GasConstant
#include <cantera/IdealGasMix.h>

namespace So = SpatialOps;
typedef So::SVolField   CellField;

namespace Cantera_CXX{ class IdealGasMix; } //location of polynomial

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

    typedef Expr::PlaceHolder      < CellField > Temp;
    typedef TemperaturePowers      < CellField > TemperaturePowers;
    typedef Expr::PlaceHolder      < CellField > MassFracs;
    typedef HeatCapacity_Cv        < CellField > HeatCapacity;
    typedef SpeciesHeatCapacity_Cv < CellField > SpeciesHeatCapacity;


    const Expr::Tag tTag ( "Temperature", Expr::STATE_NONE );
    const Expr::Tag yiTag ( "yi", Expr::STATE_NONE );
    Expr::TagList yiTags;
    for( n=0; n<nSpec; ++n ){
      std::ostringstream name;
      name << "yi_" << n;
      yiTags.push_back( Expr::Tag( name.str(), Expr::STATE_NONE ) );
    }
    Expr::Tag cvMixTag ("cv mix", Expr::STATE_NONE);
    Expr::TagList cvTags;
    for( n=0; n<nSpec; ++n ){
      cvTags.push_back( Expr::Tag( "cv" + boost::lexical_cast<std::string>(n), Expr::STATE_NONE ) );
    }

    Expr::ExpressionFactory exprFactory;

    exprFactory.register_expression( new Temp::Builder(tTag) );
    exprFactory.register_expression( new TemperaturePowers::Builder(tTag) );
    BOOST_FOREACH( Expr::Tag yiTag, yiTags){
      exprFactory.register_expression( new MassFracs::Builder (yiTag) );
    }
    const Expr::ExpressionID cvMixture_id = exprFactory.register_expression( new HeatCapacity::Builder( cvMixTag, tTag, yiTag ));
    std::set< Expr::ExpressionID > cvSpecies_id;
    for( n=0; n<nSpec; ++n ){
      cvSpecies_id.insert(exprFactory.register_expression( new SpeciesHeatCapacity::Builder(cvTags[n], tTag, n) ));
    }

    Expr::ExpressionTree mixtureTree( cvMixture_id, exprFactory, 0 );
    Expr::ExpressionTree speciesTree( cvSpecies_id, exprFactory, 0 );

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

      mixtureTree.register_fields( fml );
      speciesTree.register_fields( fml );

      fml.allocate_fields( Expr::FieldAllocInfo( npts, 0, 0, false, false, false ) );

      mixtureTree.bind_fields( fml );
      speciesTree.bind_fields( fml );

      using namespace SpatialOps;
      Expr::FieldMgrSelector<CellField>::type& cellFM = fml.field_manager< CellField>();

      CellField& temp = cellFM.field_ref(tTag);
      temp <<= 500.0 + 1000.0*xcoord;

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

      mixtureTree.lock_fields( fml );  // prevent fields from being deallocated so that we can get them after graph execution.
      speciesTree.lock_fields( fml );

#     ifdef TIMINGS
      std::cout << std::endl << setup.inputFile << " - " << *ptit << std::endl;
#     endif
      boost::timer mixTimer;
      mixtureTree.execute_tree();
#     ifdef TIMINGS
      std::cout << "PoKiTT mixture cv time  " << mixTimer.elapsed() << std::endl;
#     endif
      CellField& cv = fml.field_manager<CellField>().field_ref(cvMixTag);
#     ifdef ENABLE_CUDA
      cv.add_device(CPU_INDEX);
#     endif

      boost::timer specTimer;
      speciesTree.execute_tree();
#     ifdef TIMINGS
      std::cout << "PoKiTT species cv time  " << specTimer.elapsed() << std::endl;
#     endif

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

      std::vector< std::vector<double> >::iterator imass = massfracs.begin();
      std::vector<double>::const_iterator itemp = tVec.begin();
      SpatFldPtr<CellField> canteraResult  = SpatialFieldStore::get<CellField>(temp);
      boost::timer hMixtimer;
      for(CellField::iterator icant = canteraResult->begin(); icant!=canteraResult->end(); ++itemp, ++imass, ++icant){
        gasMix->setState_TPY( *itemp, refPressure, &(*imass)[0]);
        *icant=gasMix->cv_mass();
      }
#     ifdef TIMINGS
      std::cout << "Cantera mixture cv time " << hMixtimer.elapsed() << std::endl;
#     endif
      status( field_equal(cv, *canteraResult, 1e-14), "mix");

      std::vector< SpatFldPtr<CellField> > canteraResults;
      for( n=0; n < nSpec; ++n){
        canteraResults.push_back(SpatialFieldStore::get<CellField>(temp));
      }
      i=0;
      itemp = tVec.begin();
      imass = massfracs.begin();
      std::vector<double> cp_result(nSpec,0.0);
      boost::timer cTimer;
      for( i=0; i<*ptit+2; ++itemp, ++imass, ++i){
        gasMix->setState_TPY( *itemp, refPressure, &(*imass)[0]);
        gasMix->getPartialMolarCp(&cp_result[0]);
        for( n=0; n<nSpec; ++n){
          (*canteraResults[n])[i] = cp_result[n];
        }
      }

      for( n=0; n<nSpec; ++n){
        *canteraResults[n] <<= ( *canteraResults[n] - Cantera::GasConstant ) / molecularWeights[n];
        CellField& cv = cellFM.field_ref(cvTags[n]);
        status( field_equal(cv, *canteraResults[n], 1e-14), n);
      }

    } // number of points

    if( status.ok() ){
      std::cout << "PASS\n";
      return 0;
    }

  } // try
  catch( Cantera::CanteraError& ){
    Cantera::showErrors();
  }
  catch( std::exception& err ){
    std::cout << err.what() << std::endl;
  }

  std::cout << "FAIL\n";
  return -1;
}
