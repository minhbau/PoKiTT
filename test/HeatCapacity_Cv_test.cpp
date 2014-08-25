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

#include <pokitt/thermo/TemperaturePowers.h>
#include <pokitt/thermo/HeatCapacity_Cv.h>

#include <expression/ExprLib.h>

#include <spatialops/structured/Grid.h>

#include <boost/lexical_cast.hpp>
#include <boost/timer.hpp>
#include <boost/foreach.hpp>

namespace So = SpatialOps;
typedef So::SVolField   CellField;

namespace Cantera_CXX{ class IdealGasMix; } //location of polynomial

int main()
{
  bool isFailed = false;
  try {
    const CanteraObjects::Setup setup( "Mix", "thermo_tester.xml", "const_cp"  );
  //const CanteraObjects::Setup setup( "Mix", "thermo_tester.xml", "shomate_cp");
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
  typedef Expr::PlaceHolder      < CellField > MassFracs;
  typedef TemperaturePowers      < CellField > TemperaturePowers;
  typedef HeatCapacity_Cv        < CellField > HeatCapacity;
  typedef SpeciesHeatCapacity_Cv < CellField > SpeciesHeatCapacity;

  Expr::ExpressionFactory exprFactory;

  const Expr::Tag tTag ( "Temperature", Expr::STATE_NONE );

  const Expr::Tag yiTag ( "yi", Expr::STATE_NONE );
  Expr::TagList yiTags;
  for( n=0; n<nSpec; ++n ){
    std::ostringstream name;
    name << "yi_" << n;
    yiTags.push_back( Expr::Tag( name.str(), Expr::STATE_NONE ) );
  }

  Expr::Tag hcMixTag ("hc mix", Expr::STATE_NONE);

  Expr::TagList hcTags;
  for( n=0; n<nSpec; ++n ){
    hcTags.push_back( Expr::Tag( "hc" + boost::lexical_cast<std::string>(n), Expr::STATE_NONE ) );
  }

  exprFactory.register_expression( new Temp::Builder(tTag) );
  exprFactory.register_expression( new TemperaturePowers::Builder(tTag) );
  BOOST_FOREACH( Expr::Tag yiTag, yiTags){
    exprFactory.register_expression( new MassFracs::Builder (yiTag) );
  }

  const Expr::ExpressionID hcMixture_id = exprFactory.register_expression( new HeatCapacity::Builder( hcMixTag, tTag, yiTag ));

  std::set< Expr::ExpressionID > hcSpecies_id;
  for( n=0; n<nSpec; ++n ){
    hcSpecies_id.insert(exprFactory.register_expression( new SpeciesHeatCapacity::Builder(hcTags[n], tTag, n) ));
  }

  Expr::ExpressionTree mixtureTree( hcMixture_id, exprFactory, 0 );
  Expr::ExpressionTree speciesTree( hcSpecies_id, exprFactory, 0 );

  std::vector<int> ptvec;

# ifdef TIMINGS
  ptvec.push_back(8*8*8);
  ptvec.push_back(16*16*16);
  ptvec.push_back(32*32*32);
  ptvec.push_back(64*64*64);
# ifdef ENABLE_CUDA
  ptvec.push_back(128*128*128);
# endif // ENABLE_CUDA
# else // TIMINGS
  ptvec.push_back(10);
# endif // TIMINGS

  for( std::vector<int>::iterator ptit = ptvec.begin(); ptit!= ptvec.end(); ++ptit){

      size_t i;

      So::IntVec npts(*ptit,1,1);
      std::vector<double> length(3,1.0);
      So::Grid grid( npts, length );

      const So::BoundaryCellInfo cellBCInfo = So::BoundaryCellInfo::build<CellField>(false,false,false);
      const So::GhostData cellGhosts(1);
      const So::MemoryWindow vwindow( So::get_window_with_ghost(npts,cellGhosts,cellBCInfo) );

      CellField xcoord( vwindow, cellBCInfo, cellGhosts, NULL );
      grid.set_coord<SpatialOps::XDIR>( xcoord );
#     ifdef ENABLE_CUDA
      xcoord.add_device( GPU_INDEX );
#     endif

      Expr::FieldManagerList fml;

      {
        std::ofstream mixture( "HeatCapacityMixture.dot" );
        std::ofstream species( "HeatCapacitySpecies.dot" );

        mixtureTree.write_tree( mixture );
        speciesTree.write_tree( species ) ;
      }

      mixtureTree.register_fields( fml );
      speciesTree.register_fields( fml );

      fml.allocate_fields( Expr::FieldAllocInfo( npts, 0, 0, false, false, false ) );

      mixtureTree.bind_fields( fml );
      speciesTree.bind_fields( fml );

      using namespace SpatialOps;
      CellField& temp = fml.field_manager<CellField>().field_ref(tTag);
      temp <<= 500.0 + 1000 * xcoord;

      SpatFldPtr<CellField> sum  = SpatialFieldStore::get<CellField>(temp);
      *sum<<=0.0;
      for( n=0; n<nSpec; ++n ){
        CellField& yi = fml.field_manager<CellField>().field_ref(yiTags[n]);
        yi <<= n + 1 + xcoord;
        *sum <<= *sum + yi;
      }
      for( n=0; n<nSpec; ++n){
        CellField& yi = fml.field_manager<CellField>().field_ref(yiTags[n]);
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
      std::cout << "PoKiTT mixture cp time  " << mixTimer.elapsed() << std::endl;
#     endif

      boost::timer specTimer;
      speciesTree.execute_tree();
#     ifdef TIMINGS
      std::cout << "PoKiTT species cp time  " << specTimer.elapsed() << std::endl;
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

      std::vector<double>::const_iterator itemp;
      std::vector<double>::const_iterator itend = tVec.end();
      std::vector< double > hcMix_result(*ptit+2);
      i=0;
      boost::timer hcMixtimer;
      for( itemp = tVec.begin(); itemp!=itend; ++itemp,++i ){
        gasMix->setState_TPY(*itemp,refPressure,&massfracs[i][0]);
        hcMix_result[i]=gasMix->cv_mass();
      }
#     ifdef TIMINGS
      std::cout << "Cantera mixture cp time " << hcMixtimer.elapsed() << std::endl;
#     endif

      CellField& hc = fml.field_manager<CellField>().field_ref(hcMixTag);
#     ifdef ENABLE_CUDA
      hc.add_device(CPU_INDEX);
#     endif

      i=0;
      for( CellField::const_iterator ihc=hc.begin(); ihc!=hc.end(); ++ihc, ++i){
        const double diff=(*ihc-hcMix_result[i])/hcMix_result[i];
        if( fabs(diff) >= 1.0e-14){
          isFailed = true;
        }
      }

      i=0;
      std::vector< std::vector< double > > hc_results (*ptit + 2);
      std::vector<double> hc_result(nSpec,0.0);
      boost::timer hctime;
      for( itemp = tVec.begin(); itemp!=itend; ++itemp,++i){
        gasMix->setState_TPY(*itemp,refPressure,&massfracs[i][0]);
        gasMix->getPartialMolarCp(&hc_result[0]);
        for( std::vector<double>::iterator dit = hc_result.begin(); dit!= hc_result.end(); ++dit) *dit-=Cantera::GasConstant;
        hc_results[i]=hc_result;
      }
#     ifdef TIMINGS
      std::cout << "Cantera species cp time " << hctime.elapsed() << std::endl;
#     endif

      n=0;
      for(Expr::TagList::iterator ihcs=hcTags.begin(); ihcs<hcTags.end(); ++ihcs, ++n){
        CellField& hc = fml.field_manager<CellField>().field_ref(*ihcs);
#       ifdef ENABLE_CUDA
        hc.add_device(CPU_INDEX);
#       endif
        i=0;
        for( CellField::const_iterator ihc=hc.begin(); ihc!=hc.end(); ++ihc, ++i){
          hc_results[i][n] = hc_results[i][n] / molecularWeights[n];
          const double diff=(*ihc-hc_results[i][n])/hc_results[i][n];
          if( fabs(diff) >= 1.0e-14){
            isFailed = true;
          }
        }
      }

    }

  }
  catch( Cantera::CanteraError& ){
    Cantera::showErrors();
  }

  if( isFailed ) return -1;
  return 0;
}
