/*
 * DiffusionCoeff_test.cpp
 *
 *  Created on: June 27, 2014
 *      Author: nate
 */

//#define TIMINGS

#include <iostream>
#include <stdio.h>
#include <fstream>
#include "TestHelper.h"

#include <pokitt/transport/DiffusionCoeffMix.h>
#include <pokitt/MixtureMolWeight.h>

#include <expression/ExprLib.h>

#include <spatialops/structured/Grid.h>
#include <spatialops/structured/FieldComparisons.h>

#include <boost/lexical_cast.hpp>
#include <boost/timer.hpp>
#include <boost/foreach.hpp>

#include <cantera/kernel/MixTransport.h>

namespace So = SpatialOps;
typedef So::SVolField   CellField;

int main()
{
  try {
    TestHelper status( true );
    const CanteraObjects::Setup setup( "Mix", "h2o2.xml",          "ohmech"    );
    //const CanteraObjects::Setup setup( "Mix", "gri30.xml",         "gri30_mix" );
    //const CanteraObjects::Setup setup( "Mix", "ethanol_mech.xml",  "gas"       );
    CanteraObjects::setup_cantera(setup);

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

    typedef Expr::PlaceHolder <CellField> Temperature;
    typedef Expr::PlaceHolder <CellField> Pressure;
    typedef Expr::PlaceHolder <CellField> MassFracs;
    typedef MixtureMolWeight  <CellField> MixtureMolWeight;
    typedef DiffusionCoeff    <CellField> DiffusionCoeffMix;

    const Expr::Tag tTag ( "Temperature"   , Expr::STATE_NONE);
    const Expr::Tag pTag ( "Pressure"   , Expr::STATE_NONE);
    const Expr::Tag yiTag ( "yi", Expr::STATE_NONE );
    Expr::TagList yiTags;
    for( n=0; n<nSpec; ++n ){
      std::ostringstream name;
      name << yiTag.name() << "_" << n;
      yiTags.push_back( Expr::Tag(name.str(),yiTag.context()) );
    }
    const Expr::Tag mmwTag( "mmw", Expr::STATE_NONE);
    Expr::TagList diffusionCoeffMixTags;
    for( n=0; n<nSpec; ++n )
      diffusionCoeffMixTags.push_back( Expr::Tag( "Di" + boost::lexical_cast<std::string>(n), Expr::STATE_NONE ) );

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

      Expr::ExpressionFactory exprFactory;

      exprFactory.register_expression( new Temperature::Builder (tTag) );
      exprFactory.register_expression( new Pressure ::Builder (pTag                 ) );
      BOOST_FOREACH( Expr::Tag yiTag, yiTags){
        exprFactory.register_expression( new MassFracs::Builder (yiTag) );
      }
      exprFactory.register_expression( new MixtureMolWeight::Builder( mmwTag, yiTag, molecularWeights));
      const Expr::ExpressionID diffCoeffMix_id = exprFactory.register_expression( new DiffusionCoeffMix::Builder (diffusionCoeffMixTags, tTag, pTag ,yiTag, mmwTag) );

      Expr::ExpressionTree tree( diffCoeffMix_id, exprFactory, 0 );

      {
        std::ofstream fout( "DiffusionCoeff.dot" );
        tree.write_tree(fout);
      }

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
      tree.register_fields( fml );
      fml.allocate_fields( Expr::FieldAllocInfo( npts, 0, 0, false, false, false ) );
      tree.bind_fields( fml );

      using namespace SpatialOps;
      Expr::FieldMgrSelector<CellField>::type& cellFM = fml.field_manager< CellField>();

      CellField& temp = cellFM.field_ref(tTag);
      temp <<= 500.0 + 1000.0*xcoord;

      CellField& p = cellFM.field_ref(pTag);
      p <<= refPressure;

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

      tree.lock_fields(fml);  // prevent fields from being deallocated so that we can get them after graph execution.

      std::cout<<setup.inputFile<<" - "<<*ptit<<std::endl;
      boost::timer treetimer;
      tree.execute_tree();
      std::cout << "tree time " << treetimer.elapsed() << std::endl;

#     ifdef ENABLE_CUDA
      for( n=0; n<nSpec; ++n){
        CellField& D = cellFM.field_ref(diffusionCoeffMixTags[n]);
        D.add_device(CPU_INDEX);
      }
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

      std::vector< SpatFldPtr<CellField> > canteraResults;
      for( n=0; n < nSpec; ++n){
        canteraResults.push_back(SpatialFieldStore::get<CellField>(xcoord));
#       ifdef ENABLE_CUDA
        canteraResults[n]->add_device( CPU_INDEX );
#       endif
      }
      i=0;
      std::vector<double>::const_iterator itemp = tVec.begin();
      std::vector< std::vector<double> >::iterator imass = massfracs.begin();
      std::vector<double> d_result(nSpec,0.0);
      boost::timer cTimer;
      for( i=0; i<*ptit+2; ++itemp, ++imass, ++i){
        mixTrans->thermo().setState_TPY( *itemp, refPressure, &(*imass)[0]);
        mixTrans->getMixDiffCoeffsMass(&d_result[0]);
        for( n=0; n<nSpec; ++n){
          (*canteraResults[n])[i] = d_result[n];
        }
      }

      for( n=0; n<nSpec; ++n ){
        CellField& D = cellFM.field_ref(diffusionCoeffMixTags[n]);
        status( field_equal(D, *canteraResults[n], 1e-12), n);
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
