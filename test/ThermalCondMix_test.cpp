/*
 * ThermalConductivity_test.cpp
 *
 *  Created on: June 16, 2014
 *      Author: Nathan Yonkee
 */

#include <iostream>
#include <stdio.h>
#include <fstream>
#include "TestHelper.h"

#include <pokitt/transport/ThermalCondMix.h>
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

int main(){
  TestHelper status( true );
  try {
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

    typedef Expr::PlaceHolder   <CellField> Temp;
    typedef Expr::PlaceHolder   <CellField> MassFracs;
    typedef MixtureMolWeight    <CellField> MixtureMolWeight;
    typedef ThermalConductivity <CellField> ThermalConductivityMix;

    Expr::ExpressionFactory exprFactory;

    const Expr::Tag tTag ( "Temperature"   , Expr::STATE_NONE);
    const Expr::Tag yiTag ( "yi", Expr::STATE_NONE );
    Expr::TagList yiTags;
    for( n=0; n<nSpec; ++n ){
      std::ostringstream name;
      name << yiTag.name() << "_" << n;
      yiTags.push_back( Expr::Tag(name.str(),yiTag.context()) );
    }
    const Expr::Tag mmwTag( "mmw", Expr::STATE_NONE);
    const Expr::Tag tCondMixTag ( "Thermal Conductivity Mix", Expr::STATE_NONE);

    exprFactory.register_expression( new Temp ::Builder (tTag                 ) );
    BOOST_FOREACH( Expr::Tag yiTag, yiTags){
      exprFactory.register_expression( new MassFracs::Builder (yiTag) );
    }
    exprFactory.register_expression( new MixtureMolWeight::Builder( mmwTag, yiTag, molecularWeights));
    const Expr::ExpressionID tCondMix_id = exprFactory.register_expression( new ThermalConductivityMix::Builder (tCondMixTag ,tTag ,yiTag, mmwTag) );

    Expr::ExpressionTree tree( tCondMix_id, exprFactory, 0 );

    {
      std::ofstream fout( "ThermalConductivity.dot" );
      tree.write_tree(fout);
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
      tree.register_fields( fml );
      fml.allocate_fields( Expr::FieldAllocInfo( npts, 0, 0, false, false, false ) );
      tree.bind_fields( fml );

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

      tree.lock_fields(fml);  // prevent fields from being deallocated so that we can get them after graph execution.

      std::cout<<setup.inputFile<<" - "<<*ptit<<std::endl;
      boost::timer treetimer;
      tree.execute_tree();
      std::cout << "tree time " << treetimer.elapsed() << std::endl;

      CellField& tCondMix = cellFM.field_ref(tCondMixTag);
#     ifdef ENABLE_CUDA
      tCondMix.add_device(CPU_INDEX);
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
      SpatFldPtr<CellField> canteraResult  = SpatialFieldStore::get<CellField>(tCondMix);

      for(CellField::iterator icant = canteraResult->begin(); icant!=canteraResult->end(); ++itemp, ++imass, ++icant){
        mixTrans->thermo().setState_TPY( *itemp, refPressure, &(*imass)[0]);
        *icant=mixTrans->thermalConductivity();
      }
      status( field_equal(tCondMix, *canteraResult, 1e-12), "thermal conductivity");

    } // number of points

  } // try
  catch( Cantera::CanteraError& ){
    Cantera::showErrors();
  }

  if( status.ok() ) return 0;
    return -1;
}
