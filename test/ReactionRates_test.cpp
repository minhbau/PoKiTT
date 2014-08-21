/*
 * ReactionRate_test.cpp
 *
 *  Created on: Mar 31, 2014
 *      Author: Nathan Yonkee
 */

#include <iostream>
#include <stdio.h>
#include <fstream>

#include <pokitt/kinetics/ReactionRates.h>
#include <pokitt/MixtureMolWeight.h>

#include <expression/ExprLib.h>

#include <spatialops/structured/Grid.h>

#include <boost/lexical_cast.hpp>
#include <boost/timer.hpp>
#include <boost/foreach.hpp>

#include <cantera/kernel/ct_defs.h> // contains value of Cantera::GasConstant
#include <cantera/kernel/speciesThermoTypes.h> // contains definitions for which polynomial is being used
#include <cantera/IdealGasMix.h>

namespace Cantera_CXX{ class IdealGasMix; }

namespace So = SpatialOps;
typedef So::SVolField   CellField;

int main(){
  bool isFailed = false;

  try {
    //const CanteraObjects::Setup setup( "Mix", "h2o2.xml",          "ohmech"    );
    //const CanteraObjects::Setup setup( "Mix", "gri30.xml",         "gri30_mix" );
    const CanteraObjects::Setup setup( "Mix", "ethanol_mech.xml",  "gas"       );

    CanteraObjects::setup_cantera(setup);
    Cantera_CXX::IdealGasMix* const gasMix = CanteraObjects::get_gasmix();

    const int nSpec=gasMix->nSpecies();
    size_t n;
    const double refPressure=gasMix->pressure();
    const std::vector<double>& molecularWeights = gasMix->molecularWeights();

    typedef Expr::PlaceHolder <CellField> Temp;
    typedef Expr::PlaceHolder <CellField> Pressure;
    typedef Expr::PlaceHolder <CellField> MassFracs;
    typedef ReactionRates <CellField> ReactionRate;
    typedef MixtureMolWeight <CellField> MixtureMolWeight;

    Expr::ExpressionFactory exprFactory;

    Expr::TagList rTags;
    for( n=0; n<nSpec; ++n )
      rTags.push_back( Expr::Tag( "ri" + boost::lexical_cast<std::string>(n), Expr::STATE_NONE ) );
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

    exprFactory.register_expression( new Temp ::Builder (tTag                 ) );
    exprFactory.register_expression( new Pressure ::Builder (pTag                 ) );
    for( n=0; n<nSpec; ++n)
      exprFactory.register_expression( new MassFracs::Builder (yiTags[n]) );
    exprFactory.register_expression( new MixtureMolWeight::Builder( mmwTag, yiTag, molecularWeights));

    Expr::ExpressionID rRate_id = exprFactory.register_expression( new ReactionRate::Builder (rTags, tTag, pTag, yiTag, mmwTag) );

    Expr::ExpressionTree tree( rRate_id, exprFactory, 0 );
    std::vector<int> ptvec;

#ifdef TIMINGS
    ptvec.push_back(8*8*8);
    ptvec.push_back(16*16*16);
    ptvec.push_back(32*32*32);
    ptvec.push_back(64*64*64);
#ifdef ENABLE_CUDA
    ptvec.push_back(128*128*128);
#endif
#else
    ptvec.push_back(10);
#endif

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
        std::ofstream fout( "ReactionRate.dot" );
        tree.write_tree(fout);
      }


      tree.register_fields( fml );

      fml.allocate_fields( Expr::FieldAllocInfo( npts, 0, 0, false, false, false ) );

      tree.bind_fields( fml );

      using namespace SpatialOps;
      Expr::FieldMgrSelector<CellField>::type& cellFM = fml.field_manager< CellField>();
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

      std::vector<double>::const_iterator itemp = tVec.begin();
      std::vector<double>::const_iterator itend = tVec.end();

      SpatFldPtr<CellField> mmw  = SpatialFieldStore::get<CellField>(temp);
      *mmw <<= 0.0;
      for( n=0; n<nSpec; ++n ){
        CellField& yi = cellFM.field_ref(yiTags[n]);
        *mmw <<= *mmw + yi / molecularWeights[n];
      }
      *mmw <<= 1 / *mmw;

      CellField& p = cellFM.field_ref(pTag);
      p <<= refPressure;

      tree.lock_fields(fml);  // prevent fields from being deallocated so that we can get them after graph execution.

      std::cout<<setup.inputFile<<" - "<<*ptit<<std::endl;

      boost::timer treetimer;
      tree.execute_tree();
      std::cout << "tree time " << treetimer.elapsed() << std::endl;

#ifdef ENABLE_CUDA
      for( n=0; n<nSpec; ++n){
        CellField& r = fml.field_manager<CellField>().field_ref(rTags[n]);

        r.add_device_sync(CPU_INDEX);
        r.set_field_loc_active(CPU_INDEX);
      }
#endif


      std::vector<double> nratevec(nSpec,0.0);
      std::vector<std::vector<double> > nratevecs(*ptit+2);
      i=0;
      boost::timer cTimer;
      for( itemp = tVec.begin(); itemp!=itend; ++itemp,++i){
        gasMix->setState_TPY(*itemp,refPressure,&massfracs[i][0]);
        gasMix->getNetProductionRates(&nratevec[0]);

        nratevecs[i]=nratevec;
      }

      n=0;
      for ( Expr::TagList::iterator tagit = rTags.begin(); tagit!= rTags.end(); ++tagit, ++n){
        i=0;
        CellField& r = cellFM.field_ref(*tagit);
        for ( CellField::iterator ir= r.begin(); ir!= r.end(); ++ir,++i){
          double diff=(*ir - nratevecs[i][n] * molecularWeights[n]);
          if( fabs(diff) > 1e-11) diff = diff / (nratevecs[i][n] * molecularWeights[n]);
          if( fabs(diff) > 1e-8){
            std::cout << "diff " << diff << std::endl;
            isFailed = true;}
        }
      }
      std::cout<<std::endl;
    }
  }
  catch( Cantera::CanteraError& ){
    Cantera::showErrors();
  }

  if( isFailed ) return -1;
  return 0;
}


