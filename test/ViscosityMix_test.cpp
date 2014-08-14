/*
 * Viscosity_test.cpp
 *
 *  Created on: June 16, 2014
 *      Author: nate
 */

//#define AURORA
//#define FINDME
#include<iostream>
#include<stdio.h>
#include<fstream>

#include <boost/timer.hpp>

#include"print.h"
#include"CanteraObjects.h"
#include"myViscosityMixMass.h"

#ifdef AURORA
#include "/scratch/local/aurora/yonkee/Cantera/build/include/cantera/IdealGasMix.h"
#else
#include<cantera/transport.h>
#endif

#include <expression/ExprLib.h>

#include <spatialops/structured/Grid.h>

namespace SS = SpatialOps;
typedef SS::SVolField   CellField;

int main(){

  std::vector<double> vtime;
  std::vector<double> vdiff;
  std::vector<int> ptvec;
  ptvec.push_back(8*8*8);
  ptvec.push_back(16*16*16);
  ptvec.push_back(32*32*32);
  ptvec.push_back(64*64*64);
#ifdef NFIELDS
#ifdef ENABLE_CUDA
  ptvec.push_back(128*128*128);
#endif
#endif

  std::ofstream myfile ("../timings.txt", std::ios::app);
  myfile<<"\nViscosity times\n";
  myfile<<"ExprLib\n";
  myfile.close();
  for( std::vector<int>::iterator ptit = ptvec.begin(); ptit!= ptvec.end(); ++ptit){
    try{

//          const CanteraObjects::Setup setup =CanteraObjects::Setup("Mix","gri30.xml","gri30_mix");
      //        const CanteraObjects::Setup setup =CanteraObjects::Setup("Mix","air.xml","air");
//          const CanteraObjects::Setup setup =CanteraObjects::Setup("Mix","h2o2.xml","ohmech");
      const CanteraObjects::Setup setup =CanteraObjects::Setup("Mix","ethanol_mech.xml","gas");

      CanteraObjects& cantera = CanteraObjects::self();
      cantera.setup_cantera(setup);

      Cantera::Transport* transport = cantera.get_transport();
      Cantera::MixTransport* mixTrans;

      if( transport->model() ==210 | transport->model()==211)
        mixTrans = dynamic_cast<Cantera::MixTransport*>( transport );
      else {print("error, transport not mixture\ntransport model is ",transport->model()); return -1;}
      const int nSpec=mixTrans->thermo().nSpecies();
      //    print("nspecies",nSpec);
      const double pressure = mixTrans->thermo().pressure();
      std::vector<double> molecularWeights(nSpec);
      mixTrans->thermo().getMolecularWeights(&molecularWeights[0]);

      size_t i;
      size_t n;

      SS::IntVec npts(*ptit,1,1);

      std::vector<double> length(3,1.0);

      typedef Expr::PlaceHolder <CellField > Temperature;
      typedef Expr::PlaceHolder <CellField> MassFracs;
      typedef Viscosity <CellField> ViscosityMix;

      Expr::ExpressionFactory exprFactory;

      const Expr::Tag visMixTag ( "Viscosity Mix", Expr::STATE_NONE);
      const Expr::Tag yiTag ( "yi", Expr::STATE_NONE );
      Expr::TagList yiTags;
      for( n=0; n<nSpec; ++n ){
        std::ostringstream name;
        name << yiTag.name() << "_" << n;
        yiTags.push_back( Expr::Tag(name.str(),yiTag.context()) );
      }
      const Expr::Tag tTag ( "Temperature"   , Expr::STATE_NONE);

      exprFactory.register_expression( new Temperature::Builder (tTag) );
      for( n=0; n<nSpec; ++n)
              exprFactory.register_expression( new MassFracs::Builder (yiTags[n]) );
            const Expr::ExpressionID visMix_id = exprFactory.register_expression( new ViscosityMix::Builder (visMixTag ,tTag ,yiTag) );

      const SS::BoundaryCellInfo cellBCInfo = SS::BoundaryCellInfo::build<CellField>(false,false,false);
      const SS::GhostData cellGhosts(1);
      const SS::MemoryWindow vwindow( SS::get_window_with_ghost(npts,cellGhosts,cellBCInfo) );
      CellField xcoord( vwindow, cellBCInfo, cellGhosts, NULL );
      CellField ycoord( vwindow, cellBCInfo, cellGhosts, NULL );
      CellField zcoord( vwindow, cellBCInfo, cellGhosts, NULL );
      SS::Grid grid( npts, length );
      grid.set_coord<SS::XDIR>( xcoord );
      grid.set_coord<SS::YDIR>( ycoord );
      grid.set_coord<SS::ZDIR>( zcoord );

# ifdef ENABLE_CUDA
      xcoord.add_device( GPU_INDEX );
      ycoord.add_device( GPU_INDEX );
      zcoord.add_device( GPU_INDEX );
# endif

      Expr::ExpressionTree tree( visMix_id, exprFactory, 0 );

      {
        std::ofstream fout( "Viscosity.dot" );
        tree.write_tree(fout);
      }

      Expr::FieldManagerList fml;
      tree.register_fields( fml );

      //
      // allocate all fields on the patch for this tree
      //
      fml.allocate_fields( Expr::FieldAllocInfo( npts, 0, 0, false, false, false ) );

      //
      // bind fields to expressions
      //

      tree.bind_fields( fml );

      using namespace SpatialOps;
      Expr::FieldMgrSelector<CellField>::type& cellFM = fml.field_manager< CellField>();

      CellField& temp = cellFM.field_ref(tTag);

      SpatFldPtr<CellField> sum  = SpatialFieldStore::get<CellField>(temp);
      *sum<<=0.0;
      for( n=0; n<nSpec; ++n ){
        CellField& xi = cellFM.field_ref(yiTags[n]);
        xi <<= n + 1 + xcoord;
        *sum <<= *sum + xi;
      }
      for( n=0; n<nSpec; ++n){
        CellField& xi = cellFM.field_ref(yiTags[n]);
        xi <<= xi / *sum;
      }

      temp <<= 500.0 + 1000.0*(xcoord);

      tree.lock_fields(fml);  // prevent fields from being deallocated so that we can get them after graph execution.

      std::cout<<setup.inputFile<<" - "<<*ptit<<std::endl;

      tree.execute_tree();

      CellField& visMix = cellFM.field_ref(visMixTag);

#ifdef ENABLE_CUDA
      visMix.add_device(CPU_INDEX);
      visMix.set_field_loc_active(CPU_INDEX);
#endif

      std::vector<std::vector<double> > concentrations;
      for( i=0; i<*ptit+2; ++i){
        std::vector<double> concentration;
        for( n=0; n<nSpec; ++n)concentration.push_back(1 + n + (i-0.5)/ *ptit);
        std::vector<double> molefrac;
        double sum = 0.0;
        for( n=0; n<nSpec; ++n) sum+=concentration[n];
        for( n=0; n<nSpec; ++n)
          molefrac.push_back(concentration[n]/sum);
        concentrations.push_back(molefrac);
      }

      std::vector<double> tVec;
      for( i=0; i<*ptit+2; ++i)
        tVec.push_back( 500.0 + 1000.0 * (i-0.5)/ *ptit);

      std::vector<double> results(*ptit+2);
      i=0;
      find(1);
      std::vector<double>::const_iterator itend = tVec.end();
      find(2);
      std::vector<double>::const_iterator itemp;
      //    boost::timer stateTimer;
      //    double stateTime = 0.0;
      //    double t;
      boost::timer cTimer;
      for(itemp = tVec.begin(); itemp!=itend; ++itemp, ++i){
        //      t = stateTimer.elapsed();
        mixTrans->thermo().setState_TPY( *itemp, pressure, &concentrations[i][0]);
        //      stateTime+=stateTimer.elapsed()-t;
        results[i]=mixTrans->viscosity();
      }
      vtime.push_back(cTimer.elapsed());
      //    print("cantera time",cTimer.elapsed());
      //    print("state time",stateTime);

      find(4);
      std::vector<double>::iterator rit = results.begin();
      itemp = tVec.begin();
      double maxerror = 0.0;
      for ( CellField::const_iterator ivis= visMix.begin(); ivis!= visMix.end(); ++rit, ++ivis, ++itemp){
        find(5);
        const double err = (*rit-*ivis)/ *rit;
        if( fabs(err) > maxerror) maxerror = fabs(err);
//        if(std::abs(err) >= 1e-4) {
//          print("test failed @ T",*itemp);
//          print("Cantera's value",*rit);
//          print("my value",*ivis);
//          print("err",err);
//        }
        find(6);
      }

      vdiff.push_back(maxerror);
      //    print("max error",maxerror);

      std::cout<<std::endl;

    }
    catch( Cantera::CanteraError& ){
      Cantera::showErrors();
    }

  }

  std::ofstream myfileapp ("../timings.txt", std::ios::app);
  myfileapp<<"\nCantera\n";
  for( std::vector<double>::iterator itime = vtime.begin(); itime!=vtime.end(); ++itime)
    myfileapp<<*itime<<",";
  myfileapp.close();
  print("cantera time",vtime);
  print("max diff",vdiff);
  return 0;
}
