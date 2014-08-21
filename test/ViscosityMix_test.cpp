/*
 * Viscosity_test.cpp
 *
 *  Created on: June 16, 2014
 *      Author: nate
 */

#include <iostream>
#include <fstream>

#include <boost/timer.hpp>

#include <pokitt/transport/ViscosityMix.h>

#include <expression/ExprLib.h>

#include <spatialops/structured/Grid.h>

namespace SS = SpatialOps;
typedef SS::SVolField   CellField;

int main(){
  bool isFailed = false;

  try {
    //const CanteraObjects::Setup setup( "Mix", "h2o2.xml",          "ohmech"    );
    //const CanteraObjects::Setup setup( "Mix", "gri30.xml",         "gri30_mix" );
    const CanteraObjects::Setup setup( "Mix", "ethanol_mech.xml",  "gas"       );
    CanteraObjects::setup_cantera(setup);

    Cantera::Transport* transport = CanteraObjects::get_transport();
    Cantera::MixTransport* mixTrans;

    if( transport->model() ==210 || transport->model()==211)
      mixTrans = dynamic_cast<Cantera::MixTransport*>( transport );
    else {
      std::cout<<"error, transport not mixture\ntransport model is " << transport->model() << std::endl;
      return -1;}

    const int nSpec=mixTrans->thermo().nSpecies();
    size_t n;
    const double refPressure=mixTrans->thermo().pressure();
    const std::vector<double>& molecularWeights = mixTrans->thermo().molecularWeights();

  std::vector<double> vtime;
  std::vector<double> vdiff;
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

  std::ofstream myfile ("../timings.txt", std::ios::app);
  myfile<<"\nViscosity times\n";
  myfile<<"ExprLib\n";
  myfile.close();
  for( std::vector<int>::iterator ptit = ptvec.begin(); ptit!= ptvec.end(); ++ptit){


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
      std::vector<double>::const_iterator itend = tVec.end();
      std::vector<double>::const_iterator itemp;
      boost::timer cTimer;
      for(itemp = tVec.begin(); itemp!=itend; ++itemp, ++i){
        mixTrans->thermo().setState_TPY( *itemp, refPressure, &concentrations[i][0]);
        results[i]=mixTrans->viscosity();
      }
      vtime.push_back(cTimer.elapsed());

      std::vector<double>::iterator rit = results.begin();
      itemp = tVec.begin();
      double maxerror = 0.0;
      for ( CellField::const_iterator ivis= visMix.begin(); ivis!= visMix.end(); ++rit, ++ivis, ++itemp){
        const double err = (*rit-*ivis)/ *rit;
        if(fabs(err) >= 1e-12) {
          std::cout << "error " << err << std::endl;
          isFailed = true;
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
