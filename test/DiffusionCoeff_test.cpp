/*
 * DiffusionCoeff_test.cpp
 *
 *  Created on: June 27, 2014
 *      Author: nate
 */

#include <iostream>
#include <stdio.h>

#include <boost/timer.hpp>

#include <pokitt/transport/DiffusionCoeffMix.h>
#include <pokitt/MixtureMolWeightExpr.h>


#include <expression/ExprLib.h>

#include <spatialops/structured/Grid.h>

namespace SS = SpatialOps;
typedef SS::SVolField   CellField;

int main(){

  std::vector<double> dtime;
  std::vector<double> ddiff;
  std::vector<int> ptvec;
  //    ptvec.push_back(3);
  ptvec.push_back(8*8*8);
  ptvec.push_back(16*16*16);
  ptvec.push_back(32*32*32);
  ptvec.push_back(64*64*64);
#ifdef ENABLE_CUDA
  ptvec.push_back(128*128*128);
#endif

  std::ofstream myfile ("../timings.txt", std::ios::app);
  myfile<<"\nDiffusion Coefficient times\n";
  myfile<<"ExprLib\n";
  myfile.close();

  for( std::vector<int>::iterator ptit = ptvec.begin(); ptit!= ptvec.end(); ++ptit){
    try{

      //      const CanteraObjects::Setup setup =CanteraObjects::Setup("None","cp_tester.xml","shomate_cp");
      //      const CanteraObjects::Setup setup =CanteraObjects::Setup("None","cp_tester.xml","const_cp");
//                const CanteraObjects::Setup setup =CanteraObjects::Setup("Mix","gri30.xml","gri30_mix");
      //      const CanteraObjects::Setup setup =CanteraObjects::Setup("Mix","air.xml","air");
//                const CanteraObjects::Setup setup =CanteraObjects::Setup("Mix","h2o2.xml","ohmech");
      const CanteraObjects::Setup setup =CanteraObjects::Setup("Mix","ethanol_mech.xml","gas");


      CanteraObjects& cantera = CanteraObjects::self();
      cantera.setup_cantera(setup);

      Cantera::Transport* transport = cantera.get_transport();
      Cantera::MixTransport* mixTrans;

      if( transport->model() ==210 | transport->model()==211)
        mixTrans = dynamic_cast<Cantera::MixTransport*>( transport );
      else {print("error, transport not mixture\ntransport model is ",transport->model()); return -1;}

      const int nSpec=mixTrans->thermo().nSpecies();
      const double pressure = mixTrans->thermo().pressure();
      std::vector<double> molecularWeights(nSpec);
      mixTrans->thermo().getMolecularWeights(&molecularWeights[0]);

      size_t i;
      size_t n;

      SS::IntVec npts(*ptit,1,1);

      std::vector<double> length(3,1.0);

      typedef Expr::PlaceHolder <CellField > Temperature;
      typedef Expr::PlaceHolder <CellField> Pressure;
      typedef Expr::PlaceHolder <CellField> MassFracs;
      typedef DiffusionCoeff <CellField> DiffusionCoeffMix;

      Expr::ExpressionFactory exprFactory;

      Expr::TagList diffusionCoeffMixTags;
      for( n=0; n<nSpec; ++n )
        diffusionCoeffMixTags.push_back( Expr::Tag( "Di" + boost::lexical_cast<std::string>(n), Expr::STATE_NONE ) );
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

      exprFactory.register_expression( new Temperature::Builder (tTag) );
      exprFactory.register_expression( new Pressure ::Builder (pTag                 ) );
      for( n=0; n<nSpec; ++n)
              exprFactory.register_expression( new MassFracs::Builder (yiTags[n]) );
            exprFactory.register_expression( new MixtureMolWeightExpr::Builder( mmwTag, yiTag, molecularWeights));
      const Expr::ExpressionID diffCoeffMix_id = exprFactory.register_expression( new DiffusionCoeffMix::Builder (diffusionCoeffMixTags, tTag, pTag ,yiTag, mmwTag) );

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

      Expr::ExpressionTree tree( diffCoeffMix_id, exprFactory, 0 );

      {
        std::ofstream fout( "DiffusionCoeff.dot" );
        tree.write_tree(fout);
      }

      Expr::FieldManagerList fml;
      tree.register_fields( fml );

      //
      // allocate all fields on the patch for this tree
      //
      fml.allocate_fields( Expr::FieldAllocInfo( npts, 0, 0, false, false, false ) );
      //fml.dump_fields( std::cout );
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

      temp <<= 500.0 + 1000.0*xcoord;

      CellField& p = cellFM.field_ref(pTag);
      p <<= mixTrans->thermo().pressure();

      tree.lock_fields(fml);  // prevent fields from being deallocated so that we can get them after graph execution.

      std::cout<<setup.inputFile<<" - "<<*ptit<<std::endl;

      tree.execute_tree();

#ifdef ENABLE_CUDA
      for( n=0; n<nSpec; ++n){
        CellField& d = fml.field_manager<CellField>().field_ref(diffusionCoeffMixTags[n]);

        d.add_device(CPU_INDEX);
      }
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

      std::vector< std::vector<double> > d_results(*ptit+2);
      i=0;
      std::vector<double>::const_iterator itemp;
      std::vector<double>::const_iterator itend = tVec.end();
      std::vector<double> d_result(nSpec,0.0);
      boost::timer cTimer;
      for(itemp = tVec.begin(); itemp!=itend; ++itemp, ++i){
        mixTrans->thermo().setState_TPY( *itemp, pressure, &concentrations[i][0]);
        mixTrans->getMixDiffCoeffsMass(&d_result[0]);
        d_results[i]=d_result;
      }
      dtime.push_back(cTimer.elapsed());
      //      print("cantera time",cTimer.elapsed());

      n=0;
      double maxerror = 0.0;
      for(Expr::TagList::iterator itag=diffusionCoeffMixTags.begin(); itag!=diffusionCoeffMixTags.end(); ++itag, ++n){
        const CellField& d = fml.field_manager<CellField>().field_ref(*itag);
        i=0;
        itemp = tVec.begin();
        for( CellField::const_iterator idiff=d.begin(); idiff!=d.end(); ++idiff, ++i, ++itemp){
          const double diff=(*idiff-d_results[i][n])/d_results[i][n];
          if( fabs(diff) > maxerror ) maxerror = fabs(diff);
//          if( std::abs(diff) >= 1.0e-8){
//            print("my value", *idiff);
//            print("Cantera's value", d_results[i][n]);
//            std::cout<<"diffusion coefficient failed for species "<<n<<" at T="<<*itemp<<std::endl;
//            return -1;
//          }
        }
      }
      ddiff.push_back(maxerror);
      //      print("maxerror",maxerror);
      std::cout<<std::endl;

    }
    catch( Cantera::CanteraError& ){
      Cantera::showErrors();
    }

  }

  std::ofstream myfileapp ("../timings.txt", std::ios::app);
  myfileapp<<"\nCantera\n";
  for( std::vector<double>::iterator itime = dtime.begin(); itime!=dtime.end(); ++itime)
    myfileapp<<*itime<<",";
  myfileapp.close();
  print("cantera time",dtime);
  print("max diff",ddiff);
  return 0;
}
