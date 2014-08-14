/*
 * Enthalpy_test.cpp
 *
 *  Created on: July 7, 2014
 *      Author: nate
 */

//#define MIX
//#define FIND
#include"print.h"

#include<iostream>
#include<stdio.h>
#include<fstream>

#include "CanteraObjects.h"
#include "myEnthalpyMass.h"
#include "myEnthalpyMixMass.h"
#include "MixtureMolWeightExpr.h"

#ifdef AURORA
#include "/scratch/local/aurora/yonkee/Cantera/build/include/cantera/IdealGasMix.h"
#else
#include <cantera/IdealGasMix.h>
#endif

#include <expression/ExprLib.h>

#include <spatialops/structured/Grid.h>

#include <boost/lexical_cast.hpp>
#include <boost/timer.hpp>

#define GASCONSTANT 8314.47215

namespace SS = SpatialOps;
typedef SS::SVolField   CellField;

namespace Cantera_CXX{ class IdealGasMix; } //location of polynomial

int main(){

  std::vector<double> htimes;
  std::vector<double> hmixtimes;
  std::vector<double> hdiff;
  std::vector<double> hmixdiff;
  std::vector<int> ptvec;
  ptvec.push_back(8*8*8);
  ptvec.push_back(16*16*16);
  ptvec.push_back(32*32*32);
  ptvec.push_back(64*64*64);
#ifdef ENABLE_CUDA
  ptvec.push_back(128*128*128);
#endif

  std::ofstream myfile ("../timings.txt", std::ios::app);
  myfile<<"\Enthalpy times\n";
  myfile<<"ExprLib\n";
  myfile.close();
  for( std::vector<int>::iterator ptit = ptvec.begin(); ptit!= ptvec.end(); ++ptit){
    try {

      int i;
      int n;

      find(0);
//            const CanteraObjects::Setup setup =CanteraObjects::Setup("Mix","gri30.xml","gri30_mix");
//            const CanteraObjects::Setup setup =CanteraObjects::Setup("Mix","h2o2.xml","ohmech");
//                  const CanteraObjects::Setup setup =CanteraObjects::Setup("Mix","cp_tester.xml","const_cp");
//                      const CanteraObjects::Setup setup =CanteraObjects::Setup("Mix","cp_tester.xml","shomate_cp");
      const CanteraObjects::Setup setup =CanteraObjects::Setup("Mix","ethanol_mech.xml","gas");

      find(0.05);
      CanteraObjects& cantera = CanteraObjects::self();
      find(0.06);
      cantera.setup_cantera(setup);
      find(0.1);
      Cantera_CXX::IdealGasMix* const gasMix=cantera.get_gasmix();
      find(0.2);
      const int nSpec=gasMix->nSpecies();
      const double refPressure=gasMix->pressure();
      std::vector<double> molecularWeights(nSpec);
      gasMix->getMolecularWeights(&molecularWeights[0]);

      SS::IntVec npts(*ptit,1,1);

      std::vector<double> length(3,1.0);

      typedef Expr::PlaceHolder <CellField > Temp;
      typedef Expr::PlaceHolder <CellField > MassFracs;
      typedef SpeciesEnthalpy <CellField> Enthalpy;
      typedef EnthalpyMix < CellField > EnthalpyMix;

      find(1);

      Expr::ExpressionFactory exprFactory;
      Expr::ExpressionFactory exprFactory_mix;

      const Expr::Tag tTag ( "Temperature", Expr::STATE_NONE );
      const Expr::Tag mmwTag ( "Mixture Mol Weight", Expr::STATE_NONE );
      const Expr::Tag yiTag ( "yi", Expr::STATE_NONE );
      Expr::TagList yiTags;
      for( size_t n=0; n<nSpec; ++n ){
        std::ostringstream name;
        name << yiTag.name() << "_" << n;
        yiTags.push_back( Expr::Tag(name.str(),yiTag.context()) );
      }

      Expr::TagList hTags;
      Expr::Tag hMixTag ("enthalpy mix", Expr::STATE_NONE);
      for( n=0; n<nSpec; ++n )
        hTags.push_back( Expr::Tag( "enthalpy" + boost::lexical_cast<std::string>(n), Expr::STATE_NONE ) );


      exprFactory.register_expression( new Temp::Builder(tTag) );
      exprFactory_mix.register_expression( new Temp::Builder(tTag) );
      for( size_t n=0; n<nSpec; ++n)
        exprFactory_mix.register_expression( new MassFracs::Builder (yiTags[n]) );
      exprFactory_mix.register_expression( new MixtureMolWeightExpr::Builder ( mmwTag, yiTag, molecularWeights));

      const Expr::ExpressionID h_id   = exprFactory.register_expression( new Enthalpy   ::Builder(hTags,tTag) );
      const Expr::ExpressionID hMix_id = exprFactory_mix.register_expression( new EnthalpyMix::Builder( hMixTag, tTag, yiTag));

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
      find(2);
# ifdef ENABLE_CUDA
      xcoord.add_device( GPU_INDEX );
      ycoord.add_device( GPU_INDEX );
      zcoord.add_device( GPU_INDEX );
# endif

      find(3);
      Expr::ExpressionTree tree( h_id, exprFactory, 0 );
      Expr::ExpressionTree tree_mix( hMix_id, exprFactory_mix, 0);
      find(3.1);
      {
        std::ofstream fout( "Enthalpy.dot" );
        tree.write_tree(fout);
      }
      find(3.2);

      {
        std::ofstream fout( "EnthalpyMix.dot" );
        tree_mix.write_tree(fout);
      }

      find(3.3);
      Expr::FieldManagerList fml;
      Expr::FieldManagerList fml_mix;
      find(3.4);
      tree.register_fields( fml );
      tree_mix.register_fields( fml_mix );
      find(3);
      //
      // allocate all fields on the patch for this tree
      //
      fml.allocate_fields( Expr::FieldAllocInfo( npts, 0, 0, false, false, false ) );
      fml_mix.allocate_fields( Expr::FieldAllocInfo( npts, 0, 0, false, false, false ) );
      //  fml.dump_fields(std::cout);

      //
      // bind fields to expressions
      //

      tree.bind_fields( fml );
      tree_mix.bind_fields( fml_mix );

      find(4);


      using namespace SpatialOps;
      find(5.8);
      CellField& temp = fml.field_manager<CellField>().field_ref(tTag);
      find(5.9);
      temp <<= 500.0 + 1000*(xcoord);

      CellField& temp_mix = fml_mix.field_manager<CellField>().field_ref(tTag);
      find(5.9);
      temp_mix <<= 500.0 + 1000*(xcoord);

      SpatFldPtr<CellField> sum  = SpatialFieldStore::get<CellField>(temp_mix);
      *sum<<=0.0;
      for( n=0; n<nSpec; ++n ){
        CellField& xi = fml_mix.field_manager<CellField>().field_ref(yiTags[n]);
        xi <<= n + 1 + xcoord;
        *sum <<= *sum + xi;
      }
      for( n=0; n<nSpec; ++n){
        CellField& xi = fml_mix.field_manager<CellField>().field_ref(yiTags[n]);
        xi <<= xi / *sum;
      }

      find(6);

      std::cout<<setup.inputFile<<" - "<<*ptit<<std::endl;
find(6.1);
#ifdef MIX
      tree_mix.lock_fields(fml_mix);  // prevent fields from being deallocated so that we can get them after graph execution.

      tree_mix.execute_tree();
      find(6.4);
#else
      tree.lock_fields(fml);  // prevent fields from being deallocated so that we can get them after graph execution.

      tree.execute_tree();
      find(6.2);
#endif


      CellField& h_mix = fml_mix.field_manager<CellField>().field_ref(hMixTag);

#ifdef ENABLE_CUDA
      for( n=0; n<nSpec; ++n){
        CellField& h = fml.field_manager<CellField>().field_ref(hTags[n]);

        h.add_device_sync(CPU_INDEX);
      }

      h_mix.add_device(CPU_INDEX);
      h_mix.set_field_loc_active(CPU_INDEX);

#endif

      std::vector<double> tVec;
      for( i=0; i<*ptit+2; ++i)
        tVec.push_back( 500.0 + 1000.0 * (i-0.5)/ *ptit);
      std::vector<double>::const_iterator itemp;
      std::vector<double>::const_iterator itend = tVec.end();

#ifdef MIX
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

      std::vector< double > h_mix_results(*ptit+2);
      i=0;
      //      boost::timer stateTimer;
      //      double stateTime;
      boost::timer hmixtimer;
      for( itemp = tVec.begin(); itemp!=itend; ++itemp,++i ){
        //        double t=stateTimer.elapsed();
        gasMix->setState_TPY(*itemp,refPressure,&concentrations[i][0]);
        //        stateTime+=stateTimer.elapsed()-t;
        h_mix_results[i]=gasMix->enthalpy_mass();
      }
      hmixtimes.push_back(hmixtimer.elapsed());
      //      print("Cantera cp mix time",cpmixtimer.elapsed());
      //      print("state time",stateTime);


      double h_mix_maxerror = 0.0;
      i=0;
      itemp = tVec.begin();
      for( CellField::const_iterator ihm=h_mix.begin(); ihm!=h_mix.end(); ++ihm, ++i, ++itemp){
        const double diff=(*ihm-h_mix_results[i])/h_mix_results[i];
        if( fabs(diff) > h_mix_maxerror) h_mix_maxerror = fabs(diff);
        if( fabs(diff) >= 5.0e-8){
          std::cout<<"h mix failed at T="<<*itemp<<std::endl;
          print("my h mix",*ihm);
          print("cantera's",h_mix_results[i]);
          print("diff",diff);
        }
      }
      hmixdiff.push_back(h_mix_maxerror);
#else
      std::vector< double > concentrations(nSpec,1.0);
      std::vector< std::vector<double> > h_results(*ptit+2);
      i=0;
      std::vector<double> h_result(nSpec,0.0);
      boost::timer htimer;
      for( itemp = tVec.begin(); itemp!=itend; ++itemp,++i){
        gasMix->setState_TPX(*itemp,refPressure,&concentrations[0]);
        gasMix->getPartialMolarEnthalpies(&h_result[0]);
        h_results[i]=h_result;
      }
      htimes.push_back(htimer.elapsed());
      //      print("Cantera cp time",cptime.elapsed());

      n=0;
      double h_maxerror = 0.0;
      for(Expr::TagList::iterator ihs=hTags.begin(); ihs<hTags.end(); ++ihs, ++n){
        const CellField& h = fml.field_manager<CellField>().field_ref(*ihs);
        i=0;
        itemp = tVec.begin();
        for( CellField::const_iterator ih=h.begin(); ih!=h.end(); ++ih, ++i, ++itemp){
          const double diff=(*ih-h_results[i][n] / molecularWeights[n]) / (h_results[i][n]/molecularWeights[n]);
          if( fabs(diff) > h_maxerror ) h_maxerror = fabs(diff);
          if( fabs(diff) >= 5.0e-8){
//            std::cout<<"h failed at T="<<*itemp<<" for species n="<<n<<std::endl;
//            print("my h",*ih);
//            print("cantera's",h_results[i][n]);
//            print("diff",diff);
          }
        }
      }
      hdiff.push_back(h_maxerror);

#endif

      //      print("cp max error",cp_maxerror);
      //      print("cv max error",cv_maxerror);
      //      print("cp mix max error",cp_mix_maxerror);
      std::cout<<std::endl;

    }
    catch( Cantera::CanteraError& ){
      Cantera::showErrors();
    }
  }
#ifdef MIX
  std::ofstream myfileapp ("../timings.txt", std::ios::app);
  myfileapp<<"\nCantera\n";
  for( std::vector<double>::iterator itime = hmixtimes.begin(); itime!=hmixtimes.end(); ++itime)
    myfileapp<<*itime<<",";
  myfileapp.close();
  print("cantera h mix",hmixtimes);
  print("h mix max diff",hmixdiff);
#else
  print("cantera h",htimes);
  print("h max diff",hdiff);
#endif

  return 0;
}
