/*
 * HeatCapacity_test.cpp
 *
 *  Created on: Mar 31, 2014
 *      Author: nate
 */

#include "print.h"
#include<iostream>
#include<stdio.h>
#include<fstream>

#include <pokitt/thermo/TemperaturePowers.h>
#include <pokitt/thermo/HeatCapacity_Cp.h>

#include <cantera/IdealGasMix.h>

#include <expression/ExprLib.h>

#include <spatialops/structured/Grid.h>

#include <boost/lexical_cast.hpp>
#include <boost/timer.hpp>

#define GASCONSTANT 8314.47215

namespace SS = SpatialOps;
typedef SS::SVolField   CellField;

namespace Cantera_CXX{ class IdealGasMix; } //location of polynomial

int main(){

  std::vector<double> hctimes;
  std::vector<double> hcdiff;
  std::vector<int> ptvec;
//  ptvec.push_back(4);
  ptvec.push_back(8*8*8);
  ptvec.push_back(16*16*16);
  ptvec.push_back(32*32*32);
  ptvec.push_back(64*64*64);
#ifdef ENABLE_CUDA
  ptvec.push_back(128*128*128);
#endif

  std::ofstream myfile ("../timings.txt", std::ios::app);
  myfile<<"\nHeat Capacity times\n";
  myfile<<"ExprLib\n";
  myfile.close();
  for( std::vector<int>::iterator ptit = ptvec.begin(); ptit!= ptvec.end(); ++ptit){
    try {

      size_t i;
      size_t n;

//            const CanteraObjects::Setup setup =CanteraObjects::Setup("Mix","gri30.xml","gri30_mix");
//            const CanteraObjects::Setup setup =CanteraObjects::Setup("Mix","h2o2.xml","ohmech");
//                  const CanteraObjects::Setup setup =CanteraObjects::Setup("Mix","cp_tester.xml","const_cp");
//                      const CanteraObjects::Setup setup =CanteraObjects::Setup("Mix","cp_tester.xml","shomate_cp");
      const CanteraObjects::Setup setup =CanteraObjects::Setup("Mix","ethanol_mech.xml","gas");


      CanteraObjects& cantera = CanteraObjects::self();
      cantera.setup_cantera(setup);
      Cantera_CXX::IdealGasMix* const gasMix=cantera.get_gasmix();
      const int nSpec=gasMix->nSpecies();

      const double refPressure=gasMix->pressure();
      std::vector<double> molecularWeights(nSpec);
      gasMix->getMolecularWeights(&molecularWeights[0]);
      std::vector< double > concentrations(nSpec,1.0);

      SS::IntVec npts(*ptit,1,1);

      std::vector<double> length(3,1.0);

      typedef Expr::PlaceHolder <CellField > Temp;
      typedef Expr::PlaceHolder <CellField > MassFracs;
      typedef HeatCapacity_Cp < CellField > HeatCapacity;

      Expr::ExpressionFactory exprFactory;

      const Expr::Tag tTag ( "Temperature", Expr::STATE_NONE );

#ifdef MIX
      const Expr::Tag yiTag ( "yi", Expr::STATE_NONE );
      Expr::TagList yiTags;
      for( size_t n=0; n<nSpec; ++n ){
        std::ostringstream name;
        name << yiTag.name() << "_" << n;
        yiTags.push_back( Expr::Tag(name.str(),yiTag.context()) );
      }
      Expr::Tag hcMixTag ("hc mix", Expr::STATE_NONE);
      Expr::TagList tPowTags;
        tPowTags.push_back(Expr::Tag("t2",Expr::STATE_NONE));
        tPowTags.push_back(Expr::Tag("t3",Expr::STATE_NONE));
        tPowTags.push_back(Expr::Tag("t4",Expr::STATE_NONE));
        tPowTags.push_back(Expr::Tag("t5",Expr::STATE_NONE));
        tPowTags.push_back(Expr::Tag("tRecip",Expr::STATE_NONE));
        tPowTags.push_back(Expr::Tag("tRecipRecip",Expr::STATE_NONE));

        typedef TemperaturePowers <CellField> TemperaturePowers;
        exprFactory.register_expression( new TemperaturePowers::Builder( tPowTags,
                                                                                                         tTag));
      exprFactory.register_expression( new Temp::Builder(tTag) );
      for( size_t n=0; n<nSpec; ++n)
        exprFactory.register_expression( new MassFracs::Builder (yiTags[n]) );
      const Expr::ExpressionID hc_id = exprFactory.register_expression( new HeatCapacity::Builder( hcMixTag, tTag, tPowTags, yiTag ));
#else
      Expr::TagList hcTags;
      for( n=0; n<nSpec; ++n )
        hcTags.push_back( Expr::Tag( "hc" + boost::lexical_cast<std::string>(n), Expr::STATE_NONE ) );
      exprFactory.register_expression( new Temp::Builder(tTag) );
      const Expr::ExpressionID hc_id   = exprFactory.register_expression( new HeatCapacity   ::Builder(hcTags,tTag) );
#endif

      const SS::BoundaryCellInfo cellBCInfo = SS::BoundaryCellInfo::build<CellField>(false,false,false);
      const SS::GhostData cellGhosts(1);
      const SS::MemoryWindow vwindow( SS::get_window_with_ghost(npts,cellGhosts,cellBCInfo) );
      CellField xcoord( vwindow, cellBCInfo, cellGhosts, NULL );
      CellField ycoord( vwindow, cellBCInfo, cellGhosts, NULL );
      CellField zcoord( vwindow, cellBCInfo, cellGhosts, NULL );
      SS::Grid grid( npts, length );
      grid.set_coord<SpatialOps::XDIR>( xcoord );
      grid.set_coord<SpatialOps::YDIR>( ycoord );
      grid.set_coord<SpatialOps::ZDIR>( zcoord );
# ifdef ENABLE_CUDA
      xcoord.add_device_sync( GPU_INDEX );
      ycoord.add_device_sync( GPU_INDEX );
      zcoord.add_device_sync( GPU_INDEX );
# endif

      Expr::ExpressionTree tree( hc_id, exprFactory, 0 );
#ifdef MIX
      {
        std::ofstream fout( "HeatCapacity_mix.dot" );
        tree.write_tree(fout);
      }
#else
      {
        std::ofstream fout( "HeatCapacity.dot" );
        tree.write_tree(fout);
      }
#endif

      Expr::FieldManagerList fml;

      tree.register_fields( fml );
      //
      // allocate all fields on the patch for this tree
      //
      fml.allocate_fields( Expr::FieldAllocInfo( npts, 0, 0, false, false, false ) );
      //  fml.dump_fields(std::cout);

      //
      // bind fields to expressions
      //

      tree.bind_fields( fml );



      using namespace SpatialOps;
      CellField& temp = fml.field_manager<CellField>().field_ref(tTag);
      temp <<= 500.0 + 1000*(xcoord);

#ifdef MIX
      SpatFldPtr<CellField> sum  = SpatialFieldStore::get<CellField>(temp);
      *sum<<=0.0;
      for( n=0; n<nSpec; ++n ){
        CellField& xi = fml.field_manager<CellField>().field_ref(yiTags[n]);
        xi <<= n + 1 + xcoord;
        *sum <<= *sum + xi;
      }
      for( n=0; n<nSpec; ++n){
        CellField& xi = fml.field_manager<CellField>().field_ref(yiTags[n]);
        xi <<= xi / *sum;
      }
#endif

      std::cout<<setup.inputFile<<" - "<<*ptit<<std::endl;

      tree.lock_fields(fml);  // prevent fields from being deallocated so that we can get them after graph execution.

      tree.execute_tree();
#ifdef MIX
      CellField& hc = fml.field_manager<CellField>().field_ref(hcMixTag);

#ifdef ENABLE_CUDA
      hc.add_device_sync(CPU_INDEX);
#endif
#else
#ifdef ENABLE_CUDA
      for( n=0; n<nSpec; ++n){
        CellField& hc = fml.field_manager<CellField>().field_ref(hcTags[n]);

        hc.add_device_sync(CPU_INDEX);
      }
#endif
#endif

      std::vector<double> tVec;
      for( i=0; i<*ptit+2; ++i)
        tVec.push_back( 500.0 + 1000.0 * (i-0.5)/ *ptit);
      std::vector<double>::const_iterator itemp;
      std::vector<double>::const_iterator itend = tVec.end();

#ifdef MIX
      std::vector<std::vector<double> > molefracs;
      for( i=0; i<*ptit+2; ++i){
        std::vector<double> concentration;
        for( n=0; n<nSpec; ++n)concentration.push_back(1 + n + (i-0.5)/ *ptit);
        std::vector<double> massfrac;
        double sum = 0.0;
        for( n=0; n<nSpec; ++n) sum+=concentration[n];
        for( n=0; n<nSpec; ++n)
          massfrac.push_back(concentration[n]/sum);
        molefracs.push_back(massfrac);
      }

      std::vector< double > hc_results(*ptit+2);
      i=0;
      boost::timer hctimer;
      for( itemp = tVec.begin(); itemp!=itend; ++itemp,++i ){
        gasMix->setState_TPY(*itemp,refPressure,&molefracs[i][0]);
#ifdef CONSTVOLUME
        hc_results[i]=gasMix->cv_mass();
#else
        hc_results[i]=gasMix->cp_mass();
#endif
      }
      hctimes.push_back(hctimer.elapsed());

      double hc_maxerror = 0.0;
      i=0;
      itemp = tVec.begin();
      for( CellField::const_iterator ihc=hc.begin(); ihc!=hc.end(); ++ihc, ++i, ++itemp){
        const double diff=(*ihc-hc_results[i])/hc_results[i];
        if( fabs(diff) > hc_maxerror) hc_maxerror = fabs(diff);
        if( fabs(diff) >= 1.0e-8){
          std::cout<<"hc mix failed at T="<<*itemp<<std::endl;
          print("my hc mix",*ihc);
          print("cantera's",hc_results[i]);
          print("diff",diff);
        }
      }
      hcdiff.push_back(hc_maxerror);
#else
      std::vector< std::vector<double> > hc_results(*ptit+2);
      i=0;
      std::vector<double> hc_result(nSpec,0.0);
      boost::timer hctime;
      for( itemp = tVec.begin(); itemp!=itend; ++itemp,++i){
        gasMix->setState_TPX(*itemp,refPressure,&concentrations[0]);
        gasMix->getPartialMolarCp(&hc_result[0]);
#ifdef CONSTVOLUME
        for( std::vector<double>::iterator dit = hc_result.begin(); dit!= hc_result.end(); ++dit) *dit-=GASCONSTANT;
#endif
        hc_results[i]=hc_result;
      }
      hctimes.push_back(hctime.elapsed());

      n=0;
      double hc_maxerror = 0.0;
      for(Expr::TagList::iterator ihcs=hcTags.begin(); ihcs<hcTags.end(); ++ihcs, ++n){
        const CellField& hc = fml.field_manager<CellField>().field_ref(*ihcs);
        i=0;
        itemp = tVec.begin();
        for( CellField::const_iterator ihc=hc.begin(); ihc!=hc.end(); ++ihc, ++i, ++itemp){
          const double diff=(*ihc-hc_results[i][n])/hc_results[i][n];
          if( fabs(diff) > hc_maxerror ) hc_maxerror = fabs(diff);
          if( fabs(diff) >= 1.0e-8){
            std::cout<<"cp failed at T= "<<*itemp<<" for species n= "<<n<<std::endl;
            print("my hc",*ihc);
            print("cantera's",hc_results[i]);
            print("diff",diff);
          }
        }
      }
      hcdiff.push_back(hc_maxerror);

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
print("cantera time",hctimes);
#ifdef MIX
  print("hc mix max diff",hcdiff);
#else
  print("hc max diff",hcdiff);
#endif
  return 0;
}
