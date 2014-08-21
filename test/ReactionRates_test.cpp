/*
 * ReactionRate_test.cpp
 *
 *  Created on: Mar 31, 2014
 *      Author: Nathan Yonkee
 */

#include <string>
#include <iostream>
#include <cmath>
#include <numeric>
#include <stdlib.h>

#include <boost/timer.hpp>

#include <expression/ExprLib.h>

#include <spatialops/structured/Grid.h>

#include <pokitt/kinetics/ReactionRates.h>
#include <pokitt/MixtureMolWeight.h>

#define GASCONSTANT 8314.47215

namespace Cantera_CXX{ class IdealGasMix; }

namespace SS = SpatialOps;
typedef SS::SVolField   CellField;

int main(){
bool isFailed = false;

  std::vector<double> rtimes;
  std::vector<double> rdiff;
  std::vector<int> ptvec;
      ptvec.push_back(10);

    //            const CanteraObjects::Setup setup =CanteraObjects::Setup("Mix","gri30.xml","gri30_mix");
//                const CanteraObjects::Setup setup =CanteraObjects::Setup("Mix","h2o2.xml","ohmech");
                const CanteraObjects::Setup setup =CanteraObjects::Setup("Mix","ethanol_mech.xml","gas");

           CanteraObjects::setup_cantera(setup);
           Cantera_CXX::IdealGasMix* const gasMix = CanteraObjects::get_gasmix();
          const double pressure = gasMix->pressure();
          const int nSpec=gasMix->nSpecies();
          std::vector<double> molecularWeights(nSpec);
          gasMix->getMolecularWeights(&molecularWeights[0]);

  for( std::vector<int>::iterator ptit = ptvec.begin(); ptit!= ptvec.end(); ++ptit){
    try{




      size_t i;
      size_t n;


      std::vector<double> tVec;
      for( i=0; i<*ptit+2; ++i)
        tVec.push_back( 500.0 + 1000.0 * (i-0.5)/ *ptit);

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

      std::vector<double>::const_iterator itemp;
      std::vector<double>::const_iterator itend = tVec.end();
      const std::vector<Cantera::ReactionData> m_datavec = gasMix->getReactionData();
      int NRXNS = m_datavec.size();
#ifdef CANTERATEST
      i=0;
      for( itemp = tVec.begin(); itemp!=itend; ++itemp,++i){
        gasMix->setState_TPX(*itemp,pressure,&concentrations[i][0]);
        const std::vector<Cantera::ReactionData> m_datavec = gasMix->getReactionData();
        NRXNS = m_datavec.size();
        std::vector<std::vector<double> > mcalcs = rates(gasMix);
        std::vector<double> kfwdvec=mcalcs[0];
        std::vector<double> fROPvec=mcalcs[1];
        std::vector<double> rROPvec=mcalcs[2];
        std::vector<double> nROPvec=mcalcs[3];
        std::vector<double> rateverify=mcalcs[4];

        std::vector<double> kfwd(NRXNS,0.0);
        gasMix->getFwdRateConstants(&kfwd[0]);

        std::vector<double> netROPvec(NRXNS,0.0);
        gasMix->getNetRatesOfProgress(&netROPvec[0]);

        std::vector<double> revROPvec(NRXNS,0.0);
        gasMix->getRevRatesOfProgress(&revROPvec[0]);

        std::vector<double> fwdROP(NRXNS,0.0);
        gasMix->getFwdRatesOfProgress(&fwdROP[0]);

        std::vector<double> netprodvec(nSpec,0.0);
        gasMix->getNetProductionRates(&netprodvec[0]);

        for( int r=0; r<NRXNS; ++r){
          double diff = (nROPvec[r]-netROPvec[r])/netROPvec[r];
          if(fabs(diff)>1e-8){
          }
        }
        double diffmaxv=0.0;
        for( int n=0; n<nSpec; ++n){
          double diff = (rateverify[n]-netprodvec[n])/netprodvec[n];
          if(fabs(diff) > diffmaxv) diffmaxv = fabs(diff);
          if(fabs(diff)>1e-8){
          }
        }
      }
#endif //CANTERATEST

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


      SS::IntVec npts(*ptit,1,1);

      std::vector<double> length(3,1.0);

      Expr::ExpressionFactory exprFactory;

      typedef Expr::PlaceHolder <CellField> Temp;
      typedef Expr::PlaceHolder <CellField> Pressure;
      typedef Expr::PlaceHolder <CellField> MassFracs;
      typedef ReactionRates <CellField> ReactionRate;
      typedef MixtureMolWeight <CellField> MixtureMolWeight;

      exprFactory.register_expression( new Temp ::Builder (tTag                 ) );
      exprFactory.register_expression( new Pressure ::Builder (pTag                 ) );
      for( n=0; n<nSpec; ++n)
        exprFactory.register_expression( new MassFracs::Builder (yiTags[n]) );
      exprFactory.register_expression( new MixtureMolWeight::Builder( mmwTag, yiTag, molecularWeights));

      Expr::ExpressionID rRate_id = exprFactory.register_expression( new ReactionRate::Builder (rTags, tTag, pTag, yiTag, mmwTag) );

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


      Expr::ExpressionTree tree( rRate_id, exprFactory, 0 );
      Expr::FieldManagerList fml;

      {
        std::ofstream fout( "ReactionRate.dot" );
        tree.write_tree(fout);
      }


      tree.register_fields( fml );
      //
      // allocate all fields on the patch for this tree
      //
      fml.allocate_fields( Expr::FieldAllocInfo( npts, 0, 0, false, false, false ) );

      //        fml.dump_fields(std::cout);

      //
      // bind fields to expressions
      //
      tree.bind_fields( fml );
      //

      using namespace SpatialOps;

      Expr::FieldMgrSelector<CellField>::type& cellFM = fml.field_manager< CellField>();

      CellField& t = cellFM.field_ref(tTag);

      SpatFldPtr<CellField> sum  = SpatialFieldStore::get<CellField>(t);
      *sum<<=0.0;
      for( n=0; n<nSpec; ++n ){
        CellField& xi = cellFM.field_ref(yiTags[n]);
        xi <<= n + 1 + xcoord;
        *sum <<= *sum + xi;
      }
      for( n=0; n<nSpec; ++n){
        CellField& xi = cellFM.field_ref(yiTags[n]);
        xi <<= xi / *sum;
        i=0;
      }

      t <<= 500.0 + 1000.0 * xcoord;

      CellField& p = cellFM.field_ref(pTag);
      p <<= gasMix->pressure();

      i=0;


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
        gasMix->setState_TPY(*itemp,pressure,&concentrations[i][0]);
        gasMix->getNetProductionRates(&nratevec[0]);

        nratevecs[i]=nratevec;
      }
      rtimes.push_back(cTimer.elapsed());


      double diffmax = 0.0;
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
      rdiff.push_back(diffmax);
      std::cout<<std::endl;

    }
    catch( Cantera::CanteraError& ){
      Cantera::showErrors();
    }

  }
if( isFailed ) return -1;
  return 0;
}


