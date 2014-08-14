/*
 * Temperature_test.cpp
 *
 *      Author: nate
 */

#include <iostream>

#include <boost/timer.hpp>
#include <boost/lexical_cast.hpp>

#include <pokitt/thermo/Temperature.h>
#include <pokitt/MixtureMolWeightExpr.h>

#include <expression/ExprLib.h>

#include <spatialops/structured/Grid.h>

#define GASCONSTANT 8314.47215 //J/kmole K - value which matches Cantera results

namespace SS = SpatialOps;
typedef SS::SVolField   CellField;

namespace Cantera_CXX{ class IdealGasMix; } //location of polynomial

void print(std::string label, SS::SVolField& field);

int main(){

  std::vector<double> times;
  std::vector<double> tdiff;
  std::vector<double> tetimes;
  std::vector<double> tediff;
  std::vector<int> ptvec;
//  ptvec.push_back(5);
    ptvec.push_back(8*8*8);
    ptvec.push_back(16*16*16);
    ptvec.push_back(32*32*32);
    ptvec.push_back(64*64*64);
  #ifdef ENABLE_CUDA
    ptvec.push_back(128*128*128);
  #endif

  std::ofstream myfile ("../timings.txt", std::ios::app);
  myfile<<"\nTemperature times\n";
  myfile<<"ExprLib\n";
  myfile.close();

  for( std::vector<int>::iterator ptit = ptvec.begin(); ptit!= ptvec.end(); ++ptit){
    try{

//      const CanteraObjects::Setup setup =CanteraObjects::Setup("None","cp_tester.xml","shomate_cp");
//                        const CanteraObjects::Setup setup =CanteraObjects::Setup("None","cp_tester.xml","const_cp");
//                  const CanteraObjects::Setup setup =CanteraObjects::Setup("Mix","gri30.xml","gri30_mix");
//            const CanteraObjects::Setup setup =CanteraObjects::Setup("Mix","h2o2.xml","ohmech");
            const CanteraObjects::Setup setup =CanteraObjects::Setup("Mix","ethanol_mech.xml","gas");

      CanteraObjects& cantera = CanteraObjects::self();
      cantera.setup_cantera(setup);

      Cantera_CXX::IdealGasMix* const gasMix=cantera.get_gasmix();
      const int nSpec=gasMix->nSpecies();
      const double pressure = gasMix->pressure();
      std::vector<double> molecularWeights(nSpec);
      gasMix->getMolecularWeights(&molecularWeights[0]);

      size_t i;
      size_t n;

      SS::IntVec npts(*ptit,1,1);

      std::vector<double> length(3,1.0);

      typedef Expr::PlaceHolder <CellField> Enthalpy;
      typedef Expr::PlaceHolder <CellField> IntEnergy;
      typedef Expr::PlaceHolder <CellField> MassFracs;
      typedef Expr::PlaceHolder <CellField> KineticEnergy;
      typedef Temperature <CellField> Temperature;
      typedef TemperatureFromE0 <CellField> TemperatureE0;

      Expr::ExpressionFactory exprFactory;
      Expr::ExpressionFactory exprFactorye0;

      const Expr::Tag tTag ( "Temperature", Expr::STATE_NONE);
      const Expr::Tag te0Tag ( "Temperature E0", Expr::STATE_NONE);
      const Expr::Tag hTag ( "Enthalpy"   , Expr::STATE_NONE);
      const Expr::Tag e0Tag ( "Internal Energy"   , Expr::STATE_NONE);
      const Expr::Tag keTag ( "kinetic Energy"   , Expr::STATE_NONE);
      const Expr::Tag mmwTag ( "Mixture Mol Weight", Expr::STATE_NONE );
      const Expr::Tag yiTag ( "yi", Expr::STATE_NONE );
      Expr::TagList yiTags;
      for( size_t n=0; n<nSpec; ++n ){
        std::ostringstream name;
        name << yiTag.name() << "_" << n;
        yiTags.push_back( Expr::Tag(name.str(),yiTag.context()) );
      }

      exprFactory.register_expression( new Enthalpy::Builder (hTag) );
      for( size_t n=0; n<nSpec; ++n)
        exprFactory.register_expression( new MassFracs::Builder (yiTags[n]) );
      exprFactorye0.register_expression( new IntEnergy::Builder (e0Tag) );
      exprFactorye0.register_expression( new KineticEnergy::Builder (keTag) );
      for( size_t n=0; n<nSpec; ++n)
        exprFactorye0.register_expression( new MassFracs::Builder (yiTags[n]) );
      exprFactory.register_expression( new MixtureMolWeightExpr::Builder ( mmwTag, yiTag, molecularWeights));
      exprFactorye0.register_expression( new MixtureMolWeightExpr::Builder ( mmwTag, yiTag, molecularWeights));
      const Expr::ExpressionID temp_id = exprFactory.register_expression( new Temperature ::Builder (tTag, yiTag, hTag) );
      const Expr::ExpressionID tempe0_id = exprFactorye0.register_expression( new TemperatureE0 ::Builder (te0Tag, yiTag, e0Tag, keTag) );

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
      xcoord.add_device_sync( GPU_INDEX );
      ycoord.add_device_sync( GPU_INDEX );
      zcoord.add_device_sync( GPU_INDEX );
# endif

      Expr::ExpressionTree tree( temp_id, exprFactory, 0 );
      Expr::ExpressionTree treee0( tempe0_id, exprFactorye0, 0 );

      {
        std::ofstream fout( "Temperature.dot" );
        tree.write_tree(fout);
      }

      {
        std::ofstream fout( "Temperaturee0.dot" );
        treee0.write_tree(fout);
      }

      Expr::FieldManagerList fml;
      Expr::FieldManagerList fmle0;
      tree.register_fields( fml );
      treee0.register_fields( fmle0 );

      //
      // allocate all fields on the patch for this tree
      //
      fml.allocate_fields( Expr::FieldAllocInfo( npts, 0, 0, false, false, false ) );
      fmle0.allocate_fields( Expr::FieldAllocInfo( npts, 0, 0, false, false, false ) );

      //
      // bind fields to expressions
      //

      tree.bind_fields( fml );
      treee0.bind_fields( fmle0 );


      using namespace SpatialOps;
      Expr::FieldMgrSelector<CellField>::type& cellFM = fml.field_manager< CellField>();

      CellField& temp = cellFM.field_ref(tTag);
      CellField& enthalpy = cellFM.field_ref(hTag);

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

      temp <<= 500.0 + 1000* xcoord;

      Expr::FieldMgrSelector<CellField>::type& cellFMe0 = fmle0.field_manager< CellField>();

      CellField& tempe0 = cellFMe0.field_ref(te0Tag);

      SpatFldPtr<CellField> sume0  = SpatialFieldStore::get<CellField>(tempe0);
      *sume0 <<= 0.0;
      for( n=0; n<nSpec; ++n ){
        CellField& xi = cellFMe0.field_ref(yiTags[n]);
        xi <<= n + 1 + xcoord;
        *sume0 <<= *sume0 + xi;
      }
      for( n=0; n<nSpec; ++n){
        CellField& xi = cellFMe0.field_ref(yiTags[n]);
        xi <<= xi / *sume0;
      }

      tempe0 <<= 500.0 + 1000* xcoord;

      find(2);

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

      std::vector<double>::const_iterator itemp = tVec.begin();
      std::vector<double>::const_iterator itend = tVec.end();

      SpatFldPtr<CellField> mmw  = SpatialFieldStore::get<CellField>(temp);
      *mmw <<= 0.0;
      for( n=0; n<nSpec; ++n ){
        CellField& yi = cellFMe0.field_ref(yiTags[n]);
        *mmw <<= *mmw + yi / molecularWeights[n];
      }
*mmw <<= 1 / *mmw;

      enthalpy <<= 0.0;
      const Cantera::SpeciesThermo* spThermo = &gasMix->speciesThermo();
      for( n=0; n<nSpec; ++n){
        int type; //type of correlation used
        std::vector<double> c(15,0); //vector of Cantera's coefficients, at most 15 for NASA7
        double minTemp;
        double maxTemp;
        double refPressure;

        spThermo->reportParams(n, type, &c.at(0), minTemp, maxTemp, refPressure);

        CellField& xi = cellFM.field_ref(yiTags[n]);
        double minTempScaled = minTemp/1000;
        double maxTempScaled = maxTemp/1000;
        if( type == SIMPLE )
          enthalpy <<= enthalpy + xi*(c[1] + c[3]*(temp-c[0]))
          / molecularWeights[n];
        else if (type == NASA2)
          enthalpy<<= enthalpy + GASCONSTANT * xi* cond( temp <= c[0] && temp >= minTemp, c[ 6] + c[1] * temp + c[2]/2 * temp*temp + c[ 3]/3 * temp*temp*temp + c[ 4]/4 * temp*temp*temp*temp + c[ 5]/5 * temp*temp*temp*temp*temp) // if low temp
        ( temp >  c[0] && temp <= maxTemp, c[13] + c[8] * temp + c[9]/2 * temp*temp + c[10]/3 * temp*temp*temp + c[11]/4 * temp*temp*temp*temp+ c[12]/5 * temp*temp*temp*temp*temp)  // else if high range
        ( temp < minTemp, c[1]*temp + c[2]*minTemp*(temp-minTemp/2) + c[ 3]*minTemp*minTemp*(temp-2*minTemp/3) + c[ 4]*pow(minTemp,3)*(temp-3*minTemp/4) + c[ 5]*pow(minTemp,4)*(temp-4*minTemp/5) + c[ 6]) // else if out of bounds - low
        (                 c[8]*temp + c[9]*maxTemp*(temp-maxTemp/2) + c[10]*maxTemp*maxTemp*(temp-2*maxTemp/3) + c[11]*pow(maxTemp,3)*(temp-3*maxTemp/4) + c[12]*pow(maxTemp,4)*(temp-4*maxTemp/5) + c[13]) // else out of bounds - high
        / molecularWeights[n];
        else if ( type == SHOMATE2 )
          enthalpy<<= enthalpy + xi * 1e6 * cond( temp <= c[0] && temp >= minTemp, c[ 6] + c[1]*temp*1e-3 + c[2]/2 * temp*temp*1e-6 + c[ 3]/3 * temp*temp*temp*1e-9 + c[ 4]/4 * temp*temp*temp*temp*1e-12 - c[ 5] *1e3/temp ) // if low temp
        ( temp >  c[0] && temp <= maxTemp, c[13] + c[8]*temp*1e-3 + c[9]/2 * temp*temp*1e-6 + c[10]/3 * temp*temp*temp*1e-9 + c[11]/4 * temp*temp*temp*temp*1e-12 - c[12] * 1e3/temp )  // else if high range
        ( temp < minTemp, c[1]*temp*1e-3 + c[2] *minTempScaled* ( temp*1e-3 - minTempScaled/2 ) + c[ 3] * minTempScaled*minTempScaled * ( temp*1e-3 - 2*minTempScaled/3 ) + c[ 4] * pow(minTempScaled,3) * ( temp*1e-3 - 3*minTempScaled/4 ) - c[ 5] * pow(minTempScaled,-1) * (-temp*1e-3/minTempScaled + 2 ) + c[ 6] )
        (                 c[8]*maxTempScaled + c[9]/2 * maxTempScaled * maxTempScaled + c[10]/3 * pow(maxTempScaled,3) + c[11]/4 * pow(maxTempScaled,4) - c[12] * pow(maxTempScaled,-1) + c[13] + (temp*1e-3 - maxTempScaled)*(c[8] + c[9] * maxTempScaled + c[10] * maxTempScaled * maxTempScaled + c[11] * pow(maxTempScaled,3) + c[12] * pow(maxTempScaled,-2)))
        / molecularWeights[n];
        else
          std::cout<<"Thermo model not supported\nModel = "<<type<<std::endl;
      }

      std::vector<double> enthalpy_massVec;
      i=0;
      for( itemp = tVec.begin(); itemp != itend; ++itemp, ++i){
        gasMix->setState_TPY(*itemp,pressure, &concentrations[i][0]);
        enthalpy_massVec.push_back( gasMix -> enthalpy_mass());
      }

      CellField& e0 = cellFMe0.field_ref(e0Tag);
      e0 <<= 0.0;
      for( n=0; n<nSpec; ++n){
        int type; //type of correlation used
        std::vector<double> c(15,0); //vector of Cantera's coefficients, at most 15 for NASA7
        double minTemp;
        double maxTemp;
        double refPressure;

        spThermo->reportParams(n, type, &c.at(0), minTemp, maxTemp, refPressure);
        double minTempScaled = minTemp/1000;
        double maxTempScaled = maxTemp/1000;
        CellField& xi = cellFMe0.field_ref(yiTags[n]);
        if( type == SIMPLE )
          e0<<= e0 + xi*(c[1] + c[3]*(temp-c[0]) - GASCONSTANT*temp)
        / molecularWeights[n];
        else if (type == NASA2){
          c[1]-=1;
          c[8]-=1;

          e0<<= e0 + GASCONSTANT * xi* cond( temp <= c[0] && temp >= minTemp, c[ 6] + c[1] * temp + c[2]/2 * temp*temp + c[ 3]/3 * temp*temp*temp + c[ 4]/4 * temp*temp*temp*temp + c[ 5]/5 * temp*temp*temp*temp*temp) // if low temp
          ( temp >  c[0] && temp <= maxTemp, c[13] + c[8] * temp + c[9]/2 * temp*temp + c[10]/3 * temp*temp*temp + c[11]/4 * temp*temp*temp*temp+ c[12]/5 * temp*temp*temp*temp*temp)  // else if high range
          ( temp < minTemp, c[1]*temp + c[2]*minTemp*(temp-minTemp/2) + c[ 3]*minTemp*minTemp*(temp-2*minTemp/3) + c[ 4]*pow(minTemp,3)*(temp-3*minTemp/4) + c[ 5]*pow(minTemp,4)*(temp-4*minTemp/5) + c[ 6]) // else if out of bounds - low
          (                 c[8]*temp + c[9]*maxTemp*(temp-maxTemp/2) + c[10]*maxTemp*maxTemp*(temp-2*maxTemp/3) + c[11]*pow(maxTemp,3)*(temp-3*maxTemp/4) + c[12]*pow(maxTemp,4)*(temp-4*maxTemp/5) + c[13])
           / molecularWeights[n];
        }// else out of bounds - high
        else if ( type == SHOMATE2 )
        {
          e0<<= e0 + xi *(-GASCONSTANT*temp + 1e6 * cond( temp <= c[0] && temp >= minTemp, c[ 6] + c[1]*temp*1e-3 + c[2]/2 * temp*temp*1e-6 + c[ 3]/3 * temp*temp*temp*1e-9 + c[ 4]/4 * temp*temp*temp*temp*1e-12 - c[ 5] *1e3/temp ) // if low temp
          ( temp >  c[0] && temp <= maxTemp, c[13] + c[8]*temp*1e-3 + c[9]/2 * temp*temp*1e-6 + c[10]/3 * temp*temp*temp*1e-9 + c[11]/4 * temp*temp*temp*temp*1e-12 - c[12] * 1e3/temp )  // else if high range
          ( temp < minTemp, c[1]*temp*1e-3 + c[2] *minTempScaled* ( temp*1e-3 - minTempScaled/2 ) + c[ 3] * minTempScaled*minTempScaled * ( temp*1e-3 - 2*minTempScaled/3 ) + c[ 4] * pow(minTempScaled,3) * ( temp*1e-3 - 3*minTempScaled/4 ) - c[ 5] * pow(minTempScaled,-1) * (-temp*1e-3/minTempScaled + 2 ) + c[ 6] )
          (                 c[8]*maxTempScaled + c[9]/2 * maxTempScaled * maxTempScaled + c[10]/3 * pow(maxTempScaled,3) + c[11]/4 * pow(maxTempScaled,4) - c[12] * pow(maxTempScaled,-1) + c[13] + (temp*1e-3 - maxTempScaled)*(c[8] + c[9] * maxTempScaled + c[10] * maxTempScaled * maxTempScaled + c[11] * pow(maxTempScaled,3) + c[12] * pow(maxTempScaled,-2))))
    / molecularWeights[n];
        }
        else
          std::cout<<"Thermo model not supported\nModel = "<<type<<std::endl;
      }
      std::vector<double> e0_massVec;
      i=0;
      for( itemp = tVec.begin(); itemp != itend; ++itemp, ++i){
        gasMix->setState_TPY(*itemp,pressure, &concentrations[i][0]);
        e0_massVec.push_back( gasMix -> intEnergy_mass());
      }

      CellField& ke = cellFMe0.field_ref(keTag);

      ke <<= 0.0;

      std::cout<<setup.inputFile<<" - "<<*ptit<<std::endl;

      tempe0<<=tempe0 + 25 - 50 * xcoord;
      temp<<=temp + 25 - 50 * xcoord;
      find(5);

      std::vector<double> tVecDiff;
      for( i=0; i<*ptit+2; ++i)
        tVecDiff.push_back( 525.0 + 950.0 * (i-0.5)/ *ptit);

      find(6);
      tree.lock_fields(fml);  // prevent fields from being deallocated so that we can get them after graph execution.
      treee0.lock_fields(fmle0);
      find(7);
      tree.execute_tree();
      treee0.execute_tree();
      find(8);

#ifdef ENABLE_CUDA
      temp.add_device_sync(CPU_INDEX);
      temp.set_field_loc_active(CPU_INDEX);
#endif

      std::vector<double> results(*ptit+2,0.0);
      i=0;
      boost::timer cTimer;
      find(2);

      std::vector<double>::const_iterator itempd;
      std::vector<double>::const_iterator itendd = tVecDiff.end();
      for( itempd = tVecDiff.begin(); itempd!=itendd; ++itempd, ++i){
        gasMix->setState_TPY( *itempd, pressure, &concentrations[i][0]);
        gasMix->setState_HP( enthalpy_massVec[i], pressure);
        results[i]=gasMix->temperature();
      }
      times.push_back(cTimer.elapsed());
      //      print("cantera time",cTimer.elapsed());

      double maxerror = 0.0;
      std::vector<double>::iterator rit = results.begin();
      for ( CellField::const_iterator it = temp.begin(); it!= temp.end(); ++rit, ++it){
        const double err = (*rit-*it)/ *rit;
        if (fabs(err) > maxerror) maxerror=fabs(err);
        if(fabs(err) >= 1e-6) {
          std::cout<<"test failed @T="<<*rit<<std::endl;
          std::cout<<"my T = "<<*it<<std::endl;
        }
      }
      tdiff.push_back(maxerror);


      find(2);

      std::vector<double> resultse0(*ptit+2,0.0);
      i=0;
      boost::timer cTimere0;
      for( itempd = tVecDiff.begin(); itempd!=itendd; ++itempd, ++i){
        gasMix->setState_TPY( *itempd, pressure, &concentrations[i][0]);
        double meanMW = gasMix->meanMolecularWeight();
        gasMix->setState_UV( e0_massVec[i], GASCONSTANT*tVec[i]*meanMW/pressure);
        resultse0[i]=gasMix->temperature();
      }
      tetimes.push_back(cTimere0.elapsed());
      //      print("cantera time",cTimer.elapsed());

      maxerror = 0.0;
      rit = resultse0.begin();
      for ( CellField::const_iterator it = tempe0.begin(); it!= tempe0.end(); ++rit, ++it){
        const double err = (*rit-*it)/ *rit;
        if (fabs(err) > maxerror) maxerror=fabs(err);
        if(fabs(err) >= 1e-6) {
          std::cout<<"test failed @Te0="<<*rit<<std::endl;
          std::cout<<"my Te0 = "<<*it<<std::endl;
        }
      }
      tediff.push_back(maxerror);

      std::cout<<std::endl;

    }
    catch( Cantera::CanteraError& ){
      Cantera::showErrors();
    }
  }

  std::ofstream myfileapp ("../timings.txt", std::ios::app);
  myfileapp<<"\nCantera\n";
  for( std::vector<double>::iterator itime = times.begin(); itime!=times.end(); ++itime)
    myfileapp<<*itime<<",";
  myfileapp.close();

  print("cantera t time",times);
  print("max t diff",tdiff);
  print("cantera te0 time",tetimes);
  print("max te0 diff",tediff);
  return 0;
}

void print(std::string label, SS::SVolField& field){
  int n=0;
  for( SS::SVolField::iterator it=field.begin(); it!=field.end(); ++it,++n){
    std::cout<<label<<n<<"= "<<*it<<std::endl;
  }
}
