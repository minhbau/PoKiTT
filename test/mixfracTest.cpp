#include <pokitt/MixtureFraction.h>

#include <cantera/IdealGasMix.h>
#include <cantera/thermo/ConstCpPoly.h>

#include <iostream>

using std::cout;
using std::endl;
using std::vector;

bool test_mixfrac()
{
  cout << "  mixture fraction object...";
  try{
    Cantera::IdealGasMix gas;

    // hard-wire some stuff for testing.  DO NOT CHANGE.
    gas.addElement("C",12.0112);
    gas.addElement("H",1.00797);
    gas.addElement("O",15.9994);
    gas.addElement("N",14.0067);

    // junk values to allow us to create thermo objects on species.
    const double tlow  = 298;
    const double thigh = 300;
    const double pref  = 101325;
    const double coefs[4] = {1,1,1,1};
    Cantera::shared_ptr<Cantera::ConstCpPoly> thermo( new Cantera::ConstCpPoly(tlow,thigh,pref,coefs) );

    typedef Cantera::Species Spec;
    typedef Cantera::shared_ptr<Spec> SpecPtr;
    Cantera::Composition comp;

    comp.clear(); comp["C"]=1; comp["H"]=4; SpecPtr ch4( new Spec("CH4", comp) );
    comp.clear(); comp["O"]=2;              SpecPtr  o2( new Spec("O2",  comp) );
    comp.clear(); comp["N"]=2;              SpecPtr  n2( new Spec("N2",  comp) );
    comp.clear(); comp["H"]=2;              SpecPtr  h2( new Spec("H2",  comp) );
    comp.clear(); comp["C"]=1; comp["O"]=2; SpecPtr co2( new Spec("CO2", comp) );
    comp.clear(); comp["H"]=2; comp["O"]=1; SpecPtr h2o( new Spec("H2O", comp) );

    ch4->thermo = thermo;
    o2 ->thermo = thermo;
    n2 ->thermo = thermo;
    h2 ->thermo = thermo;
    co2->thermo = thermo;
    h2o->thermo = thermo;

    gas.addSpecies( ch4 );
    gas.addSpecies(  o2 );
    gas.addSpecies(  n2 );
    gas.addSpecies(  h2 );
    gas.addSpecies( co2 );
    gas.addSpecies( h2o );

    // now we have set up the cantera things that are required.
    // so initialize a mixture of gases...
    const int nspec = gas.nSpecies();
    vector<double> oxid(nspec,0.0);
    vector<double> fuel(nspec,0.0);

    // mole fractions
    oxid[ gas.speciesIndex("O2") ] = 0.21;

    // mole fractions
    fuel[ gas.speciesIndex("CH4") ] = 0.221;

    PoKiTT::MixtureFraction f0( gas, oxid, fuel, false );

    // mole fractions
    oxid[ gas.speciesIndex("O2") ] = 0.21;
    oxid[ gas.speciesIndex("N2") ] = 0.79;

    // mole fractions
    fuel[ gas.speciesIndex("CH4") ] = 0.221;
    fuel[ gas.speciesIndex("H2")  ] = 0.332;
    fuel[ gas.speciesIndex("N2")  ] = 0.447;

    PoKiTT::MixtureFraction f( gas, oxid, fuel, false );

    //-----------------------------------------
    // don't change this!
    const double fst_true = 0.166925;

    bool okay = true;
    const double fst = f.stoich_mixfrac();
    if( fabs( fst - fst_true ) > 2.0e-6 )
      okay = false;

    double eqrat = f.mixfrac_to_equiv_ratio(fst);
    if( (1.0-eqrat) > 2.0e-6 )
      okay = false;

    double fst2 = f.equiv_ratio_to_mixfrac(eqrat);
    if( fabs(fst2-fst_true) > 2.0e-6 )
      okay = false;
    //-----------------------------------------

    //-----------------------------------------
    // check calculation of product composition
    // should work independently of changes above.
    vector<double> prodComp(nspec);
    double fprod;
    f.estimate_product_comp( 0.2, prodComp, true );
    f.species_to_mixfrac( prodComp, fprod );
    if( fabs(fprod-0.2) > 1.0e-6 )
      okay = false;

    f.estimate_product_comp( 0.6, prodComp, true );
    f.species_to_mixfrac( prodComp, fprod );
    if( fabs(fprod-0.6) > 1.0e-6 )
      okay = false;
    //-----------------------------------------

    if( okay )  cout << "PASS" << endl;
    else        cout << "FAIL!" << endl;

    return okay;

  }
  catch (Cantera::CanteraError) {
    Cantera::showErrors(cout);
    return false;
  }
  catch( std::exception& err ){
    std::cout << err.what() << std::endl;
    return false;
  }
  catch(...){
    std::cout << "Unhandled exception!\n";
    return false;
  }
  // should not get here.
  return false;
}

bool perform_mixfrac_tests()
{
  bool okay = true;
  okay = (test_mixfrac() == true) ? okay : false;
  return okay;
}


int main()
{
  if( perform_mixfrac_tests() ) return 0;
  return -1;
}
