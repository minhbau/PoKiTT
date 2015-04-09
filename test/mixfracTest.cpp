#include <pokitt/MixtureFraction.h>

#include <cantera/IdealGasMix.h>

#include <iostream>

using std::cout;
using std::endl;
using std::vector;

bool test_mixfrac()
{
  cout << "  mixture fraction object...";
  try{

    Cantera_CXX::IdealGasMix gas;

    // hard-wire some stuff for testing.  DO NOT CHANGE.
    gas.addElement("C",12.0112);
    gas.addElement("H",1.00797);
    gas.addElement("O",15.9994);
    gas.addElement("N",14.0067);
    gas.freezeElements();

    double ch4[] = {1,4,0,0};  gas.addSpecies("CH4", ch4 );
    double o2 [] = {0,0,2,0};  gas.addSpecies("O2",  o2  );
    double n2 [] = {0,0,0,2};  gas.addSpecies("N2",  n2  );
    double h2 [] = {0,2,0,0};  gas.addSpecies("H2",  h2  );
    double co2[] = {1,0,2,0};  gas.addSpecies("CO2", co2 );
    double h2o[] = {0,2,1,0};  gas.addSpecies("H2O", h2o );
    gas.freezeSpecies();

    // now we have set up the cantera things that are required.
    // so initialize a mixture of gases...
    const int nspec = gas.nSpecies();
    vector<double> oxid(nspec);
    vector<double> fuel(nspec);

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
