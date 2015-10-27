#include "MixtureFraction.h"
#include "CanteraObjects.h"

#include <cantera/IdealGasMix.h>

#include <cassert>
#include <string>
#include <algorithm>

using std::vector;
using std::string;

namespace PoKiTT{

//====================================================================

void mole_to_mass( const vector<double> & molWt,
                   vector<double> & moleFrac,
                   vector<double> & massFrac )
{
  const int ns = (int)moleFrac.size();
        
  assert( ns == (int)massFrac.size() );
  assert( ns == (int)molWt.size()    );
        
  double mixMW = 0.0;
  int n;
  for( n=0; n<ns; n++ )
    mixMW += moleFrac[n]*molWt[n];
  for( n=0; n<ns; n++ )
    massFrac[n] = moleFrac[n] * molWt[n] / mixMW;
}
//--------------------------------------------------------------------
void mass_to_mole( const vector<double> & molWt,
                   vector<double> & massFrac,
                   vector<double> & moleFrac )
{
  const int ns = (int)moleFrac.size();
        
  assert( ns == (int)massFrac.size() );
  assert( ns == (int)molWt.size()    );
        
  double mixMW = 0.0;
  int n;
  for( n=0; n<ns; n++ ){
    mixMW += massFrac[n]/molWt[n];
  }
  mixMW = 1.0 / mixMW;
  for( n=0; n<ns; n++ ){
    moleFrac[n] = mixMW * massFrac[n] / molWt[n];
  }
}

//====================================================================
//====================================================================
//====================================================================

MixtureFraction::MixtureFraction( const vector<double> & oxidFrac,
                                  const vector<double> & fuelFrac,
                                  const bool inputMassFrac )
{
  Cantera::IdealGasMix* gas = CanteraObjects::get_gasmix();
  initialize( *gas, oxidFrac, fuelFrac, inputMassFrac );
  CanteraObjects::restore_gasmix( gas );
}

//--------------------------------------------------------------------

MixtureFraction::MixtureFraction( Cantera::IdealGasMix& gas,
                                  const vector<double> & oxidFrac,
                                  const vector<double> & fuelFrac,
                                  const bool inputMassFrac )
{
  initialize( gas, oxidFrac, fuelFrac, inputMassFrac );
}

//--------------------------------------------------------------------

MixtureFraction::~MixtureFraction()
{}

//--------------------------------------------------------------------

void
MixtureFraction::initialize( Cantera::IdealGasMix& gas,
                             const vector<double> & oxid,
                             const vector<double> & fuel,
                             const bool massFrac )
{
  nelem_ = gas.nElements();
  nspec_ = gas.nSpecies();

  assert( nelem_ > 0 );
  assert( nspec_ > 0 );

  assert( oxid.size() == nspec_ );
  assert( fuel.size() == nspec_ );

  gamma_.resize( nelem_ );
  elemMassFr_.resize( nelem_ );

  nAtom_.clear();
  for( size_t i=0; i<nspec_; ++i ){
    nAtom_.push_back( vector<int>(nelem_) );
    for( size_t m=0; m<nelem_; ++m )
      nAtom_[i][m] = gas.nAtoms(i,m);
  }

  // ensure that mass fractions sum to unity
  double oSum=0, fSum=0;
  for( int i=0; i<nspec_; i++ ){
    oSum += oxid[i];
    fSum += fuel[i];
  }

  oxidMassFrac_.resize( nspec_ );
  fuelMassFrac_.resize( nspec_ );
  for( int i=0; i<nspec_; i++ ){
    oxidMassFrac_[i] = oxid[i] / oSum;
    fuelMassFrac_[i] = fuel[i] / fSum;
  }

  // copy species and element MW into local storage.  Unfortunate to have to copy,
  // but required since Cantera won't use std::vector<>...
        
  specMolWt_.resize(nspec_);
  for( int i=0; i<nspec_; i++ )
    specMolWt_[i] = gas.molecularWeight(i);
        
  elemWt_.resize(nelem_);
  for( int i=0; i<nelem_; i++ )
    elemWt_[i] = gas.atomicWeight(i);

  // convert to mass fractions if we got mole fractions.
  if( !massFrac ){
    mole_to_mass( specMolWt_, oxidMassFrac_, oxidMassFrac_ );
    mole_to_mass( specMolWt_, fuelMassFrac_, fuelMassFrac_ );
  }
        
  // set the elemental weighting factors
  set_gammas( gas );

  // set pure stream coupling functions
  beta0_ = compute_beta( oxidMassFrac_ );
  beta1_ = compute_beta( fuelMassFrac_ );
        
  assert( beta1_ != beta0_ );
        
  // set the stoichiometric mixture fraction
  stoichMixfrac_ = compute_stoich_mixfrac();
        
  if( (beta0_ != beta1_)      &&
      (stoichMixfrac_ >= 0.0) &&
      (stoichMixfrac_ <= 1.0) )
    ready_ = true;
  else
    ready_ = false;
        
  // set stoichiometric coefficients for complete reaction
  set_stoichiometry( gas );
}

//--------------------------------------------------------------------

void
MixtureFraction::species_to_mixfrac( const vector<double> & species,
                                     double & mixFrac )
{
  assert( ready_ );
  /*
   *               beta  - beta0        beta0 evaluated in air  stream
   *    mixfrac = ---------------       beta1 evaluated in fuel stream
   *               beta1 - beta0
   *
   *    require pure stream compositions
   *    coupling functions (gammas)
   */
  const double beta = compute_beta( species );
  mixFrac = ( beta - beta0_ ) / ( beta1_ - beta0_ );
}

//--------------------------------------------------------------------

void
MixtureFraction::mixfrac_to_species( const double mixFrac,
                                     vector<double> & species ) const
{
  assert( ready_ );
  //  if( mixFrac < 0.0 || mixFrac > 1.0 ) cout << mixFrac << endl;
  assert( mixFrac >= 0.0 && mixFrac <= 1.0 );
  assert( nspec_ == (int)species.size() );
        
  vector<double>::iterator isp;
  vector<double>::const_iterator ifuel, ioxid;
  ifuel = fuelMassFrac_.begin();
  ioxid = oxidMassFrac_.begin();
  for( isp=species.begin(); isp!=species.end(); ++isp, ++ifuel, ++ioxid ){
    *isp = (mixFrac)*(*ifuel) + (1.0-mixFrac)*(*ioxid);
    //    assert( *isp >= 0.0 && *isp <= 1.0 );
  }
}

//--------------------------------------------------------------------

double
MixtureFraction::compute_stoich_mixfrac() const
{
  return -beta0_/(beta1_-beta0_);
}

//--------------------------------------------------------------------

void
MixtureFraction::set_gammas( const Cantera::IdealGasMix& gas )
{
  // set the element name vector
  const vector<string> elemName = gas.elementNames();
        
  // Bilger's mixture fraction
  for ( int i=0; i<nelem_; i++ ){
    const double elemWt = gas.atomicWeight(i);
                
    if( elemName[i] == "C" )
      gamma_[i] = 2.0 / elemWt;
                
    else if( elemName[i] == "H" )
      gamma_[i] = 1.0 / (2.0*elemWt);
                
    else if( elemName[i] == "O" )
      gamma_[i] = -1.0 / elemWt;
                
    else
      gamma_[i] = 0.0;
  }
}

//--------------------------------------------------------------------

double
MixtureFraction::compute_beta( const vector<double> & massFrac )
{
  compute_elem_mass_frac( massFrac, elemMassFr_ );
  double beta = 0.0;
  for ( int i=0; i<nelem_; i++ ){
    beta += gamma_[i] * elemMassFr_[i];
  }
  return beta;
}

//--------------------------------------------------------------------

void
MixtureFraction::compute_elem_mass_frac( const vector<double> & spec,
                                         vector<double> & elem ) const
{
  assert( nspec_ == (int)spec.size() );
  assert( nelem_ == (int)elem.size() );
        
  for( int ielm=0; ielm<nelem_; ielm++ ){
    assert( elem[ielm] >= 0.0 );
    elem[ielm]=0.0;
    const double & eWt = elemWt_[ielm];
    for( int isp=0; isp<nspec_; isp++ ){
      //assert( spec[isp] >= 0.0 );
      elem[ielm] += nAtom_[isp][ielm] * eWt * spec[isp] / specMolWt_[isp];
    }
  }
}

//--------------------------------------------------------------------

double
MixtureFraction::mixfrac_to_equiv_ratio( const double mixFrac ) const
{
  assert( ready_ );
  assert( mixFrac < 1.0 && mixFrac >= 0.0 );
  return mixFrac*(1.0-stoichMixfrac_) / (stoichMixfrac_*(1.0-mixFrac));
}

//--------------------------------------------------------------------

double
MixtureFraction::equiv_ratio_to_mixfrac( const double eqrat ) const
{
  assert( ready_ );
  assert( eqrat >= 0.0 );
  return (stoichMixfrac_*eqrat) / (1.0+stoichMixfrac_*(eqrat-1.0));
}

//--------------------------------------------------------------------

void
MixtureFraction::estimate_product_comp( const double mixFrac,
                                        vector<double> & massFrac,
                                        const bool calcMassFrac ) const
{
  if( mixFrac > stoichMixfrac_ ){  // fuel in excess
    const double fac = (mixFrac - stoichMixfrac_) / (1.0-stoichMixfrac_);
    for( int i=0; i<nspec_; i++ ){
      massFrac[i] = stoichProdMassFrac_[i]*(1.0-fac) + fuelMassFrac_[i]*fac;
    }
  }
  else if( mixFrac < stoichMixfrac_ ){ // oxidizer in excess
    const double fac = mixFrac / stoichMixfrac_;
    for( int i=0; i<nspec_; i++ ){
      massFrac[i] = oxidMassFrac_[i]*(1.0-fac) + stoichProdMassFrac_[i]*fac;
    }
  }
  else{  // stoichiometric
    for( int i=0; i<nspec_; i++ )
      massFrac[i] = stoichProdMassFrac_[i];
  }
        
  // convert to mole fractions if requested
  if( !calcMassFrac )  mass_to_mole( specMolWt_, massFrac, massFrac );
        
  double invYsum = 1.0/accumulate( massFrac.begin(), massFrac.end(), 0.0 );
  for( int i=0; i<nspec_; i++ ) massFrac[i] *= invYsum;
  assert( invYsum < 1.001  && invYsum > 0.999 );
}
//--------------------------------------------------------------------

void
MixtureFraction::set_stoichiometry( const Cantera::IdealGasMix& gas )
{
  // set stoichiometric coefficients assuming that the products are
  //    CO2  H2O  N2  AR
  // Reactants have negative coefficients while Products have positive coefficients
        
  std::vector<double> phi_reactant;   phi_reactant.assign( nspec_, 0.0 );
  std::vector<double> phi_product ;   phi_product.assign(  nspec_, 0.0 );
        
  vector<double> elemMoles_rx( nelem_, 0.0 );
        
  //
  // set the reactant composition (mole fractions) at stoichiometric conditions
  // this is also the stoichiometric coefficients for these species.
  //
  mixfrac_to_species( stoichMixfrac_, phi_reactant );
  mass_to_mole( specMolWt_, phi_reactant, phi_reactant );
        
  //
  // get the elemental mole fractions for the reactants
  // this gives the stoichiometry for the reactants
  //
  for( int ielm=0; ielm<nelem_; ielm++ ){
    elemMoles_rx[ielm] = 0.0;
    for( int isp=0; isp<nspec_; isp++ ){
      elemMoles_rx[ielm] += nAtom_[isp][ielm] * phi_reactant[isp];
    }
  }
        
  // now we can do the elemental balances by solving a system of equations:
        
  // Carbon balance to get phi[iCO2], assuming CO2 is the only product containing C
  const int iCO2 = gas.speciesIndex("CO2");
  const int iC   = gas.elementIndex("C");
  if( iCO2 >= 0 )
    phi_product[iCO2] = elemMoles_rx[iC] + phi_reactant[iCO2];
        
  // Hydrogen balance to get phi[iH2O], assuming H2O is the only product containing H
  const int iH2O = gas.speciesIndex("H2O");
  const int iH   = gas.elementIndex("H");
  if( iH2O >= 0 )  phi_product[iH2O] = 0.5*elemMoles_rx[iH] + phi_reactant[iH2O];
        
  // N2 balance
  const int iN2 = gas.speciesIndex("N2");
  const int iN  = gas.elementIndex("N");
  if( iN2 >= 0 ) phi_product[iN2] = 0.5*elemMoles_rx[iN];
        
  // Sulfur balanceot get phi[iSO2], assuming SO2 is the only product containing S.
  const int iSspecies = gas.speciesIndex("S");
  const int iSelem    = gas.elementIndex("S");
  const int iSO2      = gas.speciesIndex("SO2");
  if( iSO2 >= 0 ) phi_product[iSO2] = phi_reactant[iSspecies] + elemMoles_rx[iSelem];
        
  // deal with other elements
  const vector<string> & elementNames = gas.elementNames();
  for( int ielm=0; ielm<nelem_; ielm++ ){
    const string & nam = elementNames[ielm];
    if( nam != "C" && nam != "H" && nam != "O" && nam != "N" && nam != "S" ){
      // see what species this element is present in
      int n = 0;
      int ispec=-1;
      for( int isp=0; isp<nspec_; isp++ ){
        if( nAtom_[isp][ielm] > 0 ){  n++;  ispec=isp; }
      }
      // don't know what to do if it is in more than one species.
      assert( n <= 1 );
      if( n == 1 ){
        assert( ispec >= 0 );
        phi_product[ispec] = elemMoles_rx[ielm];
      }
    }
  }
        
  // normalize phi_product so that we have the product mole fractions
  // at stoichiometric conditions for complete reaction.
  stoichProdMassFrac_ = phi_product;
  const double invSum = 1.0 / accumulate( stoichProdMassFrac_.begin(), stoichProdMassFrac_.end(), 0.0 );
  for( vector<double>::iterator ispmf = stoichProdMassFrac_.begin();
       ispmf != stoichProdMassFrac_.end();
       ispmf++ )
    {
      (*ispmf) *= invSum;
    }
        
}
//--------------------------------------------------------------------

} // namespace PoKiTT
