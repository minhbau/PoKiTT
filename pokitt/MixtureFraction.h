/*
 * The MIT License
 *
 * Copyright (c) 2015 The University of Utah
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifndef MIXTURE_FRACTION
#define MIXTURE_FRACTION

#include <vector>

#include <spatialops/Nebo.h>
#include <spatialops/structured/SpatialFieldStore.h>

#include <expression/Tag.h>
#include <expression/FieldRequest.h>

namespace Cantera_CXX{
  class IdealGasMix;
}

///** test the mixture fraction class.  Returns true for pass, false for fail */
//bool perform_mixfrac_tests();

namespace PoKiTT{

//--------------------------------------------------------------------
/**
 *  @class MixtureFraction
 *  @brief Provides tools related to the mixture fraction.
 *
 *  @author James C. Sutherland
 *  @date   February, 2005, vectorized April, 2015
 *
 *  The MixtureFraction class is intended for use in two-stream mixing
 *  problems.  It requires information about species molecular weights
 *  and elemental composition.  Currently this information must be provided
 *  by Cantera.  However, this could easily be abstracted into a general
 *  object of which Cantera could be a particular case...
 *
 *  The default implementation uses Bilger's mixture fraction.  This
 *  is easily changed, however, as the set_gammas() method is virtual.
 *  Deriving from this class and overloading set_gammas() will facilitate
 *  other definitions of the mixture fraction.
 *
 *  See
 *      Sutherland, Smith, & Chen
 *      "Quantification of differential diffusion in nonpremixed systems"
 *      Combustion Theory & Modelling
 *      May, 2005 Volume 9, Number 2, p. 365
 */
class MixtureFraction
{
public:
                
  /** Constructor
   *  @brief Construct a mixture fraction object
   *
   *  @param specProps : Cantera object which provides information such as
   *                     molecular and elemental weights, species atomic
   *                     composition, etc.
   *  @param oxidFrac  : vector of oxidizer mass or mole fractions
   *  @param fuelFrac  : vector of   fuel   mass or mole fractions
   *  @param inputMassFrac : flag set to true if MASS fractions are provided.
   */
  MixtureFraction( const std::vector<double> & oxidFrac,
                   const std::vector<double> & fuelFrac,
                   const bool inputMassFrac );

  MixtureFraction( Cantera_CXX::IdealGasMix& gas,
                   const std::vector<double> & oxidFrac,
                   const std::vector<double> & fuelFrac,
                   const bool inputMassFrac );
                
  /** destructor */
  ~MixtureFraction();
                
  /**
   * @brief Set the mass fractions in the fuel stream
   *  @param fuelMassFrac : Vector of mass fractions in the fuel stream.
   */
  inline void set_fuel_mass_frac( std::vector<double> & fuelMassFrac ){
    fuelMassFrac_ = fuelMassFrac;
  }
                
  /**
   *  @brief Set the mass fractions in the oxidizer stream
   *  @param oxidMassFrac : Vector of mass fractions in the oxidizer stream.
   */
  inline void set_oxid_mass_frac( std::vector<double> & oxidMassFrac ){
    oxidMassFrac_ = oxidMassFrac;
  }
                
  /**
   *  @brief Initialize everything.
   *  @param oxid : Oxidizer composition vector
   *  @param fuel : Fuel composition vector
   *  @param massFrac : true if composition is in mass fractions,
   *                    false for mole fraction.
   *  If the short constructor is used, then the user must set the
   *  oxidizer and fuel compositions and then call the initialize()
   *  method to complete initialization.
   */
  void initialize( Cantera_CXX::IdealGasMix& gas,
                   const std::vector<double> & oxid,
                   const std::vector<double> & fuel,
                   const bool massFrac );
                
  /** @brief Check to see if the mixture fraction object is ready for use. */
  inline bool is_ready() const{ return ready_; }
                
  /** @brief Return the stoichiometric mixture fraction */
  inline double stoich_mixfrac() const{ assert(ready_); return stoichMixfrac_; }
                
  /** Get fuel and oxidizer composition as double pointers */
  inline const double* fuel_massfr() const { return &(fuelMassFrac_[0]); }
  inline const double* oxid_massfr() const { return &(oxidMassFrac_[0]); }
                
  /** Get fuel and oxidizer composition as std::vector<double> */
  inline const std::vector<double> & fuel_massfr_vec() const { return fuelMassFrac_; }
  inline const std::vector<double> & oxid_massfr_vec() const { return oxidMassFrac_; }
                
  /**
   *  @brief Convert species composition to mixture fraction
   * 
   *  @param massFrac : Vector of species MASS fractions
   *  @param mixFrac  : The resulting mixture fraction
   */
  template< typename FieldT >
  inline void
  species_to_mixfrac( std::vector< boost::shared_ptr< const Expr::FieldRequest<FieldT> > >& massFrac,
                      FieldT& mixFrac ) const
  {
    using namespace SpatialOps;
    assert( ready_ );
    assert( nspec_ == (int)massFrac.size() );
    /*
     *   Z are elemental mass fractions
     *   a_{m,i} is the number of atoms m in species i.
     *   W are molecular weights
     *
     *   f = \frac{\beta-\beta_0}{\beta_1-\beta_0}
     *
     *   \beta &= \sum_{m=1}^{N_e} \gamma_m Z_m \\
     *         &= \sum_{m=1}^{N_e} \gamma_m \sum_{i=1}^{N_s}\frac{ a_{m,i}W_m Y_i }{W_i}  \\
     *         &= \sum_{i=1}^{N_s} \frac{Y_i}{W_i} \sum_{m=1}^{N_e}\gamma_m W_m a_{m,i}
     *
     *  This can be rewritten as:
     *
     *    f = \frac{-\beta_0}{\beta_1-\beta_0} + \sum_{i=1}^{N_s} \frac{Y_i}{W_i} \sum_{m=1}^{N_e}\frac{\gamma_m W_m a_{m,i}}{\beta_1-\beta_0}
     *
     *  which is what we implement below.  This results in fewer nebo assignments.
     */

    mixFrac <<= -beta0_/(beta1_-beta0_);

    for( size_t i=0; i<nspec_; ++i ){
      double fac = 0;
      for( size_t m=0; m<nelem_; ++m ){
        const int ali = nAtom_[i][m];
        if( ali > 0 ) fac += gamma_[m] * elemWt_[m] * ali;
      }
      fac /= ( beta1_ - beta0_ );

      mixFrac <<= mixFrac + fac * massFrac[i]->field_ref() / specMolWt_[i];
    }
  }

  // scalar version
  void species_to_mixfrac( const std::vector<double> & species,
                           double & mixFrac );

                
  /**
   *  @brief Convert mixture fraction to unreacted species mass fractions
   *
   *  @param mixFrac  : The mixture fraction
   *  @param massFrac : The UNREACTED mixture mass fractions obtained
   *                    by pure mixing of the fuel and oxidizer streams
   */
  template< typename FieldT >
  inline void
  mixfrac_to_species( const FieldT& mixFrac,
                      std::vector< boost::shared_ptr< const Expr::FieldRequest<FieldT> > >& massFrac ) const
  {
    using namespace SpatialOps;
    assert( ready_ );
    assert( nspec_ == (int)massFrac.size() );
    for( size_t isp=0; isp!=nspec_; ++isp ){
      massFrac[isp]->field_ref() <<= mixFrac * fuelMassFrac_[isp] + (1.0-mixFrac) * oxidMassFrac_[isp];
    }
  }

  // scalar version
  void mixfrac_to_species( const double mixFrac,
                           std::vector<double> & species ) const;
                
  /**
   *  @brief Compute equivalence ratio from mixture fraction - scalar version
   *  @param mixFrac : The mixture fraction
   */
  double mixfrac_to_equiv_ratio( const double mixFrac ) const;

  /**
   *  @brief Compute equivalence ratio from mixture fraction - vectorized version
   *  @param mixFrac : The mixture fraction
   */
  template<typename T1, typename T2>
  inline void mixfrac_to_equiv_ratio( const T1& mixFrac, T2& eqRat ) const
  {
    assert( ready_ );
    using namespace SpatialOps;
    eqRat <<= mixFrac * ( 1.0 - stoichMixfrac_ ) / ( stoichMixfrac_*( 1.0 - mixFrac ) );
  }
                
  /**
   *  @brief Compute mixture fraction from equivalence ratio
   *  @param eqRat : The equivalence ratio
   */
  double equiv_ratio_to_mixfrac( const double eqRat ) const;

  /**
   *  @brief Compute mixture fraction from equivalence ratio
   *  @param eqRat The equivalence ratio
   *  @param mixFrac the mixture fraction
   *  @tparam T1 the type for eqRat.  This could either be a double or it could be a field of the same type as mixFrac.
   *  @tparam T2 the type for mixFrac.  This is a SpatialField.
   */
  template<typename T1, typename T2>
  inline void
  equiv_ratio_to_mixfrac( const T1& eqRat, T2& mixFrac ) const
  {
    using namespace SpatialOps;
    assert( ready_ );
    mixFrac <<= ( stoichMixfrac_ * eqRat ) / ( 1.0 + stoichMixfrac_ * ( eqRat - 1.0 ) );
  }

                
  /**
   *  @brief Estimate the products of COMPLETE combustion
   *  @param mixFrac  : The mixture fraction.
   *  @param massFrac : Product composition for COMPLETE combustion
   *  @param calcMassFrac : set true to return mass fractions,
   *                        false to return mole fractions.
   *
   *  The Burke-Schumann approximation for nonpremixed combustion
   *  assumes complete and infinitely fast combustion with products
   *  of CO2 and H2O.  This method calculates the products of
   *  complete combustion.  If the composition is rich (i.e. mixture
   *  fraction greater than stoichiometric) then the product
   *  composition includes unburnt fuel, while if the composition
   *  is lean, the resulting product composition includes oxidizer.
   *  Species compositions are thus piecewise linear functions of the
   *  mixture fraction.
   */
  void estimate_product_comp( const double mixFrac,
                              std::vector<double> & massFrac,
                              const bool calcMassFrac ) const;


  /**
   *  @brief Estimate the products of COMPLETE combustion - vector form
   *  @param mixFrac  : The mixture fraction.
   *  @param massFrac : Product composition for COMPLETE combustion
   *
   *  The Burke-Schumann approximation for nonpremixed combustion
   *  assumes complete and infinitely fast combustion with products
   *  of CO2 and H2O.  This method calculates the products of
   *  complete combustion.  If the composition is rich (i.e. mixture
   *  fraction greater than stoichiometric) then the product
   *  composition includes unburnt fuel, while if the composition
   *  is lean, the resulting product composition includes oxidizer.
   *  Species compositions are thus piecewise linear functions of the
   *  mixture fraction.
   */
  template< typename FieldT >
  void estimate_product_comp( const FieldT& mixFrac,
                              std::vector<FieldT*> & massFrac ) const
  {
    using namespace SpatialOps;
    SpatFldPtr<FieldT> rich = SpatialFieldStore::get<FieldT>( mixFrac );
    SpatFldPtr<FieldT> lean = SpatialFieldStore::get<FieldT>( mixFrac );

    *rich <<= (mixFrac - stoichMixfrac_) / (1.0-stoichMixfrac_);
    *lean <<= mixFrac / stoichMixfrac_;
    for( size_t i=0; i<nspec_; ++i ){
      *massFrac[i] <<= cond( mixFrac > stoichMixfrac_,
                             stoichProdMassFrac_[i]*( 1.0 - *rich ) + fuelMassFrac_[i] * *rich )
                           ( mixFrac < stoichMixfrac_,
                               oxidMassFrac_[i] * ( 1.0 - *lean ) + stoichProdMassFrac_[i] * *lean )
                           ( stoichProdMassFrac_[i] );
    }
  }
                
  //------------------------------------------------------------------
                
protected:
                
  void set_gammas( const Cantera_CXX::IdealGasMix& gas );
                
  /** @brief Calculate the stoichiometric mixture fraction and return its value. */
  double compute_stoich_mixfrac() const;
                
  double compute_beta( const std::vector<double> & massFrac );

  template< typename FieldT >
  void compute_beta( const std::vector<const FieldT*> & massFrac, FieldT& beta ) const;
                
  void compute_elem_mass_frac( const std::vector<double> & spec,
                               std::vector<double> & elem ) const;
                
  void set_stoichiometry( const Cantera_CXX::IdealGasMix& gas );

  int nelem_, nspec_;
                
  double stoichMixfrac_;  // stoichiometric mixture fraction
  double beta0_, beta1_;  // coupling functions in fuel (1) and oxidizer (0)
                
  bool ready_;
                
  std::vector<double> gamma_;       // (nelem)       elemental weighting factors
  std::vector<double> elemMassFr_;  // work space
                
  std::vector<double> fuelMassFrac_, oxidMassFrac_;
                
  std::vector<double> specMolWt_, elemWt_;
                
  std::vector<std::vector<int> > nAtom_;

  // normalized stoichiometric coefficients of products and reactants,
  // used for estimate_product_comp()
  std::vector<double> stoichProdMassFrac_;

  size_t iCO2_, iC, iH2O_, uH_, iN2_, iN_;  // species indices
                
  //------------------------------------------------------------------
private:
                
  MixtureFraction( const MixtureFraction & );             // no copying
  MixtureFraction & operator = (const MixtureFraction& ); // no assignment
                
};

} // namespace PoKiTT

#endif
