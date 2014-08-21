#ifndef HeatCapacity_Cp_Expr_h
#define HeatCapacity_Cp_Expr_h

#include <expression/Expression.h>

#include <pokitt/CanteraObjects.h> //include cantera wrapper
#include <pokitt/thermo/TemperaturePowers.h>

#include <cantera/kernel/ct_defs.h> // contains value of Cantera::GasConstant
#include <cantera/kernel/speciesThermoTypes.h> // contains definitions for which polynomial is being used
#include <cantera/IdealGasMix.h>

namespace Cantera_CXX{ class IdealGasMix; } // location of polynomial coefficients

/**
 *  \class  HeatCapacity_Cp
 *  \author Nate Yonkee
 *  \date July, 2014
 *
 *  \brief Calculates the constant pressure heat capacity of a mixture
 *  using either NASA7 or Shomate polynomials with 2 temperature ranges.
 *  Units of J/kg/K
 *
 *  Five coefficients \f$(a_0,\dots,a_4)\f$ are used to represent
 *  \f$ c_p^0(T)\f$ of species i as a polynomial in \f$ T \f$
 *
 * NASA7:
 *
 * \f[
 * \frac{c_{p,i}(T)}{R} = a_0 + a_1 T + a_2 T^2 + a_3 T^3 + a_4 T^4
 * \f]
 *
 * Shomate:
 *
 * \f[
 * \frac{c_{p,i}(T)}{1000} = a_0 + a_1 t + a_2 t^2 + a_3 t^3 + \frac{a_4}{t^2}
 * \f]
 *
 * Where \f[ t = T (K)/1000 \f]
 *
 * The heat capacity is a weighted average of the pure species \f$ c_p^0(T)\f$,
 *
 * \f[
 * c_p(T) = \sum_{n=0}^{nSpec} Y_n c_{p,n}
 * \f]
 */

template< typename FieldT >
class HeatCapacity_Cp
 : public Expr::Expression<FieldT>
{
  const Expr::Tag tTag_;
  const Expr::TagList tPowerTags_;
  Expr::TagList massFracTags_;
  const FieldT* t_;
  std::vector<const FieldT*> tPowers_;
  std::vector<const FieldT*> massFracs_;
  Cantera_CXX::IdealGasMix* const gasMix_;
  const int nSpec_; // number of species

  HeatCapacity_Cp( const Expr::Tag& tTag,
                   const Expr::Tag& massFracTag );

public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a HeatCapacity_Cp expression
     *  @param resultTag the mixture heat capacity
     *  @param tTag tag for temperature
     *  @param massFracTag tag for mass fractions of each species, ordering is consistent with Cantera input
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::Tag& tTag,
             const Expr::Tag& massFracTag );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag tTag_;
    const Expr::Tag massFracTag_;
  };

  ~HeatCapacity_Cp();
  void advertise_dependents( Expr::ExprDeps& exprDeps );
  void bind_fields( const Expr::FieldManagerList& fml );
  void evaluate();
};

/**
 *  \class  SpeciesHeatCapacity_Cp
 *  \author Nate Yonkee
 *  \date May, 2014
 *
 *  \brief Calculates the constant pressure heat capacity of a single species
 *  using either NASA7 or Shomate polynomials with 2 temperature ranges.
 *  Units of J/kg/K
 *
 *  Five coefficients \f$(a_0,\dots,a_4)\f$ are used to represent
 * \f$ c_p^0(T)\f$ as a polynomial in \f$ T \f$
 *
 * NASA7:
 *
 * \f[
 * \frac{c_p(T)}{R} = a_0 + a_1 T + a_2 T^2 + a_3 T^3 + a_4 T^4
 * \f]
 *
 * Shomate:
 *
 * \f[
 * \frac{c_p(T)}{1000} = a_0 + a_1 t + a_2 t^2 + a_3 t^3 + \frac{a_4}{t^2}
 * \f]
 *
 * Where \f[ t = T (K)/1000 \f]
 *
 */

template< typename FieldT >
class SpeciesHeatCapacity_Cp
 : public Expr::Expression<FieldT>
{
  const Expr::Tag tTag_;
  const Expr::TagList tPowerTags_;
  const FieldT* t_;
  std::vector<const FieldT*> tPowers_;
  Cantera_CXX::IdealGasMix* const gasMix_;
  const int n_; //index of species to be evaluated

  SpeciesHeatCapacity_Cp( const Expr::Tag& tTag,
                          const int n );
public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a HeatCapacity expression
     *  @param resultTag tag for the pure species heat capacity
     *  @param tTag tag for temperature
     *  @param n species index
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::Tag& tTag,
             const int n);

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag tTag_;
    const int n_;
  };

  ~SpeciesHeatCapacity_Cp();
  void advertise_dependents( Expr::ExprDeps& exprDeps );
  void bind_fields( const Expr::FieldManagerList& fml );
  void evaluate();
};


// ###################################################################
//
//                          Implementation
//
// ###################################################################



template< typename FieldT >
HeatCapacity_Cp<FieldT>::
HeatCapacity_Cp( const Expr::Tag& tTag,
                 const Expr::Tag& massFracTag )
  : Expr::Expression<FieldT>(),
    tTag_( tTag ),
    tPowerTags_( TemperaturePowers<FieldT>::temperature_powers_tags() ),
    gasMix_( CanteraObjects::get_gasmix() ),
    nSpec_( gasMix_->nSpecies() )
{
  this->set_gpu_runnable( true );

  massFracTags_.clear();
  for( size_t n=0; n<nSpec_; ++n ){
    std::ostringstream name;
    name << massFracTag.name() << "_" << n;
    massFracTags_.push_back( Expr::Tag(name.str(),massFracTag.context()) );
  }

}

//--------------------------------------------------------------------

template< typename FieldT >
HeatCapacity_Cp<FieldT>::
~HeatCapacity_Cp()
{
  CanteraObjects::restore_gasmix(gasMix_);
}

//--------------------------------------------------------------------

template< typename FieldT >
void
HeatCapacity_Cp<FieldT>::
advertise_dependents( Expr::ExprDeps& exprDeps )
{
  exprDeps.requires_expression( tTag_ );
  exprDeps.requires_expression( tPowerTags_ );
  exprDeps.requires_expression( massFracTags_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
HeatCapacity_Cp<FieldT>::
bind_fields( const Expr::FieldManagerList& fml )
{
  const typename Expr::FieldMgrSelector<FieldT>::type& fm = fml.field_manager<FieldT>();

  t_ = &fm.field_ref( tTag_ );

  tPowers_.clear();
  BOOST_FOREACH( const Expr::Tag& tag, tPowerTags_ ){
    tPowers_.push_back( &fm.field_ref(tag) );
  }

  massFracs_.clear();
  BOOST_FOREACH( const Expr::Tag& tag, massFracTags_ ){
    massFracs_.push_back( &fm.field_ref(tag) );
  }
}

//--------------------------------------------------------------------

template< typename FieldT >
void
HeatCapacity_Cp<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  using namespace Cantera;
  FieldT& cp = this->value();

  const FieldT& t2 = *tPowers_[0]; // t^2
  const FieldT& t3 = *tPowers_[1]; // t^3
  const FieldT& t4 = *tPowers_[2]; // t^4
  const FieldT& recipRecipT = *tPowers_[5]; // t^-2

  std::vector<double> c(15,0); //vector of Cantera's coefficients
  int polyType;
  double minT;
  double maxT;
  double refPressure;
  const Cantera::SpeciesThermo& spThermo = gasMix_->speciesThermo();
  const std::vector<double>& molecularWeights = gasMix_->molecularWeights();

  cp <<= 0.0; // set cp to 0 before starting summation

  for( size_t n=0; n<nSpec_; ++n ){
    spThermo.reportParams(n, polyType, &c[0], minT, maxT, refPressure);
    switch (polyType) {
    case SIMPLE:
      cp <<= cp + *massFracs_[n] * c[3]
                                 / molecularWeights[n];
      break;
      /* polynomials are applicable in two temperature ranges - high and low
       * If the temperature is out of range, the value is set to the value at the min or max temp
       */
    case NASA2:
      for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic)
        *ic *= GasConstant / molecularWeights[n]; // dimensionalize coefficients
      cp <<= cp + *massFracs_[n] * cond( *t_ <= c[0] && *t_ >= minT, c[1] + c[2] * *t_ + c[ 3] * t2 + c[ 4] * t3 + c[ 5] * t4 ) // if low temp
                                       ( *t_ >  c[0] && *t_ <= maxT, c[8] + c[9] * *t_ + c[10] * t2 + c[11] * t3 + c[12] * t4 )  // else if high temp
                                       ( *t_ < minT, c[1] + c[2] * minT + c[ 3] * minT * minT + c[ 4] * pow(minT,3) + c[ 5] * pow(minT,4) )  // else if out of bounds - low
                                       (             c[8] + c[9] * maxT + c[10] * maxT * maxT + c[11] * pow(maxT,3) + c[12] * pow(maxT,4) ); // else out of bounds - high

      break;
    case SHOMATE2:
      double minTScaled = minT/1000;
      double maxTScaled = maxT/1000;
      for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic )
        *ic *= 1e3 / molecularWeights[n]; // scale the coefficients to keep units consistent
      cp <<= cp + *massFracs_[n] * cond( *t_ <= c[0] && *t_ >= minT, c[1] + c[2] * *t_*1e-3 + c[ 3] * t2*1e-6 + c[ 4] * t3*1e-9 + c[ 5] * recipRecipT*1e6 ) // if low temp
                                       ( *t_ >  c[0] && *t_ <= maxT, c[8] + c[9] * *t_*1e-3 + c[10] * t2*1e-6 + c[11] * t3*1e-9 + c[12] * recipRecipT*1e6 )  // else if high temp
                                       ( *t_ < minT, c[1] + c[2] * minTScaled + c[ 3] * minTScaled * minTScaled + c[ 4] * pow(minTScaled,3) + c[ 5] * pow(minTScaled,-2) )  // else if out of bounds - low
                                       (             c[8] + c[9] * maxTScaled + c[10] * maxTScaled * maxTScaled + c[11] * pow(maxTScaled,3) + c[12] * pow(maxTScaled,-2) ); // else out of bounds - high

      break;
    }
  }
}

//--------------------------------------------------------------------

template< typename FieldT >
HeatCapacity_Cp<FieldT>::
Builder::Builder( const Expr::Tag& resultTag,
                  const Expr::Tag& tTag,
                  const Expr::Tag& massFracTag )
: ExpressionBuilder( resultTag ),
  tTag_( tTag ),
  massFracTag_( massFracTag )
{}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
HeatCapacity_Cp<FieldT>::
Builder::build() const
{
  return new HeatCapacity_Cp<FieldT>( tTag_, massFracTag_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
SpeciesHeatCapacity_Cp<FieldT>::
SpeciesHeatCapacity_Cp( const Expr::Tag& tTag,
                        const int n )
  : Expr::Expression<FieldT>(),
    tTag_( tTag ),
    tPowerTags_( TemperaturePowers<FieldT>::temperature_powers_tags() ),
    gasMix_( CanteraObjects::get_gasmix() ),
    n_ ( n )
{
  this->set_gpu_runnable( true );
}

//--------------------------------------------------------------------

template< typename FieldT >
SpeciesHeatCapacity_Cp<FieldT>::
~SpeciesHeatCapacity_Cp()
{
  CanteraObjects::restore_gasmix(gasMix_);
}

//--------------------------------------------------------------------

template< typename FieldT >
void
SpeciesHeatCapacity_Cp<FieldT>::
advertise_dependents( Expr::ExprDeps& exprDeps )
{
  exprDeps.requires_expression( tTag_ );
  exprDeps.requires_expression( tPowerTags_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
SpeciesHeatCapacity_Cp<FieldT>::
bind_fields( const Expr::FieldManagerList& fml )
{
  const typename Expr::FieldMgrSelector<FieldT>::type& fm = fml.field_manager<FieldT>();

  t_ = &fm.field_ref( tTag_ );

  tPowers_.clear();
  BOOST_FOREACH( const Expr::Tag& tag, tPowerTags_ ){
    tPowers_.push_back( &fm.field_ref(tag) );
  }
}

//--------------------------------------------------------------------

template< typename FieldT >
void
SpeciesHeatCapacity_Cp<FieldT>::
evaluate()
{
  using namespace Cantera;
  using namespace SpatialOps;
  FieldT& cp = this->value();

  const FieldT& t2 = *tPowers_[0]; // t^2
  const FieldT& t3 = *tPowers_[1]; // t^3
  const FieldT& t4 = *tPowers_[2]; // t^4
  const FieldT& recipRecipT = *tPowers_[5]; // t^-2

  std::vector<double> c(15,0); // vector of Cantera's coefficients
  int polyType; // type of polynomial
  double minT;  // minimum temperature where polynomial is valid
  double maxT;  // maximum temperature where polynomial is valid
  double refPressure;
  gasMix_->speciesThermo().reportParams( n_, polyType, &c[0], minT, maxT, refPressure );
  const std::vector<double>& molecularWeights = gasMix_->molecularWeights();

  switch (polyType) {
    case SIMPLE:
      cp <<= c[3] / molecularWeights[n_];
      break;
    case NASA2:
      for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic){
        *ic *= GasConstant / molecularWeights[n_]; // dimensionalize coefficients
      }
      cp <<= cond( *t_ <= c[0] && *t_ >= minT, c[1] + c[2] * *t_ + c[ 3] * t2 + c[ 4] * t3 + c[ 5] * t4 ) // if low temp
                 ( *t_ >  c[0] && *t_ <= maxT, c[8] + c[9] * *t_ + c[10] * t2 + c[11] * t3 + c[12] * t4 )  // else if high temp
                 ( *t_ < minT, c[1] + c[2] * minT + c[ 3] * minT * minT + c[ 4] * pow(minT,3) + c[ 5] * pow(minT,4) )  // else if out of bounds - low
                 (             c[8] + c[9] * maxT + c[10] * maxT * maxT + c[11] * pow(maxT,3) + c[12] * pow(maxT,4) ); // else out of bounds - high

      break;
    case SHOMATE2:
      double minTScaled = minT/1000;
      double maxTScaled = maxT/1000;
      for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic ){
        *ic *= 1e3 / molecularWeights[n_]; // scale the coefficients to keep units consistent
      }
      cp <<= cond( *t_ <= c[0] && *t_ >= minT, c[1] + c[2] * *t_*1e-3 + c[ 3] * t2*1e-6 + c[ 4] * t3*1e-9 + c[ 5] * recipRecipT*1e6 ) // if low temp
                 ( *t_ >  c[0] && *t_ <= maxT, c[8] + c[9] * *t_*1e-3 + c[10] * t2*1e-6 + c[11] * t3*1e-9 + c[12] * recipRecipT*1e6 )  // else if high temp
                 ( *t_ < minT, c[1] + c[2] * minTScaled + c[ 3] * minTScaled * minTScaled + c[ 4] * pow(minTScaled,3) + c[ 5] * pow(minTScaled,-2) )  // else if out of bounds - low
                 (             c[8] + c[9] * maxTScaled + c[10] * maxTScaled * maxTScaled + c[11] * pow(maxTScaled,3) + c[12] * pow(maxTScaled,-2) ); // else out of bounds - high

      break;
  }
}

//--------------------------------------------------------------------

template< typename FieldT >
SpeciesHeatCapacity_Cp<FieldT>::
Builder::Builder( const Expr::Tag& resultTag,
                  const Expr::Tag& tTag,
                  const int n )
  : ExpressionBuilder( resultTag ),
    tTag_( tTag ),
    n_( n )
{}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
SpeciesHeatCapacity_Cp<FieldT>::
Builder::build() const
{
  return new SpeciesHeatCapacity_Cp<FieldT>( tTag_, n_ );
}

#endif // HeatCapacity_Cp_Expr_h
