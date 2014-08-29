#ifndef HeatCapacity_Cv_Expr_h
#define HeatCapacity_Cv_Expr_h

#include <expression/Expression.h>

#include <pokitt/CanteraObjects.h> //include cantera wrapper
#include <pokitt/thermo/Temperature.h>

#include <cantera/kernel/ct_defs.h> // contains value of Cantera::GasConstant
#include <cantera/kernel/speciesThermoTypes.h> // contains definitions for which polynomial is being used
#include <cantera/IdealGasMix.h>

#include <exception>

namespace Cantera_CXX{ class IdealGasMix; } // location of polynomial coefficients

/**
 *  \class  HeatCapacity_Cv
 *  \author Nate Yonkee
 *  \date May, 2014
 *
 *  \brief Calculates the constant volume heat capacity of each species
 *  using either NASA7 or Shomate polynomials with 2 temperature ranges.
 *  Units of J/kg/K
 *
 *  Five coefficients \f$(a_0,\dots,a_4)\f$ are used to represent
 * \f$ c_v^0(T)\f$ as a polynomial in \f$ T \f$
 *
 * NASA7:
 *
 * \f[
 * \frac{c_v(T)}{R} = -1 + a_0 + a_1 T + a_2 T^2 + a_3 T^3 + a_4 T^4
 * \f]
 *
 * Shomate:
 *
 * \f[
 * \frac{c_v(T)}{1000} = \frac{-R}{1000} + a_0 + a_1 t + a_2 t^2 + a_3 t^3 + \frac{a_4}{t^2}
 * \f]
 *
 * Where \f[ t = T (K)/1000 \f]
 *
 * The heat capacity is a weighted average of the pure species \f$ c_v^0(T)\f$,
 *
 * \f[
 * c_v(T) = \sum_{n=0}^{nSpec} Y_n c_{v,n}
 * \f]
 */

template< typename FieldT >
class HeatCapacity_Cv
 : public Expr::Expression<FieldT>
{
  typedef std::vector<double> PolyVals; // values used for polynomial
  const Expr::Tag tTag_;
  const Expr::TagList tPowerTags_;
  Expr::TagList massFracTags_;
  const FieldT* t_;
  std::vector<const FieldT*> tPowers_;
  std::vector<const FieldT*> massFracs_;

  int nSpec_; // number of species
  PolyVals minTVec_; // vector of minimum temperatures for polynomial evaluations
  PolyVals maxTVec_; // vector of maximum temperatures for polynomial evaluations
  std::vector< PolyVals > cVec_; // vector of polynomial coefficients
  std::vector<int> polyTypeVec_; // vector of polynomial types

  HeatCapacity_Cv( const Expr::Tag& tTag,
                   const Expr::Tag& massFracTag );
public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a HeatCapacity_Cv expression
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

  ~HeatCapacity_Cv();
  void advertise_dependents( Expr::ExprDeps& exprDeps );
  void bind_fields( const Expr::FieldManagerList& fml );
  void evaluate();

};

/**
 *  \class  SpeciesHeatCapacity_Cv
 *  \author Nate Yonkee
 *  \date May, 2014
 *
 *  \brief Calculates the constant volume heat capacity of a single species
 *  using either NASA7 or Shomate polynomials with 2 temperature ranges.
 *  Units of J/kg/K
 *
 *  Five coefficients \f$(a_0,\dots,a_4)\f$ are used to represent
 * \f$ c_v^0(T)\f$ as a polynomial in \f$ T \f$
 *
 * NASA7:
 *
 * \f[
 * \frac{c_v(T)}{R} = a_0 + a_1 T + a_2 T^2 + a_3 T^3 + a_4 T^4
 * \f]
 *
 * Shomate:
 *
 * \f[
 * \frac{c_v(T)}{1000} = a_0 + a_1 t + a_2 t^2 + a_3 t^3 + \frac{a_4}{t^2}
 * \f]
 *
 * Where \f[ t = T (K)/1000 \f]
 *
 */

template< typename FieldT >
class SpeciesHeatCapacity_Cv
 : public Expr::Expression<FieldT>
{
  const Expr::Tag tTag_;
  const Expr::TagList tPowerTags_;
  const FieldT* t_;
  std::vector<const FieldT*> tPowers_;

  const int n_; //index of species to be evaluated
  double minT_; // minimum temperature for polynomial evaluation
  double maxT_; // maximum temperature for polynomial evaluation
  std::vector<double> c_; // vector of polynomial coefficients
  int polyType_; // polynomial type

  SpeciesHeatCapacity_Cv( const Expr::Tag& tTag,
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

  ~SpeciesHeatCapacity_Cv();
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
HeatCapacity_Cv<FieldT>::
HeatCapacity_Cv( const Expr::Tag& tTag,
                 const Expr::Tag& massFracTag )
  : Expr::Expression<FieldT>(),
    tTag_( tTag ),
    tPowerTags_( Temperature<FieldT>::temperature_powers_tags() )
{
  this->set_gpu_runnable( true );

  Cantera_CXX::IdealGasMix* const gasMix = CanteraObjects::get_gasmix();
  const Cantera::SpeciesThermo& spThermo = gasMix->speciesThermo();

  nSpec_ = gasMix->nSpecies();
  massFracTags_.clear();
  for( size_t n=0; n<nSpec_; ++n ){
    std::ostringstream name;
    name << massFracTag.name() << "_" << n;
    massFracTags_.push_back( Expr::Tag(name.str(),massFracTag.context()) );
  }

  const std::vector<double> molecularWeights = gasMix->molecularWeights();
  std::vector<double> c(15,0); //vector of Cantera's coefficients
  int polyType;
  double minT;
  double maxT;
  double refPressure;

  for( size_t n=0; n<nSpec_; ++n ){
    spThermo.reportParams(n, polyType, &c[0], minT, maxT, refPressure);
    polyTypeVec_.push_back(polyType); // vector of polynomial types
    minTVec_.push_back(minT); // vector of minimum temperatures for polynomial evaluations
    maxTVec_.push_back(maxT); // vector of maximum temperatures for polynomial evaluations
    switch (polyType) {
    case SIMPLE:
      c[3] -= Cantera::GasConstant; // change coefficients from isobaric to isometric
      c[3] /= molecularWeights[n]; // convert to mass basis
      break;
    case NASA2:
      c[1] -= 1.0; //change coefficients from isobaric to isometric
      c[8] -= 1.0;
      for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic)
        *ic *= Cantera::GasConstant / molecularWeights[n]; // dimensionalize the coefficients to mass basis
      break;
    case SHOMATE2:
      c[1] -= Cantera::GasConstant * 1e-3; //change coefficients from isobaric to isometric
      c[8] -= Cantera::GasConstant * 1e-3;
      for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic ){
        *ic *= 1e3 / molecularWeights[n]; // scale the coefficients to keep units consistent on mass basis
      }
      break;
    }
    cVec_.push_back(c); // vector of polynomial coefficients
  }

  CanteraObjects::restore_gasmix(gasMix);
}

//--------------------------------------------------------------------

template< typename FieldT >
HeatCapacity_Cv<FieldT>::
~HeatCapacity_Cv()
{

}

//--------------------------------------------------------------------

template< typename FieldT >
void
HeatCapacity_Cv<FieldT>::
advertise_dependents( Expr::ExprDeps& exprDeps )
{
  exprDeps.requires_expression( tTag_ );
  exprDeps.requires_expression( tPowerTags_ );
  exprDeps.requires_expression( massFracTags_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
HeatCapacity_Cv<FieldT>::
bind_fields( const Expr::FieldManagerList& fml )
{
  const typename Expr::FieldMgrSelector<FieldT>::type& fm = fml.field_manager<FieldT>();

  t_ = &fml.template field_ref< FieldT >( tTag_ );

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
HeatCapacity_Cv<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  FieldT& cv = this->value();

  const FieldT& t2 =          *tPowers_[0]; // t^2
  const FieldT& t3 =          *tPowers_[1]; // t^3
  const FieldT& t4 =          *tPowers_[2]; // t^4
  const FieldT& recipRecipT = *tPowers_[5]; // t^-2

  cv <<= 0.0; // set cv to 0 before starting summation

  for( size_t n=0; n<nSpec_; ++n ){
    int polyType = polyTypeVec_[n];
    std::vector<double>& c = cVec_[n];
    double minT = minTVec_[n];
    double maxT = maxTVec_[n];
    switch (polyType) {
    case SIMPLE:
      cv <<= cv + *massFracs_[n] * c[3];
      break;
    case NASA2:
      cv <<= cv + *massFracs_[n] * cond( *t_ <= c[0] && *t_ >= minT, c[1] + c[2] * *t_ + c[ 3] * t2 + c[ 4] * t3 + c[ 5] * t4) // if low temp
                                       ( *t_ >  c[0] && *t_ <= maxT, c[8] + c[9] * *t_ + c[10] * t2 + c[11] * t3 + c[12] * t4)  // else if high temp
                                       ( *t_ < minT, c[1] + c[2] * minT + c[ 3] * minT * minT + c[ 4] * pow(minT,3) + c[ 5] * pow(minT,4))  // else if out of bounds - low
                                       (             c[8] + c[9] * maxT + c[10] * maxT * maxT + c[11] * pow(maxT,3) + c[12] * pow(maxT,4)); // else out of bounds - high
    break;
  case SHOMATE2:
      cv <<= cv + *massFracs_[n] * cond( *t_ <= c[0] && *t_ >= minT, c[1] + c[2] * *t_*1e-3 + c[ 3] * t2*1e-6 + c[ 4] * t3*1e-9 + c[ 5] * recipRecipT*1e6) // if low temp
                                       ( *t_ >  c[0] && *t_ <= maxT, c[8] + c[9] * *t_*1e-3 + c[10] * t2*1e-6 + c[11] * t3*1e-9 + c[12] * recipRecipT*1e6)  // else if high temp
                                       ( *t_ < minT, c[1] + c[2] * minT*1e-3 + c[ 3] * minT*1e-3 * minT*1e-3 + c[ 4] * pow(minT*1e-3,3) + c[ 5] * pow(minT*1e-3,-2))  // else if out of bounds - low
                                       (             c[8] + c[9] * maxT*1e-3 + c[10] * maxT*1e-3 * maxT*1e-3 + c[11] * pow(maxT*1e-3,3) + c[12] * pow(maxT*1e-3,-2)); // else out of bounds - high
    }
  }
}

//--------------------------------------------------------------------

template< typename FieldT >
HeatCapacity_Cv<FieldT>::
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
HeatCapacity_Cv<FieldT>::
Builder::build() const
{
  return new HeatCapacity_Cv<FieldT>( tTag_, massFracTag_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
SpeciesHeatCapacity_Cv<FieldT>::
SpeciesHeatCapacity_Cv( const Expr::Tag& tTag,
                        const int n )
  : Expr::Expression<FieldT>(),
    tTag_( tTag ),
    tPowerTags_( Temperature<FieldT>::temperature_powers_tags() ),
    n_ ( n )
{
  this->set_gpu_runnable( true );

  Cantera_CXX::IdealGasMix* const gasMix = CanteraObjects::get_gasmix();
  const Cantera::SpeciesThermo& spThermo = gasMix->speciesThermo();

  double molecularWeight = gasMix->molecularWeight(n_);
  c_.resize(15); //vector of Cantera's coefficients
  double refPressure;
  spThermo.reportParams(n_, polyType_, &c_[0], minT_, maxT_, refPressure);
  switch (polyType_) {
  case SIMPLE:
    c_[3] -= Cantera::GasConstant; // change coefficients from isobaric to isometric
    c_[3] /= molecularWeight; // convert to mass basis
    break;
  case NASA2:
    c_[1] -= 1.0; // change coefficients from isobaric to isometric
    c_[8] -= 1.0;
    for( std::vector<double>::iterator ic = c_.begin() + 1; ic!=c_.end(); ++ic){
      *ic *= Cantera::GasConstant / molecularWeight; // dimensionalize the coefficients to mass basis
    }
    break;
  case SHOMATE2:
    c_[1] -= Cantera::GasConstant * 1e-3; // change coefficients from isobaric to isometric
    c_[8] -= Cantera::GasConstant * 1e-3;
    for( std::vector<double>::iterator ic = c_.begin() + 1; ic!=c_.end(); ++ic ){
      *ic *= 1e3 / molecularWeight; // scale the coefficients to keep units consistent on mass basis
    }
    break;
  }

  CanteraObjects::restore_gasmix(gasMix);
}

//--------------------------------------------------------------------

template< typename FieldT >
SpeciesHeatCapacity_Cv<FieldT>::
~SpeciesHeatCapacity_Cv()
{

}

//--------------------------------------------------------------------

template< typename FieldT >
void
SpeciesHeatCapacity_Cv<FieldT>::
advertise_dependents( Expr::ExprDeps& exprDeps )
{
  exprDeps.requires_expression( tTag_ );
  exprDeps.requires_expression( tPowerTags_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
SpeciesHeatCapacity_Cv<FieldT>::
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
SpeciesHeatCapacity_Cv<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  FieldT& cv = this->value();

  const FieldT& t2 =          *tPowers_[0]; // t^2
  const FieldT& t3 =          *tPowers_[1]; // t^3
  const FieldT& t4 =          *tPowers_[2]; // t^4
  const FieldT& recipRecipT = *tPowers_[5]; // t^-2

  switch (polyType_) {
    case SIMPLE:
      cv <<= c_[3];
      break;
    case NASA2:
      cv <<= cond( *t_ <= c_[0] && *t_ >= minT_, c_[1] + c_[2] * *t_ + c_[ 3] * t2 + c_[ 4] * t3 + c_[ 5] * t4 ) // if low temp
                 ( *t_ >  c_[0] && *t_ <= maxT_, c_[8] + c_[9] * *t_ + c_[10] * t2 + c_[11] * t3 + c_[12] * t4 )  // else if high temp
                 ( *t_ < minT_, c_[1] + c_[2] * minT_ + c_[ 3] * minT_ * minT_ + c_[ 4] * pow(minT_,3) + c_[ 5] * pow(minT_,4) )  // else if out of bounds - low
                 (              c_[8] + c_[9] * maxT_ + c_[10] * maxT_ * maxT_ + c_[11] * pow(maxT_,3) + c_[12] * pow(maxT_,4) ); // else out of bounds - high

      break;
    case SHOMATE2:
      cv <<= cond( *t_ <= c_[0] && *t_ >= minT_, c_[1] + c_[2] * *t_*1e-3 + c_[ 3] * t2*1e-6 + c_[ 4] * t3*1e-9 + c_[ 5] * recipRecipT*1e6 ) // if low temp
                 ( *t_ >  c_[0] && *t_ <= maxT_, c_[8] + c_[9] * *t_*1e-3 + c_[10] * t2*1e-6 + c_[11] * t3*1e-9 + c_[12] * recipRecipT*1e6 )  // else if high temp
                 ( *t_ < minT_, c_[1] + c_[2] * minT_*1e-3 + c_[ 3] * minT_*1e-3 * minT_*1e-3 + c_[ 4] * pow(minT_*1e-3,3) + c_[ 5] * pow(minT_*1e-3,-2) )  // else if out of bounds - low
                 (              c_[8] + c_[9] * maxT_*1e-3 + c_[10] * maxT_*1e-3 * maxT_*1e-3 + c_[11] * pow(maxT_*1e-3,3) + c_[12] * pow(maxT_*1e-3,-2) ); // else out of bounds - high

      break;
  }
}

//--------------------------------------------------------------------

template< typename FieldT >
SpeciesHeatCapacity_Cv<FieldT>::
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
SpeciesHeatCapacity_Cv<FieldT>::
Builder::build() const
{
  return new SpeciesHeatCapacity_Cv<FieldT>( tTag_, n_ );
}

#endif // HeatCapacity_Cv_Expr_h
