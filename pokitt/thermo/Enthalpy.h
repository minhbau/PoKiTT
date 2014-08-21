#ifndef Enthalpy_Expr_h
#define Enthalpy_Expr_h

#include <expression/Expression.h>

#include <cantera/kernel/ct_defs.h> // contains value of gas constant
#include <cantera/kernel/SpeciesThermoInterpType.h> // contains definitions for which polynomial is being used
#include <cantera/IdealGasMix.h>

#include <pokitt/CanteraObjects.h> //include cantera wrapper

namespace Cantera_CXX{ class IdealGasMix; } //location of polynomial

/**
 *  \class  Enthalpy
 *  \author Nate Yonkee
 *  \date July, 2014
 *
 *  \brief Calculates the enthalpy of a mixture
 *  using either NASA7 or Shomate polynomials with 2 temperature ranges.
 *  Units of J/kg
 *
 *  Five coefficients \f$(a_0,\dots,a_4)\f$ are used to represent
 * \f$ h^0(T)\f$ as a polynomial in \f$ T \f$ e.g. NASA7:
 *
 * \f[
 * \frac{h^0(T)}{R} = a_0 T + a_1 T^2/2 + a_2 T^3/3 + a_3 T^4/4 + a_4 T^5/5 + a_5
 * \f]
 *
 * The enthalpy is a weighted average of the pure species \f$ h^0(T)\f$,
 *
 * \f[
 * h(T) = \sum_{n=0}^{nSpec} x_n h^0_n
 * \f]
 */

template< typename FieldT >
class Enthalpy
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

  Enthalpy( const Expr::Tag& tTag,
            const Expr::Tag& massFracTag );
public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a Enthalpy expression
     *  @param resultTag tag for the mixture enthalpy
     *  @param tTag tag for temperature
     *  @param massFracTag massFracTag tag for mass fractions of each species, ordering is consistent with Cantera input
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::Tag& tTag,
             const Expr::Tag& massFracTag );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag tTag_;
    const Expr::Tag massFracTag_;
  };

  ~Enthalpy();
  void advertise_dependents( Expr::ExprDeps& exprDeps );
  void bind_fields( const Expr::FieldManagerList& fml );
  void evaluate();

};

/**
 *  \class  SpeciesEnthalpy
 *  \author Nate Yonkee
 *  \date July, 2014
 *
 *  \brief Calculates the enthalpy of each species
 *  using either NASA7 or Shomate polynomials with 2 temperature ranges.
 *  Units of J/kg
 *
 *  Five coefficients \f$(a_0,\dots,a_4)\f$ are used to represent
 * \f$ h(T)\f$ as a polynomial in \f$ T \f$ e.g. NASA7:
 *
 * \f[
 * \frac{h(T)}{R} = a_0 T + a_1 T^2/2 + a_2 T^3/3 + a_3 T^4/4 + a_4 T^5/5 + a_5
 * \f]
 *
 */
template< typename FieldT >
class SpeciesEnthalpy
 : public Expr::Expression<FieldT>
{
  const Expr::Tag tTag_;
  const Expr::TagList tPowerTags_;
  const FieldT* t_;
  std::vector<const FieldT*> tPowers_;
  Cantera_CXX::IdealGasMix* const gasMix_;
  const int n_; //index of species to be evaluated

  SpeciesEnthalpy( const Expr::Tag& tTag,
                   const int n );
public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a SpeciesEnthalpy expression
     *  @param resultTag tag for the pure species enthalpy
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

  ~SpeciesEnthalpy();
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
Enthalpy<FieldT>::
Enthalpy( const Expr::Tag& tTag,
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
Enthalpy<FieldT>::
~Enthalpy()
{
  CanteraObjects::restore_gasmix(gasMix_);
}

//--------------------------------------------------------------------

template< typename FieldT >
void
Enthalpy<FieldT>::
advertise_dependents( Expr::ExprDeps& exprDeps )
{
  exprDeps.requires_expression( tTag_ );
  exprDeps.requires_expression( tPowerTags_ );
  exprDeps.requires_expression( massFracTags_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
Enthalpy<FieldT>::
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
Enthalpy<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  using namespace Cantera;
  FieldT& h = this->value();

  const FieldT& t2 = *tPowers_[0]; // t^2
  const FieldT& t3 = *tPowers_[1]; // t^3
  const FieldT& t4 = *tPowers_[2]; // t^4
  const FieldT& t5 = *tPowers_[3]; // t^5
  const FieldT& recipT = *tPowers_[4]; // t^-1

  std::vector<double> c(15,0); //vector of Cantera's coefficients
  int polyType;
  double minT;
  double maxT;
  double refPressure;
  const Cantera::SpeciesThermo& spThermo = gasMix_->speciesThermo();
  const std::vector<double>& molecularWeights = gasMix_->molecularWeights();

  h <<= 0.0;

  for( size_t n=0; n<nSpec_; ++n ){
    spThermo.reportParams(n, polyType, &c[0], minT, maxT, refPressure);
    switch (polyType) {
    case SIMPLE: // constant heat capacity
      h <<= h + *massFracs_[n] * ( c[1] + c[3] * (*t_-c[0]) )
                                / molecularWeights[n];
      break;
    case NASA2:
      for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic)
        *ic *= GasConstant / molecularWeights[n]; // dimensionalize the coefficients
      h <<= h + *massFracs_[n] * cond( *t_ <= c[0] && *t_ >= minT, c[ 6] + c[1] * *t_ + c[2]/2 * t2 + c[ 3]/3 * t3 + c[ 4]/4 * t4 + c[ 5]/5 * t5) // if low temp
                                     ( *t_ >  c[0] && *t_ <= maxT, c[13] + c[8] * *t_ + c[9]/2 * t2 + c[10]/3 * t3 + c[11]/4 * t4 + c[12]/5 * t5)  // else if high temp
                                     ( *t_ < minT, c[1] * *t_ + c[2] * minT * (*t_-minT/2) + c[ 3] * minT * minT * (*t_-2*minT/3) + c[ 4] * pow(minT,3) * (*t_-3*minT/4) + c[ 5] * pow(minT,4) * (*t_-4*minT/5) + c[ 6]) // else if out of bounds - low
                                     (              c[8] * *t_ + c[9] * maxT * (*t_-maxT/2) + c[10] * maxT * maxT * (*t_-2*maxT/3) + c[11] * pow(maxT,3) * (*t_-3*maxT/4) + c[12] * pow(maxT,4) * (*t_-4*maxT/5) + c[13]); // else out of bounds - high
      break;
    case SHOMATE2:
      double minTScaled = minT/1000;
      double maxTScaled = maxT/1000;
      for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic )
        *ic *= 1e6 / molecularWeights[n]; // scale the coefficients to keep units consistent
      h <<= h + *massFracs_[n] * cond( *t_ <= c[0] && *t_ >= minT, c[ 6] + c[1]**t_*1e-3 + c[2]/2 * t2*1e-6 + c[ 3]/3 * t3*1e-9 + c[ 4]/4 * t4*1e-12 - c[ 5] * recipT*1e3 ) // if low temp
                                     ( *t_ >  c[0] && *t_ <= maxT, c[13] + c[8]**t_*1e-3 + c[9]/2 * t2*1e-6 + c[10]/3 * t3*1e-9 + c[11]/4 * t4*1e-12 - c[12] * recipT*1e3 )  // else if high temp
                                     ( *t_ < minT, c[1] * *t_*1e-3 + c[2] * minTScaled * ( *t_*1e-3 - minTScaled/2 ) + c[ 3] * minTScaled * minTScaled * ( *t_*1e-3 - 2*minTScaled/3 ) + c[ 4] * pow(minTScaled,3) * ( *t_*1e-3 - 3*minTScaled/4 ) - c[ 5] * pow(minTScaled,-1) * (-*t_*1e-3 / minTScaled + 2 ) + c[ 6] ) // else if out of bounds - low
                                     (              c[8] * *t_*1e-3 + c[9] * maxTScaled * ( *t_*1e-3 - maxTScaled/2 ) + c[10] * maxTScaled * maxTScaled * ( *t_*1e-3 - 2*maxTScaled/3 ) + c[11] * pow(maxTScaled,3) * ( *t_*1e-3 - 3*maxTScaled/4 ) - c[12] * pow(maxTScaled,-1) * (-*t_*1e-3 / maxTScaled + 2 ) + c[13] ); // else out of bounds - high
  }
}
}
//--------------------------------------------------------------------

template< typename FieldT >
Enthalpy<FieldT>::
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
Enthalpy<FieldT>::
Builder::build() const
{
  return new Enthalpy<FieldT>( tTag_, massFracTag_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
SpeciesEnthalpy<FieldT>::
SpeciesEnthalpy( const Expr::Tag& tTag,
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
SpeciesEnthalpy<FieldT>::
~SpeciesEnthalpy()
{
  CanteraObjects::restore_gasmix(gasMix_);
}

//--------------------------------------------------------------------

template< typename FieldT >
void
SpeciesEnthalpy<FieldT>::
advertise_dependents( Expr::ExprDeps& exprDeps )
{
  exprDeps.requires_expression( tTag_ );
  exprDeps.requires_expression( tPowerTags_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
SpeciesEnthalpy<FieldT>::
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
SpeciesEnthalpy<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  using namespace Cantera;
  FieldT& h = this->value();

  const FieldT& t2 = *tPowers_[0]; // t^2
    const FieldT& t3 = *tPowers_[1]; // t^3
    const FieldT& t4 = *tPowers_[2]; // t^4
    const FieldT& t5 = *tPowers_[3]; // t^5
    const FieldT& recipT = *tPowers_[4]; // t^-1

    std::vector<double> c(15,0); // vector of Cantera's coefficients
    int polyType; // type of polynomial
    double minT;  // minimum temperature where polynomial is valid
    double maxT;  // maximum temperature where polynomial is valid
    double refPressure;
    gasMix_->speciesThermo().reportParams( n_, polyType, &c[0], minT, maxT, refPressure );
    const std::vector<double>& molecularWeights = gasMix_->molecularWeights();

    switch (polyType) {
      case SIMPLE:
      h <<= ( c[1] + c[3] * (*t_-c[0]) )
                         / molecularWeights[n_];
      break;
    case NASA2:
      for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic){
        *ic *= GasConstant / molecularWeights[n_]; // dimensionalize coefficients
      }
      h <<= cond( *t_ <= c[0] && *t_ >= minT, c[ 6] + c[1] * *t_ + c[2]/2 * t2 + c[ 3]/3 * t3 + c[ 4]/4 * t4 + c[ 5]/5 * t5) // if low temp
                             ( *t_ >  c[0] && *t_ <= maxT, c[13] + c[8] * *t_ + c[9]/2 * t2 + c[10]/3 * t3 + c[11]/4 * t4 + c[12]/5 * t5)  // else if high temp
                             ( *t_ <  minT                 , c[1] * *t_ + c[2] * minT * (*t_-minT/2) + c[ 3] * minT * minT * (*t_ - 2*minT/3) + c[ 4] * pow(minT,3) * (*t_ - 3*minT/4) + c[ 5] * pow(minT,4) * (*t_ - 4*minT/5) + c[ 6]) // else if out of bounds - low
                             (                                c[8] * *t_ + c[9] * maxT * (*t_-maxT/2) + c[10] * maxT * maxT * (*t_ - 2*maxT/3) + c[11] * pow(maxT,3) * (*t_ - 3*maxT/4) + c[12] * pow(maxT,4) * (*t_ - 4*maxT/5) + c[13]); // else out of bounds - high
    break;
   case SHOMATE2:
     for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic ){
       *ic *= 1e6 / molecularWeights[n_]; // scale the coefficients to keep units consistent
     }
      double minTScaled = minT/1000;
      double maxTScaled = maxT/1000;
      h <<= cond( *t_ <= c[0] && *t_ >= minT, c[ 6] + c[1] * *t_*1e-3 + c[2]/2 * t2*1e-6 + c[ 3]/3 * t3*1e-9 + c[ 4]/4 * t4*1e-12 - c[ 5] * recipT*1e3 ) // if low temp
                             ( *t_ >  c[0] && *t_ <= maxT, c[13] + c[8] * *t_*1e-3 + c[9]/2 * t2*1e-6 + c[10]/3 * t3*1e-9 + c[11]/4 * t4*1e-12 - c[12] * recipT*1e3 )  // else if high range
                             ( *t_ <  minT, c[1] * *t_*1e-3 + c[2] * minTScaled * ( *t_*1e-3 - minTScaled/2 ) + c[ 3] * minTScaled * minTScaled * ( *t_*1e-3 - 2*minTScaled/3 ) + c[ 4] * pow(minTScaled,3) * ( *t_*1e-3 - 3*minTScaled/4 ) - c[ 5] * pow(minTScaled,-1) * ( -*t_*1e-3/minTScaled + 2 ) + c[ 6] ) // else if out of bounds - low
                             (               c[8] * *t_*1e-3 + c[9] * maxTScaled * ( *t_*1e-3 - maxTScaled/2 ) + c[10] * maxTScaled * maxTScaled * ( *t_*1e-3 - 2*maxTScaled/3 ) + c[11] * pow(maxTScaled,3) * ( *t_*1e-3 - 3*maxTScaled/4 ) - c[12] * pow(maxTScaled,-1) * ( -*t_*1e-3/maxTScaled + 2 ) + c[13] ); // else out of bounds - high
    break;
    }
  }

//--------------------------------------------------------------------

template< typename FieldT >
SpeciesEnthalpy<FieldT>::
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
SpeciesEnthalpy<FieldT>::
Builder::build() const
{
  return new SpeciesEnthalpy<FieldT>( tTag_, n_ );
}



#endif // Enthalpy_Expr_h
