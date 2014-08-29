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
 * h(T) = \sum_{n=0}^{nSpec} y_n h^0_n
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

  int nSpec_; // number of species
  std::vector<double> minTVec_; // vector of minimum temperatures for polynomial evaluations
  std::vector<double> maxTVec_; // vector of maximum temperatures for polynomial evaluations
  std::vector< std::vector<double> > cVec_; // vector of polynomial coefficients
  std::vector<int> polyTypeVec_; // vector of polynomial types

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

  const int n_; //index of species to be evaluated
  double minT_; // minimum temperature for polynomial evaluation
  double maxT_; // maximum temperature for polynomial evaluation
  std::vector<double> c_; // vector of polynomial coefficients
  int polyType_; // polynomial type


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
    tPowerTags_( TemperaturePowers<FieldT>::temperature_powers_tags() )
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
      c[1] /= molecularWeights[n]; // convert to mass basis
      c[3] /= molecularWeights[n]; // convert to mass basis
      break;
    case NASA2:
      for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic)
        *ic *= Cantera::GasConstant / molecularWeights[n]; // dimensionalize the coefficients to mass basis
      break;
    case SHOMATE2:
      for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic ){
        *ic *= 1e6 / molecularWeights[n]; // scale the coefficients to keep units consistent on mass basis
      }
      break;
    }
    cVec_.push_back(c); // vector of polynomial coefficients
  }

  CanteraObjects::restore_gasmix(gasMix);
}

//--------------------------------------------------------------------

template< typename FieldT >
Enthalpy<FieldT>::
~Enthalpy()
{

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

  const FieldT& t2 =     *tPowers_[0]; // t^2
  const FieldT& t3 =     *tPowers_[1]; // t^3
  const FieldT& t4 =     *tPowers_[2]; // t^4
  const FieldT& t5 =     *tPowers_[3]; // t^5
  const FieldT& recipT = *tPowers_[4]; // t^-1

  h <<= 0.0;

  for( size_t n=0; n<nSpec_; ++n ){
    int polyType = polyTypeVec_[n];
    std::vector<double>& c = cVec_[n];
    double minT = minTVec_[n];
    double maxT = maxTVec_[n];
    switch (polyType) {
      case SIMPLE: // constant heat capacity
        h <<= h + *massFracs_[n] * ( c[1] + c[3] * (*t_ - c[0]) );
        break;
      case NASA2:
        h <<= h + *massFracs_[n] * cond( *t_ <= c[0] && *t_ >= minT, c[ 6] + c[1] * *t_ + c[2]/2 * t2 + c[ 3]/3 * t3 + c[ 4]/4 * t4 + c[ 5]/5 * t5) // if low temp
                                       ( *t_ >  c[0] && *t_ <= maxT, c[13] + c[8] * *t_ + c[9]/2 * t2 + c[10]/3 * t3 + c[11]/4 * t4 + c[12]/5 * t5)  // else if high temp
                                       ( *t_ < minT, c[1] * *t_ + c[2] * minT * (*t_-minT/2) + c[ 3] * minT * minT * (*t_-2*minT/3) + c[ 4] * pow(minT,3) * (*t_-3*minT/4) + c[ 5] * pow(minT,4) * (*t_-4*minT/5) + c[ 6]) // else if out of bounds - low
                                       (             c[8] * *t_ + c[9] * maxT * (*t_-maxT/2) + c[10] * maxT * maxT * (*t_-2*maxT/3) + c[11] * pow(maxT,3) * (*t_-3*maxT/4) + c[12] * pow(maxT,4) * (*t_-4*maxT/5) + c[13]); // else out of bounds - high
        break;
      case SHOMATE2:
        h <<= h + *massFracs_[n] * cond( *t_ <= c[0] && *t_ >= minT, c[ 6] + c[1]**t_*1e-3 + c[2]/2 * t2*1e-6 + c[ 3]/3 * t3*1e-9 + c[ 4]/4 * t4*1e-12 - c[ 5] * recipT*1e3 ) // if low temp
                                       ( *t_ >  c[0] && *t_ <= maxT, c[13] + c[8]**t_*1e-3 + c[9]/2 * t2*1e-6 + c[10]/3 * t3*1e-9 + c[11]/4 * t4*1e-12 - c[12] * recipT*1e3 )  // else if high temp
                                       ( *t_ < minT, c[1] * *t_*1e-3 + c[2] * minT*1e-3 * ( *t_*1e-3 - minT*1e-3/2 ) + c[ 3] * minT*1e-3 * minT*1e-3 * ( *t_*1e-3 - 2*minT*1e-3/3 ) + c[ 4] * pow(minT*1e-3,3) * ( *t_*1e-3 - 3*minT*1e-3/4 ) - c[ 5] * pow(minT*1e-3,-1) * (-*t_*1e-3 / minT*1e-3 + 2 ) + c[ 6] ) // else if out of bounds - low
                                       (             c[8] * *t_*1e-3 + c[9] * maxT*1e-3 * ( *t_*1e-3 - maxT*1e-3/2 ) + c[10] * maxT*1e-3 * maxT*1e-3 * ( *t_*1e-3 - 2*maxT*1e-3/3 ) + c[11] * pow(maxT*1e-3,3) * ( *t_*1e-3 - 3*maxT*1e-3/4 ) - c[12] * pow(maxT*1e-3,-1) * (-*t_*1e-3 / maxT*1e-3 + 2 ) + c[13] ); // else out of bounds - high
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
    c_[1] /= molecularWeight; // convert to mass basis
    c_[3] /= molecularWeight; // convert to mass basis
    break;
  case NASA2:
    for( std::vector<double>::iterator ic = c_.begin() + 1; ic!=c_.end(); ++ic){
      *ic *= Cantera::GasConstant / molecularWeight; // dimensionalize the coefficients to mass basis
    }
    break;
  case SHOMATE2:
    for( std::vector<double>::iterator ic = c_.begin() + 1; ic!=c_.end(); ++ic ){
      *ic *= 1e6 / molecularWeight; // scale the coefficients to keep units consistent on mass basis
    }
    break;
  }

  CanteraObjects::restore_gasmix(gasMix);
}

//--------------------------------------------------------------------

template< typename FieldT >
SpeciesEnthalpy<FieldT>::
~SpeciesEnthalpy()
{

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

  const FieldT& t2 =     *tPowers_[0]; // t^2
  const FieldT& t3 =     *tPowers_[1]; // t^3
  const FieldT& t4 =     *tPowers_[2]; // t^4
  const FieldT& t5 =     *tPowers_[3]; // t^5
  const FieldT& recipT = *tPowers_[4]; // t^-1

  switch (polyType_) {
    case SIMPLE:
      h <<= c_[1] + c_[3] * (*t_ - c_[0]);
      break;
    case NASA2:
      h <<= cond( *t_ <= c_[0] && *t_ >= minT_, c_[ 6] + c_[1] * *t_ + c_[2]/2 * t2 + c_[ 3]/3 * t3 + c_[ 4]/4 * t4 + c_[ 5]/5 * t5) // if low temp
                ( *t_ >  c_[0] && *t_ <= maxT_, c_[13] + c_[8] * *t_ + c_[9]/2 * t2 + c_[10]/3 * t3 + c_[11]/4 * t4 + c_[12]/5 * t5)  // else if high temp
                ( *t_ <  minT_                 , c_[1] * *t_ + c_[2] * minT_ * (*t_-minT_/2) + c_[ 3] * minT_ * minT_ * (*t_ - 2*minT_/3) + c_[ 4] * pow(minT_,3) * (*t_ - 3*minT_/4) + c_[ 5] * pow(minT_,4) * (*t_ - 4*minT_/5) + c_[ 6]) // else if out of bounds - low
                (                                c_[8] * *t_ + c_[9] * maxT_ * (*t_-maxT_/2) + c_[10] * maxT_ * maxT_ * (*t_ - 2*maxT_/3) + c_[11] * pow(maxT_,3) * (*t_ - 3*maxT_/4) + c_[12] * pow(maxT_,4) * (*t_ - 4*maxT_/5) + c_[13]); // else out of bounds - high
      break;
    case SHOMATE2:
      h <<= cond( *t_ <= c_[0] && *t_ >= minT_, c_[ 6] + c_[1] * *t_*1e-3 + c_[2]/2 * t2*1e-6 + c_[ 3]/3 * t3*1e-9 + c_[ 4]/4 * t4*1e-12 - c_[ 5] * recipT*1e3 ) // if low temp
                ( *t_ >  c_[0] && *t_ <= maxT_, c_[13] + c_[8] * *t_*1e-3 + c_[9]/2 * t2*1e-6 + c_[10]/3 * t3*1e-9 + c_[11]/4 * t4*1e-12 - c_[12] * recipT*1e3 )  // else if high range
                ( *t_ <  minT_, c_[1] * *t_*1e-3 + c_[2] * minT_*1e-3 * ( *t_*1e-3 - minT_*1e-3/2 ) + c_[ 3] * minT_*1e-3 * minT_*1e-3 * ( *t_*1e-3 - 2*minT_*1e-3/3 ) + c_[ 4] * pow(minT_*1e-3,3) * ( *t_*1e-3 - 3*minT_*1e-3/4 ) - c_[ 5] * pow(minT_*1e-3,-1) * ( -*t_*1e-3/minT_*1e-3 + 2 ) + c_[ 6] ) // else if out of bounds - low
                (               c_[8] * *t_*1e-3 + c_[9] * maxT_*1e-3 * ( *t_*1e-3 - maxT_*1e-3/2 ) + c_[10] * maxT_*1e-3 * maxT_*1e-3 * ( *t_*1e-3 - 2*maxT_*1e-3/3 ) + c_[11] * pow(maxT_*1e-3,3) * ( *t_*1e-3 - 3*maxT_*1e-3/4 ) - c_[12] * pow(maxT_*1e-3,-1) * ( -*t_*1e-3/maxT_*1e-3 + 2 ) + c_[13] ); // else out of bounds - high
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
