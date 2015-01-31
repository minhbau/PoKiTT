#ifndef InternalEnergy_Expr_h
#define InternalEnergy_Expr_h

#include <expression/Expression.h>

#include <pokitt/CanteraObjects.h> //include cantera wrapper

#include <cantera/kernel/ct_defs.h> // contains value of gas constant
#include <cantera/kernel/SpeciesThermoInterpType.h> // contains definitions for which polynomial is being used

namespace pokitt{

/**
 *  \class  InternalEnergy
 *  \author Nathan Yonkee
 *  \date February, 2015
 *
 *  \brief Calculates the internal energy of a mixture
 *  using either NASA7 or Shomate polynomials with 2 temperature ranges.
 *  Units of J/kg
 *
 *  Five coefficients \f$(a_0,\dots,a_4)\f$ are used to represent
 * \f$ e^0(T)\f$ as a polynomial in \f$ T \f$ e.g. NASA7:
 *
 * \f[
 * \frac{e^0(T)}{R} = a_0 T + (a_1-1) T^2/2 + a_2 T^3/3 + a_3 T^4/4 + a_4 T^5/5 + a_5
 * \f]
 *
 * The InternalEnergy is a weighted average of the pure species \f$ e^0(T)\f$,
 *
 * \f[
 * e(T) = \sum_{n=0}^{nSpec} y_n e^0_n
 * \f]
 */

template< typename FieldT >
class InternalEnergy
    : public Expr::Expression<FieldT>
{
  typedef std::vector<double> PolyVals; // values used for polynomial
  DECLARE_FIELD( FieldT, t_ )
  DECLARE_VECTOR_OF_FIELDS( FieldT, massFracs_ )

  int nSpec_; // number of species
  PolyVals minTVec_; // vector of minimum temperatures for polynomial evaluations
  PolyVals maxTVec_; // vector of maximum temperatures for polynomial evaluations
  std::vector< PolyVals > cVec_; // vector of polynomial coefficients
  std::vector<int> polyTypeVec_; // vector of polynomial types
  bool shomateFlag_; // true if any polynomial is shomate

  InternalEnergy( const Expr::Tag& tTag,
                  const Expr::TagList& massFracTags );
public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a InternalEnergy expression
     *  @param resultTag tag for the mixture internal energy
     *  @param tTag tag for temperature
     *  @param massFracTags tags for the mass fractions: must be consistent with Cantera
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::Tag& tTag,
             const Expr::TagList& massFracTags );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag tTag_;
    const Expr::TagList massFracTags_;
  };

  ~InternalEnergy(){}
  void evaluate();
};

/**
 *  \class  SpeciesInternalEnergy
 *  \author Nathan Yonkee
 *  \date February, 2015
 *
 *  \brief Calculates the internal energy of each species
 *  using either NASA7 or Shomate polynomials with 2 temperature ranges.
 *  Units of J/kg
 *
 *  Five coefficients \f$(a_0,\dots,a_4)\f$ are used to represent
 * \f$ e(T)\f$ as a polynomial in \f$ T \f$ e.g. NASA7:
 *
 * \f[
 * \frac{e(T)}{R} = a_0 T + (a_1-1) T^2/2 + a_2 T^3/3 + a_3 T^4/4 + a_4 T^5/5 + a_5
 * \f]
 *
 */
template< typename FieldT >
class SpeciesInternalEnergy
 : public Expr::Expression<FieldT>
{
  DECLARE_FIELD( FieldT, t_ )

  const int n_; //index of species to be evaluated
  double minT_; // minimum temperature for polynomial evaluation
  double maxT_; // maximum temperature for polynomial evaluation
  std::vector<double> c_; // vector of polynomial coefficients
  int polyType_; // polynomial type

  SpeciesInternalEnergy( const Expr::Tag& tTag,
                         const int n );
public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a SpeciesInternalEnergy expression
     *  @param resultTag tag for the pure species internal energy
     *  @param tTag tag for temperature
     *  @param n species index used by Cantera
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::Tag& tTag,
             const int n);

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag tTag_;
    const int n_;
  };

  ~SpeciesInternalEnergy(){}
  void evaluate();
};



// ###################################################################
//
//                          Implementation
//
// ###################################################################

template< typename FieldT >
InternalEnergy<FieldT>::
InternalEnergy( const Expr::Tag& tTag,
                const Expr::TagList& massFracTags )
  : Expr::Expression<FieldT>(),
    shomateFlag_( false )
{
  this->set_gpu_runnable( true );

  t_ = this->template create_field_request<FieldT>(tTag);
  this->template create_field_vector_request<FieldT>( massFracTags, massFracs_ );

  Cantera_CXX::IdealGasMix* const gasMix = CanteraObjects::get_gasmix();
  const Cantera::SpeciesThermo& spThermo = gasMix->speciesThermo();

  nSpec_ = gasMix->nSpecies();

  const std::vector<double> molecularWeights = gasMix->molecularWeights();
  std::vector<double> c(15,0); //vector of Cantera's coefficients
  int polyType; // type of polynomial - const_cp, shomate, or NASA
  double minT; // minimum temperature polynomial is valid
  double maxT; // maximum temperature polynomial is valid
  double refPressure;

  for( size_t n=0; n<nSpec_; ++n ){
    spThermo.reportParams(n, polyType, &c[0], minT, maxT, refPressure);
    polyTypeVec_.push_back(polyType); // vector of polynomial types
    minTVec_.push_back(minT); // vector of minimum temperatures for polynomial evaluations
    maxTVec_.push_back(maxT); // vector of maximum temperatures for polynomial evaluations
    switch (polyType) {
    case SIMPLE:
      c[3] -= Cantera::GasConstant; // change coefficients from isobaric to isometric
      c[1] -= Cantera::GasConstant * c[0];
      c[1] /= molecularWeights[n]; // convert to mass basis
      c[3] /= molecularWeights[n]; // convert to mass basis
      break;
    case NASA2:
      c[1] -= 1.0; //change coefficients from isobaric to isometric
      c[8] -= 1.0;
      for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic){
        *ic *= Cantera::GasConstant / molecularWeights[n]; // dimensionalize the coefficients to mass basis
      }
      c[ 2] /= 2;
      c[ 3] /= 3;
      c[ 4] /= 4;
      c[ 5] /= 5;
      c[ 9] /= 2;
      c[10] /= 3;
      c[11] /= 4;
      c[12] /= 5;
      break;
    case SHOMATE2:
      c[1] -= Cantera::GasConstant * 1e-3; //change coefficients from isobaric to isometric
      c[8] -= Cantera::GasConstant * 1e-3;
      shomateFlag_ = true;
      for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic ){
        *ic *= 1e6 / molecularWeights[n]; // scale the coefficients to keep units consistent on mass basis
      }
      c[ 2] /= 2;
      c[ 3] /= 3;
      c[ 4] /= 4;
      c[ 9] /= 2;
      c[10] /= 3;
      c[11] /= 4;
      c[ 1] *= 1e-3;
      c[ 8] *= 1e-3;
      c[ 2] *= 1e-6;
      c[ 9] *= 1e-6;
      c[ 3] *= 1e-9;
      c[10] *= 1e-9;
      c[ 4] *= 1e-12;
      c[11] *= 1e-12;
      c[ 5] *= 1e3;
      c[12] *= 1e3;
      break;
    }
    cVec_.push_back(c); // vector of polynomial coefficients
  }

  CanteraObjects::restore_gasmix(gasMix);
}

//--------------------------------------------------------------------

template< typename FieldT >
void
InternalEnergy<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  FieldT& e = this->value();

  const FieldT& temp = t_->field_ref();

  SpatFldPtr<FieldT> recipT;      // may be used for Shomate polynomial
  if( shomateFlag_ ) {
    recipT = SpatialFieldStore::get<FieldT>(temp);
    *recipT <<= 1 / temp;
  }
# ifndef ENABLE_CUDA
  const double maxTval = nebo_max(temp);
  const double minTval = nebo_min(temp);
# endif

  e <<= 0.0;

  for( size_t n=0; n<nSpec_; ++n ){
    const FieldT& yi = massFracs_[n]->field_ref();
    int polyType = polyTypeVec_[n];
    std::vector<double>& c = cVec_[n];
    double minT = minTVec_[n];
    double maxT = maxTVec_[n];
#   ifndef ENABLE_CUDA // optimization benefits only the CPU - cond performs betters with if/else than with if/elif/elif/else
    if( maxTval <= maxT && minTval >= minT){ // if true, temperature can only be either high or low
      switch (polyType) {
      case SIMPLE: // constant heat capacity
        e <<= e + yi * ( c[1] + c[3] * (temp - c[0]) );
        break;
      case NASA2:
        e <<= e + yi * cond( temp <= c[0] , c[ 6] + temp * ( c[1] + temp * ( c[2] + temp * ( c[ 3] + temp * ( c[ 4] + temp * c[ 5] )))) )  // if low temp
                           (                c[13] + temp * ( c[8] + temp * ( c[9] + temp * ( c[10] + temp * ( c[11] + temp * c[12] )))) );  // else if high temp
        break;
      case SHOMATE2:
        e <<= e + yi * cond( temp <= c[0] , c[ 6] + temp * ( c[1] + temp * ( c[2] + temp * ( c[ 3] + temp * c[ 4] ))) - c[ 5] * *recipT ) // if low temp
                           (                c[13] + temp * ( c[8] + temp * ( c[9] + temp * ( c[10] + temp * c[11] ))) - c[12] * *recipT );  // else if high range
        break;
      }
    }
    else
#   endif
    {
      switch (polyType) {
      /* polynomial can be out of bounds low, low temp, high temp, or out of bounds high
       * if out of bounds, internal energy is interpolated from min or max temp using a constant cv
       */
      case SIMPLE: // constant heat capacity
        e <<= e + yi * ( c[1] + c[3] * (temp - c[0]) );
        break;
      case NASA2:
        e <<= e + yi * cond( temp <= c[0] && temp >= minT, c[ 6] + temp * ( c[1] + temp * ( c[2] + temp * ( c[ 3] + temp * ( c[ 4] + temp * c[ 5] )))) )  // if low temp
                           ( temp >  c[0] && temp <= maxT, c[13] + temp * ( c[8] + temp * ( c[9] + temp * ( c[10] + temp * ( c[11] + temp * c[12] )))) )  // else if high temp
                           ( temp < minT,                  c[ 6] + c[1] * temp + minT * ( 2*c[2] * temp + minT * ( 3*c[ 3] * temp - c[2] + minT * ( 4*c[ 4] * temp - 2*c[ 3] + minT * ( 5*c[ 5] * temp - 3*c[ 4] + minT * -4*c[ 5] )))) )  // else if out of bounds - low
                           (                               c[13] + c[8] * temp + maxT * ( 2*c[9] * temp + maxT * ( 3*c[10] * temp - c[9] + maxT * ( 4*c[11] * temp - 2*c[10] + maxT * ( 5*c[12] * temp - 3*c[11] + maxT * -4*c[12] )))) ); // else out of bounds - high
        break;
      case SHOMATE2:
        e <<= e + yi * cond( temp <= c[0] && temp >= minT, c[ 6] + temp * ( c[1] + temp * ( c[2] + temp * ( c[ 3] + temp * c[ 4] ))) - c[ 5] * *recipT ) // if low temp
                           ( temp >  c[0] && temp <= maxT, c[13] + temp * ( c[8] + temp * ( c[9] + temp * ( c[10] + temp * c[11] ))) - c[12] * *recipT )  // else if high range
                           ( temp <  minT,                 c[ 6] + c[1] * temp + minT * ( 2*c[2] * temp + minT * ( 3*c[ 3] * temp - c[2] + minT * ( 4*c[ 4] * temp - 2*c[ 3] + minT * -3*c[ 4] ))) + ( c[ 5] * temp / minT - 2*c[ 5] ) / minT ) // else if out of bounds - low
                           (                               c[13] + c[8] * temp + maxT * ( 2*c[9] * temp + maxT * ( 3*c[10] * temp - c[9] + maxT * ( 4*c[11] * temp - 2*c[10] + maxT * -3*c[11] ))) + ( c[12] * temp / maxT - 2*c[12] ) / maxT ); // else out of bounds - high
        break;
      }
    }
  }
}
//--------------------------------------------------------------------

template< typename FieldT >
InternalEnergy<FieldT>::
Builder::Builder( const Expr::Tag& resultTag,
                  const Expr::Tag& tTag,
                  const Expr::TagList& massFracTags )
: ExpressionBuilder( resultTag ),
  tTag_( tTag ),
  massFracTags_( massFracTags )
{}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
InternalEnergy<FieldT>::
Builder::build() const
{
  return new InternalEnergy<FieldT>( tTag_, massFracTags_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
SpeciesInternalEnergy<FieldT>::
SpeciesInternalEnergy( const Expr::Tag& tTag,
                       const int n )
  : Expr::Expression<FieldT>(),
    n_ ( n )
{
  this->set_gpu_runnable( true );

  t_ = this->template create_field_request<FieldT>( tTag );

  Cantera_CXX::IdealGasMix* const gasMix = CanteraObjects::get_gasmix();
  const Cantera::SpeciesThermo& spThermo = gasMix->speciesThermo();

  double molecularWeight = gasMix->molecularWeight(n_);
  c_.resize(15); //vector of Cantera's coefficients
  double refPressure;
  spThermo.reportParams(n_, polyType_, &c_[0], minT_, maxT_, refPressure);
  switch (polyType_) {
  case SIMPLE:
    c_[3] -= Cantera::GasConstant; // change coefficients from isobaric to isometric
    c_[1] -= Cantera::GasConstant * c_[0];
    c_[1] /= molecularWeight; // convert to mass basis
    c_[3] /= molecularWeight; // convert to mass basis
    break;
  case NASA2:
    c_[1] -= 1.0; //change coefficients from isobaric to isometric
    c_[8] -= 1.0;
    for( std::vector<double>::iterator ic = c_.begin() + 1; ic!=c_.end(); ++ic){
      *ic *= Cantera::GasConstant / molecularWeight; // dimensionalize the coefficients to mass basis
    }
    c_[2] /= 2;
    c_[3] /= 3;
    c_[4] /= 4;
    c_[5] /= 5;
    c_[9] /= 2;
    c_[10] /= 3;
    c_[11] /= 4;
    c_[12] /= 5;
    break;
  case SHOMATE2:
    c_[1] -= Cantera::GasConstant * 1e-3; //change coefficients from isobaric to isometric
    c_[8] -= Cantera::GasConstant * 1e-3;
    for( std::vector<double>::iterator ic = c_.begin() + 1; ic!=c_.end(); ++ic ){
      *ic *= 1e6 / molecularWeight; // scale the coefficients to keep units consistent on mass basis
    }
    c_[ 2] /= 2;
    c_[ 3] /= 3;
    c_[ 4] /= 4;
    c_[ 9] /= 2;
    c_[10] /= 3;
    c_[11] /= 4;
    c_[ 1] *= 1e-3;
    c_[ 8] *= 1e-3;
    c_[ 2] *= 1e-6;
    c_[ 9] *= 1e-6;
    c_[ 3] *= 1e-9;
    c_[10] *= 1e-9;
    c_[ 4] *= 1e-12;
    c_[11] *= 1e-12;
    c_[ 5] *= 1e3;
    c_[12] *= 1e3;
    break;
  }

  CanteraObjects::restore_gasmix(gasMix);
}

//--------------------------------------------------------------------

template< typename FieldT >
void
SpeciesInternalEnergy<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  FieldT& e = this->value();

  const FieldT& temp = t_->field_ref();

  switch (polyType_) {
  /* polynomial can be out of bounds low, low temp, high temp, or out of bounds high
   * if out of bounds, internal energy is interpolated from min or max temp using a constant cv
   */
  case SIMPLE:
    e <<= c_[1] + c_[3] * (temp - c_[0]);
    break;
  case NASA2:
    e <<= cond( temp <= c_[0] && temp >= minT_, c_[ 6] + temp * ( c_[1] + temp * ( c_[2] + temp * ( c_[ 3] + temp * ( c_[ 4] + temp * c_[ 5] )))) )  // if low temp
              ( temp >  c_[0] && temp <= maxT_, c_[13] + temp * ( c_[8] + temp * ( c_[9] + temp * ( c_[10] + temp * ( c_[11] + temp * c_[12] )))) )  // else if high temp
              ( temp < minT_,                   c_[ 6] + c_[1] * temp + minT_ * ( 2*c_[2] * temp + minT_ * ( 3*c_[ 3] * temp - c_[2] + minT_ * ( 4*c_[ 4] * temp - 2*c_[ 3] + minT_ * ( 5*c_[ 5] * temp - 3*c_[ 4] + minT_ * -4*c_[ 5] )))) )  // else if out of bounds - low
              (                                 c_[13] + c_[8] * temp + maxT_ * ( 2*c_[9] * temp + maxT_ * ( 3*c_[10] * temp - c_[9] + maxT_ * ( 4*c_[11] * temp - 2*c_[10] + maxT_ * ( 5*c_[12] * temp - 3*c_[11] + maxT_ * -4*c_[12] )))) ); // else out of bounds - high
    break;
  case SHOMATE2:
    e <<= cond( temp <= c_[0] && temp >= minT_, c_[ 6] + temp * ( c_[1] + temp * ( c_[2] + temp * ( c_[ 3] + temp * c_[ 4] ))) - c_[ 5] / temp ) // if low temp
              ( temp >  c_[0] && temp <= maxT_, c_[13] + temp * ( c_[8] + temp * ( c_[9] + temp * ( c_[10] + temp * c_[11] ))) - c_[12] / temp )  // else if high range
              ( temp <  minT_,                  c_[ 6] + c_[1] * temp + minT_ * ( 2*c_[2] * temp + minT_ * ( 3*c_[ 3] * temp - c_[2] + minT_ * ( 4*c_[ 4] * temp - 2*c_[ 3] + minT_ * -3*c_[ 4] ))) + ( c_[ 5] * temp / minT_ - 2*c_[ 5] ) / minT_ ) // else if out of bounds - low
              (                                 c_[13] + c_[8] * temp + maxT_ * ( 2*c_[9] * temp + maxT_ * ( 3*c_[10] * temp - c_[9] + maxT_ * ( 4*c_[11] * temp - 2*c_[10] + maxT_ * -3*c_[11] ))) + ( c_[12] * temp / maxT_ - 2*c_[12] ) / maxT_ ); // else out of bounds - high
    break;
  }
}

//--------------------------------------------------------------------

template< typename FieldT >
SpeciesInternalEnergy<FieldT>::
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
SpeciesInternalEnergy<FieldT>::
Builder::build() const
{
  return new SpeciesInternalEnergy<FieldT>( tTag_, n_ );
}

} // namespace pokitt

#endif // InternalEnergy_Expr_h