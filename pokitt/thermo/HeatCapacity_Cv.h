#ifndef HeatCapacity_Cv_Expr_h
#define HeatCapacity_Cv_Expr_h

#include <expression/Expression.h>

#include <pokitt/CanteraObjects.h> //include cantera wrapper

#include <cantera/kernel/ct_defs.h> // contains value of Cantera::GasConstant
#include <cantera/kernel/speciesThermoTypes.h> // contains definitions for which polynomial is being used

#include <exception>


namespace pokitt{

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
  DECLARE_FIELD( FieldT, t_ )
  DECLARE_VECTOR_OF_FIELDS( FieldT, massFracs_ )

  int nSpec_; // number of species
  PolyVals minTVec_; // vector of minimum temperatures for polynomial evaluations
  PolyVals maxTVec_; // vector of maximum temperatures for polynomial evaluations
  std::vector< PolyVals > cVec_; // vector of polynomial coefficients
  std::vector<int> polyTypeVec_; // vector of polynomial types
  bool shomateFlag_; // true if any polynomial is shomate

  HeatCapacity_Cv( const Expr::Tag& tTag,
                   const Expr::TagList& massFracTags );
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
             const Expr::TagList& massFracTags );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag tTag_;
    const Expr::TagList massFracTags_;
  };

  ~HeatCapacity_Cv(){}
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
  DECLARE_FIELD( FieldT, t_ )

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

  ~SpeciesHeatCapacity_Cv(){}
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
                 const Expr::TagList& massFracTags )
  : Expr::Expression<FieldT>(),
    shomateFlag_( false )
{
  this->set_gpu_runnable( true );

  t_ = this->template create_field_request<FieldT>( tTag );
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
      c[3] /= molecularWeights[n]; // convert to mass basis
      break;
    case NASA2:
      c[1] -= 1.0; //change coefficients from isobaric to isometric
      c[8] -= 1.0;
      for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic)
        *ic *= Cantera::GasConstant / molecularWeights[n]; // dimensionalize the coefficients to mass basis
      break;
    case SHOMATE2:
      shomateFlag_ = true;
      c[1] -= Cantera::GasConstant * 1e-3; //change coefficients from isobaric to isometric
      c[8] -= Cantera::GasConstant * 1e-3;
      for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic ){
        *ic *= 1e3 / molecularWeights[n]; // scale the coefficients to keep units consistent on mass basis
      }
      //Shomate polynomial uses T/1000 so we multiply coefficients by 1e-3 to avoid division
      c[ 2] *= 1e-3;
      c[ 9] *= 1e-3;
      c[ 3] *= 1e-6;
      c[10] *= 1e-6;
      c[ 4] *= 1e-9;
      c[11] *= 1e-9;
      c[ 5] *= 1e6;
      c[12] *= 1e6;
      break;
    }
    cVec_.push_back(c); // vector of polynomial coefficients
  }

  CanteraObjects::restore_gasmix(gasMix);
}

//--------------------------------------------------------------------

template< typename FieldT >
void
HeatCapacity_Cv<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  FieldT& cv = this->value();

  const FieldT& temp = t_->field_ref();

  SpatFldPtr<FieldT> recipRecipT;      // may be used for Shomate polynomial
  if( shomateFlag_ ) {
    recipRecipT = SpatialFieldStore::get<FieldT>(temp);
    *recipRecipT <<= 1 / ( temp * temp );
  }
# ifndef ENABLE_CUDA
  const double maxTval = nebo_max(temp);
  const double minTval = nebo_min(temp);
# endif

  cv <<= 0.0; // set cv to 0 before starting summation

  for( size_t n=0; n<nSpec_; ++n ){
    const FieldT& yi = massFracs_[n]->field_ref();
    const int polyType = polyTypeVec_[n];
    const std::vector<double>& c = cVec_[n];
    const double minT = minTVec_[n];
    const double maxT = maxTVec_[n];
#   ifndef ENABLE_CUDA // optimization benefits only the CPU - cond performs betters with if/else than with if/elif/elif/else
    if( maxTval <= maxT && minTval >= minT){ // if true, temperature can only be either high or low
      switch (polyType) {
      case SIMPLE:
        cv <<= cv + yi * c[3];
        break;
      case NASA2:
        /* polynomials are applicable in two temperature ranges - high and low
         * If the temperature is out of range, the value is set to the value at the min or max temp
         */
        cv <<= cv + yi * cond( temp <= c[0] , c[1] + temp * ( c[2] + temp * ( c[ 3] + temp * ( c[ 4] + temp * c[ 5] ))) )  // if low temp
                             (                c[8] + temp * ( c[9] + temp * ( c[10] + temp * ( c[11] + temp * c[12] ))) );  // else if high temp
        break;
      case SHOMATE2:
        cv <<= cv + yi * cond( temp <= c[0] , c[1] + temp * ( c[2] + temp * ( c[ 3] + temp * c[ 4])) + c[ 5] * *recipRecipT )  // if low temp
                             (                c[8] + temp * ( c[9] + temp * ( c[10] + temp * c[11])) + c[12] * *recipRecipT );  // else if high temp
        break;
      }
    }
    else
#   endif
    {
      switch (polyType) {
      case SIMPLE:
        cv <<= cv + yi * c[3];
        break;
      case NASA2:
        /* polynomials are applicable in two temperature ranges - high and low
         * If the temperature is out of range, the value is set to the value at the min or max temp
         */
        cv <<= cv + yi * cond( temp <= c[0] && temp >= minT, c[1] + temp * ( c[2] + temp * ( c[ 3] + temp * ( c[ 4] + temp * c[ 5] ))) )  // if low temp
                             ( temp >  c[0] && temp <= maxT, c[8] + temp * ( c[9] + temp * ( c[10] + temp * ( c[11] + temp * c[12] ))) )  // else if high temp
                             ( temp < minT,                  c[1] + minT * ( c[2] + minT * ( c[ 3] + minT * ( c[ 4] + minT * c[ 5] ))) )  // else if out of bounds - low
                             (                               c[8] + maxT * ( c[9] + maxT * ( c[10] + maxT * ( c[11] + maxT * c[12] ))) ); // else out of bounds - high
        break;
      case SHOMATE2:
        cv <<= cv + yi * cond( temp <= c[0] && temp >= minT, c[1] + temp * ( c[2] + temp * ( c[ 3] + temp * c[ 4])) + c[ 5] * *recipRecipT )  // if low temp
                             ( temp >  c[0] && temp <= maxT, c[8] + temp * ( c[9] + temp * ( c[10] + temp * c[11])) + c[12] * *recipRecipT )  // else if high temp
                             ( temp < minT,                  c[1] + minT * ( c[2] + minT * ( c[ 3] + minT * c[ 4])) + c[ 5] / (minT*minT)  )  // else if out of bounds - low
                             (                               c[8] + maxT * ( c[9] + maxT * ( c[10] + maxT * c[11])) + c[12] / (maxT*maxT)  ); // else out of bounds - high
        break;
      }
    }
  }
}

//--------------------------------------------------------------------

template< typename FieldT >
HeatCapacity_Cv<FieldT>::
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
HeatCapacity_Cv<FieldT>::
Builder::build() const
{
  return new HeatCapacity_Cv<FieldT>( tTag_, massFracTags_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
SpeciesHeatCapacity_Cv<FieldT>::
SpeciesHeatCapacity_Cv( const Expr::Tag& tTag,
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
    //Shomate polynomial uses T/1000 so we multiply coefficients by 1e-3 to avoid division
    c_[ 2] *= 1e-3;
    c_[ 9] *= 1e-3;
    c_[ 3] *= 1e-6;
    c_[10] *= 1e-6;
    c_[ 4] *= 1e-9;
    c_[11] *= 1e-9;
    c_[ 5] *= 1e6;
    c_[12] *= 1e6;
    break;
  }

  CanteraObjects::restore_gasmix(gasMix);
}

//--------------------------------------------------------------------

template< typename FieldT >
void
SpeciesHeatCapacity_Cv<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  FieldT& cv = this->value();

  const FieldT& temp = t_->field_ref();

  switch (polyType_) {
  case SIMPLE:
    cv <<= c_[3];
    break;
  case NASA2:
    /* polynomials are applicable in two temperature ranges - high and low
     * If the temperature is out of range, the value is set to the value at the min or max temp
     */
    cv <<= cond( temp <= c_[0] && temp >= minT_, c_[1] + temp  * ( c_[2] + temp  * ( c_[ 3] + temp  * ( c_[ 4] + temp  * c_[ 5] ))) )  // if low temp
               ( temp >  c_[0] && temp <= maxT_, c_[8] + temp  * ( c_[9] + temp  * ( c_[10] + temp  * ( c_[11] + temp  * c_[12] ))) )  // else if high temp
               ( temp < minT_,                   c_[1] + minT_ * ( c_[2] + minT_ * ( c_[ 3] + minT_ * ( c_[ 4] + minT_ * c_[ 5] ))) )  // else if out of bounds - low
               (                                 c_[8] + maxT_ * ( c_[9] + maxT_ * ( c_[10] + maxT_ * ( c_[11] + maxT_ * c_[12] ))) ); // else out of bounds - high
    break;
  case SHOMATE2:
    cv <<= cond( temp <= c_[0] && temp >= minT_, c_[1] + temp  * ( c_[2] + temp  * ( c_[ 3] + temp  * c_[ 4])) + c_[ 5] / ( temp * temp ) )  // if low temp
               ( temp >  c_[0] && temp <= maxT_, c_[8] + temp  * ( c_[9] + temp  * ( c_[10] + temp  * c_[11])) + c_[12] / ( temp * temp ) )  // else if high temp
               ( temp < minT_,                   c_[1] + minT_ * ( c_[2] + minT_ * ( c_[ 3] + minT_ * c_[ 4])) + c_[ 5] / (minT_*minT_) )  // else if out of bounds - low
               (                                 c_[8] + maxT_ * ( c_[9] + maxT_ * ( c_[10] + maxT_ * c_[11])) + c_[12] / (maxT_*maxT_) ); // else out of bounds - high
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

} // namespace pokitt

#endif // HeatCapacity_Cv_Expr_h
