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

#ifndef Enthalpy_Expr_h
#define Enthalpy_Expr_h

#include <expression/Expression.h>

#include <pokitt/CanteraObjects.h> //include cantera wrapper

namespace pokitt{

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
  typedef std::vector<double> PolyVals; // values used for polynomial
  DECLARE_FIELD( FieldT, t_ )
  DECLARE_VECTOR_OF_FIELDS( FieldT, massFracs_ )

  const int nSpec_; // number of species
  bool shomateFlag_; // true if any polynomial is shomate
  std::vector< ThermData > specThermVec_;
  std::ostringstream exceptionMsg_; // generic exception to be thrown

  Enthalpy( const Expr::Tag& tTag,
            const Expr::TagList& massFracTags );
public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a Enthalpy expression
     *  @param resultTag tag for the mixture enthalpy
     *  @param tTag tag for temperature
     *  @param massFracTags mass fractions of each species, ordering is consistent with Cantera input
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::Tag& tTag,
             const Expr::TagList& massFracTags );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag tTag_;
    const Expr::TagList massFracTags_;
  };

  ~Enthalpy(){}
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
  DECLARE_FIELD( FieldT, t_ )

  const int n_; //index of species to be evaluated
  ThermData specTherm_;
  std::ostringstream exceptionMsg_; // generic exception to be thrown

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

  ~SpeciesEnthalpy(){}
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
          const Expr::TagList& massFracTags )
  : Expr::Expression<FieldT>(),
    shomateFlag_( false ),
    nSpec_( CanteraObjects::number_species() )
{
  this->set_gpu_runnable( true );

  exceptionMsg_ << "\nUnidentified polynomial type somehow evaded detection\n"
                << "This should be have been caught in Cantera Objects\n";

  t_ = this->template create_field_request<FieldT>(tTag);
  this->template create_field_vector_request<FieldT>( massFracTags, massFracs_ );

  const std::vector<double>& molecularWeights = CanteraObjects::molecular_weights();
  const double gasConstant = CanteraObjects::gas_constant();

  for( size_t n=0; n<nSpec_; ++n ){
    ThermData tData = CanteraObjects::species_thermo( n );
    std::vector<double>& c = tData.coefficients;
    const ThermoPoly type = tData.type;
    switch ( type ) {
    case CONST_POLY:
      c[1] /= molecularWeights[n]; // convert to mass basis
      c[3] /= molecularWeights[n]; // convert to mass basis
      break;
    case NASA_POLY:
      for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic){
        *ic *= gasConstant / molecularWeights[n]; // dimensionalize the coefficients to mass basis
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
    case SHOMATE_POLY:
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
    default: {
      std::ostringstream msg;
      msg << __FILE__ << " : " << __LINE__ << "\n Error for spec n = " << n << exceptionMsg_.str();
      throw std::runtime_error( msg.str() );
      }
    }
    specThermVec_.push_back( tData );
  }

}

//--------------------------------------------------------------------

template< typename FieldT >
void
Enthalpy<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  FieldT& h = this->value();

  const FieldT& temp = t_->field_ref();

  SpatFldPtr<FieldT> recipT;      // may be used for Shomate polynomial
  if( shomateFlag_ ) {
    recipT = SpatialFieldStore::get<FieldT>(temp);
    *recipT <<= 1 / temp;
  }
# ifndef ENABLE_CUDA
  const double maxTval = field_max_interior(temp);
  const double minTval = field_min_interior(temp);
# endif

  h <<= 0.0;

  for( size_t n=0; n<nSpec_; ++n ){
    const FieldT& yi = massFracs_[n]->field_ref();
    const ThermData& thermo = specThermVec_[n];
    const ThermoPoly polyType = thermo.type;
    const std::vector<double>& c = thermo.coefficients;
    const double minT = thermo.minTemp;
    const double maxT = thermo.maxTemp;
#   ifndef ENABLE_CUDA // optimization benefits only the CPU - cond performs betters with if/else than with if/elif/elif/else
    if( maxTval <= maxT && minTval >= minT){ // if true, temperature can only be either high or low
      switch (polyType) {
      /* polynomial can be out of bounds low, low temp, high temp, or out of bounds high
       * if out of bounds, enthalpy is interpolated from min or max temp using a constant cp
       */
      case CONST_POLY: // constant heat capacity
        h <<= h + yi * ( c[1] + c[3] * (temp - c[0]) );
        break;
      case NASA_POLY:
        h <<= h + yi * cond( temp <= c[0] , c[ 6] + temp * ( c[1] + temp * ( c[2] + temp * ( c[ 3] + temp * ( c[ 4] + temp * c[ 5] )))) )  // if low temp
                           (                c[13] + temp * ( c[8] + temp * ( c[9] + temp * ( c[10] + temp * ( c[11] + temp * c[12] )))) );  // else if high temp
        break;
      case SHOMATE_POLY:
        h <<= h + yi * cond( temp <= c[0] , c[ 6] + temp * ( c[1] + temp * ( c[2] + temp * ( c[ 3] + temp * c[ 4] ))) - c[ 5] * *recipT ) // if low temp
                           (                c[13] + temp * ( c[8] + temp * ( c[9] + temp * ( c[10] + temp * c[11] ))) - c[12] * *recipT );  // else if high range
        break;
      default: {
        std::ostringstream msg;
        msg << __FILE__ << " : " << __LINE__ << "\n Error for spec n = " << n << exceptionMsg_.str();
        throw std::runtime_error( msg.str() );
        }
      }
    }
    else
#   endif
    {
      switch (polyType) {
      /* polynomial can be out of bounds low, low temp, high temp, or out of bounds high
       * if out of bounds, enthalpy is interpolated from min or max temp using a constant cp
       */
      case CONST_POLY: // constant heat capacity
        h <<= h + yi * ( c[1] + c[3] * (temp - c[0]) );
        break;
      case NASA_POLY:
        h <<= h + yi * cond( temp <= c[0] && temp >= minT, c[ 6] + temp * ( c[1] + temp * ( c[2] + temp * ( c[ 3] + temp * ( c[ 4] + temp * c[ 5] )))) )  // if low temp
                           ( temp >  c[0] && temp <= maxT, c[13] + temp * ( c[8] + temp * ( c[9] + temp * ( c[10] + temp * ( c[11] + temp * c[12] )))) )  // else if high temp
                           ( temp < minT,                  c[ 6] + c[1] * temp + minT * ( 2*c[2] * temp + minT * ( 3*c[ 3] * temp - c[2] + minT * ( 4*c[ 4] * temp - 2*c[ 3] + minT * ( 5*c[ 5] * temp - 3*c[ 4] + minT * -4*c[ 5] )))) )  // else if out of bounds - low
                           (                               c[13] + c[8] * temp + maxT * ( 2*c[9] * temp + maxT * ( 3*c[10] * temp - c[9] + maxT * ( 4*c[11] * temp - 2*c[10] + maxT * ( 5*c[12] * temp - 3*c[11] + maxT * -4*c[12] )))) ); // else out of bounds - high
        break;
      case SHOMATE_POLY:
        h <<= h + yi * cond( temp <= c[0] && temp >= minT, c[ 6] + temp * ( c[1] + temp * ( c[2] + temp * ( c[ 3] + temp * c[ 4] ))) - c[ 5] * *recipT ) // if low temp
                           ( temp >  c[0] && temp <= maxT, c[13] + temp * ( c[8] + temp * ( c[9] + temp * ( c[10] + temp * c[11] ))) - c[12] * *recipT )  // else if high range
                           ( temp <  minT,                 c[ 6] + c[1] * temp + minT * ( 2*c[2] * temp + minT * ( 3*c[ 3] * temp - c[2] + minT * ( 4*c[ 4] * temp - 2*c[ 3] + minT * -3*c[ 4] ))) + ( c[ 5] * temp / minT - 2*c[ 5] ) / minT ) // else if out of bounds - low
                           (                               c[13] + c[8] * temp + maxT * ( 2*c[9] * temp + maxT * ( 3*c[10] * temp - c[9] + maxT * ( 4*c[11] * temp - 2*c[10] + maxT * -3*c[11] ))) + ( c[12] * temp / maxT - 2*c[12] ) / maxT ); // else out of bounds - high
        break;
      default: {
        std::ostringstream msg;
        msg << __FILE__ << " : " << __LINE__ << "\n Error for spec n = " << n << exceptionMsg_.str();
        throw std::runtime_error( msg.str() );
        }
      }
    }
  }
}
//--------------------------------------------------------------------

template< typename FieldT >
Enthalpy<FieldT>::
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
Enthalpy<FieldT>::
Builder::build() const
{
  return new Enthalpy<FieldT>( tTag_, massFracTags_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
SpeciesEnthalpy<FieldT>::
SpeciesEnthalpy( const Expr::Tag& tTag,
                 const int n )
  : Expr::Expression<FieldT>(),
    n_ ( n ),
    specTherm_( CanteraObjects::species_thermo( n ) )
{
  this->set_gpu_runnable( true );

  exceptionMsg_ << "\nUnidentified polynomial type somehow evaded detection\n"
                << "This should be have been caught in Cantera Objects\n";

  t_ = this->template create_field_request<FieldT>( tTag );

  const std::vector<double>& molecularWeights = CanteraObjects::molecular_weights();
  const double molecularWeight = molecularWeights[n];
  const double gasConstant = CanteraObjects::gas_constant();
  std::vector<double>& c = specTherm_.coefficients;
  switch ( specTherm_.type ) {
  case CONST_POLY:
    c[1] /= molecularWeight; // convert to mass basis
    c[3] /= molecularWeight; // convert to mass basis
    break;
  case NASA_POLY:
    for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic){
      *ic *= gasConstant / molecularWeight; // dimensionalize the coefficients to mass basis
    }
    c[2] /= 2;
    c[3] /= 3;
    c[4] /= 4;
    c[5] /= 5;
    c[9] /= 2;
    c[10] /= 3;
    c[11] /= 4;
    c[12] /= 5;
    break;
  case SHOMATE_POLY:
    for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic ){
      *ic *= 1e6 / molecularWeight; // scale the coefficients to keep units consistent on mass basis
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
  default: {
    std::ostringstream msg;
    msg << __FILE__ << " : " << __LINE__ << "\n Error for spec n = " << n_ << exceptionMsg_.str();
    throw std::runtime_error( msg.str() );
    }
  }

}

//--------------------------------------------------------------------

template< typename FieldT >
void
SpeciesEnthalpy<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  FieldT& h = this->value();

  const FieldT& temp = t_->field_ref();

  const std::vector<double>& c = specTherm_.coefficients;
  const double minT = specTherm_.minTemp;
  const double maxT = specTherm_.maxTemp;
  switch ( specTherm_.type ) {
  /* polynomial can be out of bounds low, low temp, high temp, or out of bounds high
   * if out of bounds, enthalpy is interpolated from min or max temp using a constant cp
   */
  case CONST_POLY:
    h <<= c[1] + c[3] * (temp - c[0]);
    break;
  case NASA_POLY:
    h <<= cond( temp <= c[0] && temp >= minT, c[ 6] + temp * ( c[1] + temp * ( c[2] + temp * ( c[ 3] + temp * ( c[ 4] + temp * c[ 5] )))) )  // if low temp
              ( temp >  c[0] && temp <= maxT, c[13] + temp * ( c[8] + temp * ( c[9] + temp * ( c[10] + temp * ( c[11] + temp * c[12] )))) )  // else if high temp
              ( temp < minT,                  c[ 6] + c[1] * temp + minT * ( 2*c[2] * temp + minT * ( 3*c[ 3] * temp - c[2] + minT * ( 4*c[ 4] * temp - 2*c[ 3] + minT * ( 5*c[ 5] * temp - 3*c[ 4] + minT * -4*c[ 5] )))) )  // else if out of bounds - low
              (                               c[13] + c[8] * temp + maxT * ( 2*c[9] * temp + maxT * ( 3*c[10] * temp - c[9] + maxT * ( 4*c[11] * temp - 2*c[10] + maxT * ( 5*c[12] * temp - 3*c[11] + maxT * -4*c[12] )))) ); // else out of bounds - high
    break;
  case SHOMATE_POLY:
    h <<= cond( temp <= c[0] && temp >= minT, c[ 6] + temp * ( c[1] + temp * ( c[2] + temp * ( c[ 3] + temp * c[ 4] ))) - c[ 5] / temp ) // if low temp
              ( temp >  c[0] && temp <= maxT, c[13] + temp * ( c[8] + temp * ( c[9] + temp * ( c[10] + temp * c[11] ))) - c[12] / temp )  // else if high range
              ( temp <  minT,                 c[ 6] + c[1] * temp + minT * ( 2*c[2] * temp + minT * ( 3*c[ 3] * temp - c[2] + minT * ( 4*c[ 4] * temp - 2*c[ 3] + minT * -3*c[ 4] ))) + ( c[ 5] * temp / minT - 2*c[ 5] ) / minT ) // else if out of bounds - low
              (                               c[13] + c[8] * temp + maxT * ( 2*c[9] * temp + maxT * ( 3*c[10] * temp - c[9] + maxT * ( 4*c[11] * temp - 2*c[10] + maxT * -3*c[11] ))) + ( c[12] * temp / maxT - 2*c[12] ) / maxT ); // else out of bounds - high
    break;
  default: {
    std::ostringstream msg;
    msg << __FILE__ << " : " << __LINE__ << "\n Error for spec n = " << n_ << exceptionMsg_.str();
    throw std::runtime_error( msg.str() );
    }
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

} // namespace pokitt

#endif // Enthalpy_Expr_h
