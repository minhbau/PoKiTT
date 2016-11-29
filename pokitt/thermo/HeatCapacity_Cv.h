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

#ifndef HeatCapacity_Cv_Expr_h
#define HeatCapacity_Cv_Expr_h

#include <expression/Expression.h>

#include <pokitt/CanteraObjects.h> //include cantera wrapper

namespace pokitt{

/**
 *  \class  HeatCapacity_Cv
 *  \author Nate Yonkee
 *  \date May, 2014
 *
 *  \brief Calculates the constant volume heat capacity of each species
 *  using NASA7 polynomials with 2 temperature ranges.
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
 * The heat capacity is a weighted average of the pure species \f$ c_v^0(T)\f$,
 *
 * \f[
 * c_v(T) = \sum_{n=0}^{nSpec} Y_n c{v,n}
 * \f]
 */

template< typename FieldT >
class HeatCapacity_Cv
 : public Expr::Expression<FieldT>
{
  typedef std::vector<double> PolyVals; // values used for polynomial
  DECLARE_FIELD( FieldT, t_ )
  DECLARE_VECTOR_OF_FIELDS( FieldT, massFracs_ )

  const int nSpec_; // number of species
  std::vector< ThermData > specThermVec_;

  std::ostringstream exceptionMsg_; // generic exception to be thrown

  HeatCapacity_Cv( const Expr::Tag& tTag,
                   const Expr::TagList& massFracTags );
public:
  class Builder : public Expr::ExpressionBuilder
  {
    const Expr::Tag tTag_;
    const Expr::TagList massFracTags_;
  public:
    /**
     *  @brief Build a HeatCapacity_Cv expression
     *  @param resultTag the mixture heat capacity
     *  @param tTag tag for temperature
     *  @param massFracTags mass fractions of each species, ordering is consistent with Cantera input
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::Tag& tTag,
             const Expr::TagList& massFracTags,
             const SpatialOps::GhostData nghost = DEFAULT_NUMBER_OF_GHOSTS )
    : ExpressionBuilder( resultTag, nghost ),
      tTag_( tTag ),
      massFracTags_( massFracTags )
    {}

    Expr::ExpressionBase* build() const{
      return new HeatCapacity_Cv<FieldT>( tTag_, massFracTags_ );
    }
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
 *  using NASA7 polynomials with 2 temperature ranges.
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
 */

template< typename FieldT >
class SpeciesHeatCapacity_Cv
 : public Expr::Expression<FieldT>
{
  DECLARE_FIELD( FieldT, t_ )

  const int n_; //index of species to be evaluated
  ThermData specTherm_;

  std::ostringstream exceptionMsg_; // generic exception to be thrown

  SpeciesHeatCapacity_Cv( const Expr::Tag& tTag,
                          const int n );
public:
  class Builder : public Expr::ExpressionBuilder
  {
    const Expr::Tag tTag_;
    const int n_;
  public:
    /**
     *  @brief Build a HeatCapacity expression
     *  @param resultTag tag for the pure species heat capacity
     *  @param tTag tag for temperature
     *  @param n species index
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::Tag& tTag,
             const int n,
             const SpatialOps::GhostData nghost = DEFAULT_NUMBER_OF_GHOSTS )
    : ExpressionBuilder( resultTag, nghost ),
      tTag_( tTag ),
      n_( n )
    {}

    Expr::ExpressionBase* build() const{
      return new SpeciesHeatCapacity_Cv<FieldT>( tTag_, n_ );
    }
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
    nSpec_( CanteraObjects::number_species() )
{
  this->set_gpu_runnable( true );

  exceptionMsg_ << "\nUnidentified polynomial type somehow evaded detection\n"
                << "This should be have been caught in Cantera Objects\n";

  t_ = this->template create_field_request<FieldT>( tTag );
  this->template create_field_vector_request<FieldT>( massFracTags, massFracs_ );

  const std::vector<double>& molecularWeights = CanteraObjects::molecular_weights();
  const double gasConstant = CanteraObjects::gas_constant();

  for( size_t n=0; n<nSpec_; ++n ){
    ThermData tData = CanteraObjects::species_thermo( n );
    std::vector<double>& c = tData.coefficients;
    const ThermoPoly type = tData.type;
    switch ( type ) {
    case CONST_POLY:
      c[3] -= gasConstant; // change coefficients from isobaric to isometric
      c[3] /= molecularWeights[n]; // convert to mass basis
      break;
    case NASA_POLY:
      c[1] -= 1.0; //change coefficients from isobaric to isometric
      c[8] -= 1.0;
      for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic)
        *ic *= gasConstant / molecularWeights[n]; // dimensionalize the coefficients to mass basis
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
HeatCapacity_Cv<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  FieldT& cv = this->value();

  const FieldT& temp = t_->field_ref();

# ifndef ENABLE_CUDA
  const double maxTval = field_max_interior(temp);
  const double minTval = field_min_interior(temp);
# endif

  cv <<= 0.0; // set cv to 0 before starting summation

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
      case CONST_POLY:
        cv <<= cv + yi * c[3];
        break;
      case NASA_POLY:
        /* polynomials are applicable in two temperature ranges - high and low
         * If the temperature is out of range, the value is set to the value at the min or max temp
         */
        cv <<= cv + yi * cond( temp <= c[0] , c[8] + temp * ( c[9] + temp * ( c[10] + temp * ( c[11] + temp * c[12] ))) )  // if low temp
                             (                c[1] + temp * ( c[2] + temp * ( c[ 3] + temp * ( c[ 4] + temp * c[ 5] ))) ); // else if high temp
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
      case CONST_POLY:
        cv <<= cv + yi * c[3];
        break;
      case NASA_POLY:
        /* polynomials are applicable in two temperature ranges - high and low
         * If the temperature is out of range, the value is set to the value at the min or max temp
         */
        cv <<= cv + yi * cond( temp <= c[0] && temp >= minT, c[8] + temp * ( c[9] + temp * ( c[10] + temp * ( c[11] + temp * c[12] ))) )  // if low temp
                             ( temp >  c[0] && temp <= maxT, c[1] + temp * ( c[2] + temp * ( c[ 3] + temp * ( c[ 4] + temp * c[ 5] ))) )  // else if high temp
                             ( temp < minT,                  c[8] + minT * ( c[9] + minT * ( c[10] + minT * ( c[11] + minT * c[12] ))) )  // else if out of bounds - low
                             (                               c[1] + maxT * ( c[2] + maxT * ( c[ 3] + maxT * ( c[ 4] + maxT * c[ 5] ))) ); // else out of bounds - high
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
SpeciesHeatCapacity_Cv<FieldT>::
SpeciesHeatCapacity_Cv( const Expr::Tag& tTag,
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
      c[3] -= gasConstant; // change coefficients from isobaric to isometric
      c[3] /= molecularWeight; // convert to mass basis
      break;
    case NASA_POLY:
      c[1] -= 1.0; // change coefficients from isobaric to isometric
      c[8] -= 1.0;
      for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic){
        *ic *= gasConstant / molecularWeight; // dimensionalize the coefficients to mass basis
      }
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
SpeciesHeatCapacity_Cv<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  FieldT& cv = this->value();

  const FieldT& temp = t_->field_ref();

  const std::vector<double>& c = specTherm_.coefficients;
  const double minT = specTherm_.minTemp;
  const double maxT = specTherm_.maxTemp;
  switch ( specTherm_.type ) {
  case CONST_POLY:
    cv <<= c[3];
    break;
  case NASA_POLY:
    /* polynomials are applicable in two temperature ranges - high and low
     * If the temperature is out of range, the value is set to the value at the min or max temp
     */
    cv <<= cond( temp <= c[0] && temp >= minT, c[8] + temp * ( c[9] + temp * ( c[10] + temp * ( c[11] + temp * c[12] ))) )  // if low temp
               ( temp >  c[0] && temp <= maxT, c[1] + temp * ( c[2] + temp * ( c[ 3] + temp * ( c[ 4] + temp * c[ 5] ))) )  // else if high temp
               ( temp < minT,                  c[8] + minT * ( c[9] + minT * ( c[10] + minT * ( c[11] + minT * c[12] ))) )  // else if out of bounds - low
               (                               c[1] + maxT * ( c[2] + maxT * ( c[ 3] + maxT * ( c[ 4] + maxT * c[ 5] ))) ); // else out of bounds - high
    break;
  default: {
    std::ostringstream msg;
    msg << __FILE__ << " : " << __LINE__ << "\n Error for spec n = " << n_ << exceptionMsg_.str();
    throw std::runtime_error( msg.str() );
    }
  }
}

//--------------------------------------------------------------------

} // namespace pokitt

#endif // HeatCapacity_Cv_Expr_h
