/*
 * The MIT License
 *
 * Copyright (c) 2015-2017 The University of Utah
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

#ifndef HeatCapacity_Cp_Expr_h
#define HeatCapacity_Cp_Expr_h

#include <expression/Expression.h>

#include <pokitt/CanteraObjects.h> //include cantera wrapper

namespace pokitt{

/**
 *  \class  HeatCapacity_Cp
 *  \author Nate Yonkee
 *  \date July, 2014
 *
 *  \brief Calculates the constant pressure heat capacity of a mixture
 *  using NASA7 polynomials with 2 temperature ranges.
 *  Units of J/kg/K
 *
 *  Five coefficients \f$(a_0,\dots,a_4)\f$ are used to represent
 *  \f$ c_p^0(T)\f$ of species i as a polynomial in \f$ T \f$
 *
 * NASA7:
 *
 * \f[
 * \frac{c{p,i}(T)}{R} = a_0 + a_1 T + a_2 T^2 + a_3 T^3 + a_4 T^4
 * \f]
 *
 * The heat capacity is a weighted average of the pure species \f$ c_p^0(T)\f$,
 *
 * \f[
 * c_p(T) = \sum_{n=0}^{nSpec} Y_n c{p,n}
 * \f]
 */

template< typename FieldT >
class HeatCapacity_Cp
 : public Expr::Expression<FieldT>
{
  typedef std::vector<double> PolyVals; // values used for polynomial

  DECLARE_FIELD( FieldT, t_ )
  DECLARE_VECTOR_OF_FIELDS( FieldT, massFracs_ )

  const int nSpec_; // number of species
  std::vector< ThermData > specThermVec_;

  std::ostringstream exceptionMsg_; // generic exception to be thrown

  HeatCapacity_Cp( const Expr::Tag& tTag,
                   const Expr::TagList& massFracTags );

public:
  class Builder : public Expr::ExpressionBuilder
  {
    const Expr::Tag tTag_;
    const Expr::TagList massFracTags_;
  public:
    /**
     *  @brief Build a HeatCapacity_Cp expression
     *  @param resultTag the mixture heat capacity
     *  @param tTag tag for temperature
     *  @param massFracTags species mass fraction, indexing is consistent with Cantera input
     *  @param nghost number of ghost cells
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
      return new HeatCapacity_Cp<FieldT>( tTag_, massFracTags_ );
    }
  };

  ~HeatCapacity_Cp();
  void evaluate();
  void sensitivity( const Expr::Tag& var );
};

/**
 *  \class  SpeciesHeatCapacity_Cp
 *  \author Nate Yonkee
 *  \date May, 2014
 *
 *  \brief Calculates the constant pressure heat capacity of a single species
 *  using NASA7 polynomials with 2 temperature ranges.
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
 * Where \f[ t = T (K)/1000 \f]
 *
 */

template< typename FieldT >
class SpeciesHeatCapacity_Cp
 : public Expr::Expression<FieldT>
{
  DECLARE_FIELD( FieldT, t_ )

  const int n_; //index of species to be evaluated
  ThermData specTherm_;

  std::ostringstream exceptionMsg_; // generic exception to be thrown

  SpeciesHeatCapacity_Cp( const Expr::Tag& tTag,
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
     *  @param nghost number of ghost cells
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
      return new SpeciesHeatCapacity_Cp<FieldT>( tTag_, n_ );
    }
  };

  ~SpeciesHeatCapacity_Cp();
  void evaluate();
  void sensitivity( const Expr::Tag& var );
};


// ###################################################################
//
//                          Implementation
//
// ###################################################################



template< typename FieldT >
HeatCapacity_Cp<FieldT>::
HeatCapacity_Cp( const Expr::Tag& tTag,
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
      c[3] /= molecularWeights[n]; // convert to mass basis
      break;
    case NASA_POLY:
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
HeatCapacity_Cp<FieldT>::
~HeatCapacity_Cp()
{}

//--------------------------------------------------------------------

template< typename FieldT >
void
HeatCapacity_Cp<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  using namespace Cantera;

  FieldT& cp = this->value();

  const FieldT& t = t_->field_ref();

# ifndef ENABLE_CUDA
  const double maxTval = field_max_interior(t);
  const double minTval = field_min_interior(t);
# endif

  cp <<= 0.0; // set cp to 0 before starting summation

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
        cp <<= cp + yi * c[3];
        break;
        /* polynomials are applicable in two temperature ranges - high and low
         * If the temperature is out of range, the value is set to the value at the min or max temp
         */
      case NASA_POLY:
        cp <<= cp + yi * cond( t <= c[0] , c[8] + t  * ( c[9] + t  * ( c[10] + t  * ( c[11] + t  * c[12] ))) )  // if low temp
                             (             c[1] + t  * ( c[2] + t  * ( c[ 3] + t  * ( c[ 4] + t  * c[ 5] ))) );  // else if high temp
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
        cp <<= cp + yi * c[3];
        break;
        /* polynomials are applicable in two temperature ranges - high and low
         * If the temperature is out of range, the value is set to the value at the min or max temp
         */
      case NASA_POLY:
        cp <<= cp + yi * cond( t <= c[0] && t >= minT, c[8] + t    * ( c[9] + t    * ( c[10] + t    * ( c[11] + t    * c[12] ))) )  // if low temp
                             ( t >  c[0] && t <= maxT, c[1] + t    * ( c[2] + t    * ( c[ 3] + t    * ( c[ 4] + t    * c[ 5] ))) )  // else if high temp
                             ( t < minT,               c[8] + minT * ( c[9] + minT * ( c[10] + minT * ( c[11] + minT * c[12] ))) )  // else if out of bounds - low
                             (                         c[1] + maxT * ( c[2] + maxT * ( c[ 3] + maxT * ( c[ 4] + maxT * c[ 5] ))) ); // else out of bounds - high
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
void
HeatCapacity_Cp<FieldT>::
sensitivity( const Expr::Tag& var )
{
  using namespace SpatialOps;
  using namespace Cantera;

  FieldT& dcpdv = this->sensitivity_result( var );

  const FieldT& t = t_->field_ref();
  const FieldT& dTdv = t_->sens_field_ref( var );

# ifndef ENABLE_CUDA
  const double maxTval = field_max_interior(t);
  const double minTval = field_min_interior(t);
# endif

  dcpdv <<= 0.0; // set cp to 0 before starting summation

  for( size_t n=0; n<nSpec_; ++n ){
    const FieldT& yi = massFracs_[n]->field_ref();
    const FieldT& dyidv = massFracs_[n]->sens_field_ref( var );
    SpatFldPtr<FieldT> cpiPtr = SpatialFieldStore::get<FieldT>( yi );
    FieldT& cpi = *cpiPtr;
    const ThermData& thermo = specThermVec_[n];
    const ThermoPoly polyType = thermo.type;
    const std::vector<double>& c = thermo.coefficients;
    const double minT = thermo.minTemp;
    const double maxT = thermo.maxTemp;
#   ifndef ENABLE_CUDA // optimization benefits only the CPU - cond performs betters with if/else than with if/elif/elif/else
    if( maxTval <= maxT && minTval >= minT){ // if true, temperature can only be either high or low
      switch (polyType) {
      case CONST_POLY:
        dcpdv <<= dcpdv + dyidv * c[3];
        break;
        /* polynomials are applicable in two temperature ranges - high and low
         * If the temperature is out of range, the value is set to the value at the min or max temp
         */
      case NASA_POLY:
        cpi <<= cond( t <= c[0] , c[8] + t  * ( c[9] + t  * ( c[10] + t  * ( c[11] + t  * c[12] ))) )  // if low temp
                    (             c[1] + t  * ( c[2] + t  * ( c[ 3] + t  * ( c[ 4] + t  * c[ 5] ))) );  // else if high temp
        dcpdv <<= dcpdv + dyidv * cpi + yi * cond( t <= c[0] , dTdv  * ( c[9] + t  * ( 2.0 * c[10] + t  * ( 3.0 * c[11] + t  * 4.0 * c[12] ))) )  // if low temp
                                                 (             dTdv  * ( c[2] + t  * ( 2.0 * c[ 3] + t  * ( 3.0 * c[ 4] + t  * 4.0 * c[ 5] ))) );  // else if high temp
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
        dcpdv <<= dcpdv + dyidv * c[3];
        break;
        /* polynomials are applicable in two temperature ranges - high and low
         * If the temperature is out of range, the value is set to the value at the min or max temp
         */
      case NASA_POLY:
        dcpdv <<= dcpdv + dyidv * cpi + yi * cond( t <= c[0] && t >= minT, dTdv  * ( c[9] + t  * ( 2.0 * c[10] + t  * ( 3.0 * c[11] + t  * 4.0 * c[12] ))) )  // if low temp
                                                 ( t >  c[0] && t <= maxT, dTdv  * ( c[2] + t  * ( 2.0 * c[ 3] + t  * ( 3.0 * c[ 4] + t  * 4.0 * c[ 5] ))) )  // else if high temp
                                                 ( t < minT,               0.0 )  // else if out of bounds - low
                                                 (                         0.0 ); // else out of bounds - high
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
SpeciesHeatCapacity_Cp<FieldT>::
SpeciesHeatCapacity_Cp( const Expr::Tag& tTag,
                        const int n )
 : Expr::Expression<FieldT>(),
   n_ ( n ),
   specTherm_( CanteraObjects::species_thermo( n ) )
{
  this->set_gpu_runnable( true );

  exceptionMsg_ << "\nUnidentified polynomial type somehow evaded detection\n"
                << "This should be have been caught in Cantera Objects\n";

  t_ = this->template create_field_request<FieldT>(tTag);

  const std::vector<double>& molecularWeights = CanteraObjects::molecular_weights();
  const double molecularWeight = molecularWeights[n];
  const double gasConstant = CanteraObjects::gas_constant();
  std::vector<double>& c = specTherm_.coefficients;
  switch ( specTherm_.type ) {
  case CONST_POLY:
    c[3] /= molecularWeight; // convert to mass basis
    break;
  case NASA_POLY:
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
SpeciesHeatCapacity_Cp<FieldT>::
~SpeciesHeatCapacity_Cp()
{}

//--------------------------------------------------------------------

template< typename FieldT >
void
SpeciesHeatCapacity_Cp<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  FieldT& cp = this->value();
  const FieldT t = t_->field_ref();

  const std::vector<double>& c = specTherm_.coefficients;
  const double minT = specTherm_.minTemp;
  const double maxT = specTherm_.maxTemp;
  switch ( specTherm_.type ) {
  case CONST_POLY:
    cp <<= c[3];
    break;
  case NASA_POLY:
    /* polynomials are applicable in two temperature ranges - high and low
     * If the temperature is out of range, the value is set to the value at the min or max temp
     */
    cp <<= cond( t <= c[0] && t >= minT, c[8] + t    * ( c[9] + t    * ( c[10] + t    * ( c[11] + t    * c[12] ))) )  // if low temp
               ( t >  c[0] && t <= maxT, c[1] + t    * ( c[2] + t    * ( c[ 3] + t    * ( c[ 4] + t    * c[ 5] ))) )  // else if high temp
               ( t < minT,               c[8] + minT * ( c[9] + minT * ( c[10] + minT * ( c[11] + minT * c[12] ))) )  // else if out of bounds - low
               (                         c[1] + maxT * ( c[2] + maxT * ( c[ 3] + maxT * ( c[ 4] + maxT * c[ 5] ))) ); // else out of bounds - high
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
SpeciesHeatCapacity_Cp<FieldT>::
sensitivity( const Expr::Tag& var )
{
  using namespace SpatialOps;
  FieldT& dcpdv = this->sensitivity_result( var );
  const FieldT t = t_->field_ref();
  const FieldT dTdv = t_->sens_field_ref( var );

  const std::vector<double>& c = specTherm_.coefficients;
  const double minT = specTherm_.minTemp;
  const double maxT = specTherm_.maxTemp;
  switch ( specTherm_.type ) {
  case CONST_POLY:
    dcpdv <<= 0.0;
    break;
  case NASA_POLY:
    /* polynomials are applicable in two temperature ranges - high and low
     * If the temperature is out of range, the value is set to the value at the min or max temp
     */
    dcpdv <<= cond( t <= c[0] && t >= minT, dTdv * ( c[9] + t    * ( 2.0 * c[10] + t    * ( 3.0 * c[11] + t    * 4.0 * c[12] ))) )  // if low temp
                  ( t >  c[0] && t <= maxT, dTdv * ( c[2] + t    * ( 2.0 * c[ 3] + t    * ( 3.0 * c[ 4] + t    * 4.0 * c[ 5] ))) )  // else if high temp
                  ( t < minT,               0.0 )  // else if out of bounds - low
                  (                         0.0 ); // else out of bounds - high
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

#endif // HeatCapacity_Cp_Expr_h
