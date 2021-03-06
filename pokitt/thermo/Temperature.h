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

#ifndef Temperature_Expr_h
#define Temperature_Expr_h

#include <expression/Expression.h>

#include <pokitt/CanteraObjects.h> // include Cantera wrapper

#include <cmath> // isnan error checking

namespace pokitt{

/**
 *  \class  Temperature
 *  \author Nate Yonkee
 *  \date July, 2014
 *
 *  \brief Calculates the temperature of a mixture from
 *  enthalpy and mass fractions.
 *
 *  This class uses Newton's method to solve for the temperature of
 *  a mixture. The enthalpy of a mixture is a nonlinear
 *  function of temperature. Typically, this is represented by
 *  a six coefficient polynomial e.g. NASA7,
 *
 * \f[
 * \frac{h(T)}{R} = a_0 T + a_1 T^2/2 + a_2 T^3/3 + a_3 T^4/4 + a_4 T^5/5 + a_5
 * \f]
 *
 * Newton's method is used to solve this non-linear equation. For a
 * given \f$ h_0 \f$ and an initial guess \f$ T_i \f$ Newton's method calculates
 * a new guess,
 *
 * \f[
 * T_{i+1} = T_{i} + \frac{h_0-h(T_i)}{c_p(T_i)}
 * \f]
 *
 * This iteration continues until a convergence criteria is met:
 *
 * \f[
 * \mid T_{i+1} - T_{i} \mid < 10^{-3}
 * \f]
 *
 * Sensitivity:
 * \f[
 * dT = dh/cp - \sum_{i=1}^{i=nspec-1} (h_i-h_{nspec})/cp * dY_i
 * \f]
 * Here, \f$nspec \f$ is the number of species.
 * In our present codes, e.g., Zodiac, ODT, we choose \f$\rho \f$, \f$T \f$, and \f$Y_i \f$ as the primitive variables. As a result, we want \f$dT/dY_i=0 \f$ and \f$dT/d\rho=0 \f$,
 * and we do not want to use chain rule based on the tree to get \f$dT/dY_i = - \sum_{i=1}^{i=nspec-1} (h_i-h_{nspec})/cp \f$.
 * Therefore, override_sensitivity is used here. If the sensitivity variable is \f$T \f$, the sensitivity is set to be one \f$dT/dT = 1 \f$.
 * If the sensitivity variable is \f$h \f$, the sensitivity is \f$dT/dh = 1/cp \f$.
 * Otherwise, the sensitivity is set to be zero. \f$dT/d\rho = 0 \f$. \f$dT/dY_i = 0 \f$.
 * For some other sensitivities, like \f$dT/d(\rho h) \f$, it needs to be given directly by hand, not by this class.
 */

template< typename FieldT >
class Temperature
    : public Expr::Expression<FieldT>
{
  typedef std::vector<FieldT*> SpecT;
  const Expr::Tag enthTag_;
  const Expr::TagList massFracTags_;

  DECLARE_FIELDS( FieldT, enth_, tguess_ )
  DECLARE_VECTOR_OF_FIELDS( FieldT, massFracs_ )

  std::vector< ThermData > specThermVec_;
  std::vector< std::vector<double> > cFracVec_; // vector of polynomial coefficients divided by integers carried from integration
  const size_t nSpec_; // number of species to iterate over

  const double tol_; // tolerance for Newton's method
  const double maxTemp_; //temperatures above this will throw an exception, default is 5000K
  const int maxIterations_; // number of iterations before exception is thrown, default is 20

  std::ostringstream exceptionMsg_; // generic exception to be thrown

  Temperature( const Expr::TagList& massFracTags,
               const Expr::Tag& enthTag,
               const Expr::Tag& temperatureGuessTag,
               const double tol,
               const double maxTemp,
               const int maxIterations);
public:

  class Builder : public Expr::ExpressionBuilder
  {
    const Expr::Tag enthTag_, tempGuessTag_;
    const Expr::TagList massFracTags_;
    const double tol_, maxTemp_;
    const int maxIterations_;
  public:
    /**
     *  @brief Build a Temperature expression
     *  @param resultTag the tag for the temperature of the mixture
     *  @param massFracTags tag for the mass fraction of each species, ordering is consistent with Cantera input file
     *  @param enthTag tag for the enthalpy of the mixture
     *  @param temperatureGuessTag if supplied, this field will be used as the
     *         guess for temperature.  If empty, we assume that the result field
     *         is pre-populated with the guess value.
     *  @param tol     tolerance for the non-linear solve
     *  @param maxTemp the maximum allowable temperature
     *  @param maxIterations the maximum iteration count
     *  @param nghost number of ghost cells
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::TagList& massFracTags,
             const Expr::Tag& enthTag,
             const Expr::Tag& temperatureGuessTag,
             const double tol = 1e-3,
             const double maxTemp = 5000,
             const int maxIterations = 20,
             const SpatialOps::GhostData nghost = DEFAULT_NUMBER_OF_GHOSTS )
    : ExpressionBuilder( resultTag, nghost ),
      massFracTags_( massFracTags ),
      enthTag_( enthTag ),
      tempGuessTag_( temperatureGuessTag ),
      tol_( tol ),
      maxTemp_( maxTemp ),
      maxIterations_( maxIterations )
    {
      if( resultTag == temperatureGuessTag ){
        throw std::runtime_error( "Error! Do not provide temperature's tag as the guess tag! Use an empty tag to use T as the guess." );
      }
    }

    Expr::ExpressionBase* build() const{
      return new Temperature<FieldT>( massFracTags_, enthTag_, tempGuessTag_, tol_, maxTemp_, maxIterations_ );
    }

  };

  ~Temperature(){}
  void evaluate();
  bool override_sensitivity() const{return true;}
  void sensitivity( const Expr::Tag& );
  void find_bad_points( std::ostringstream& msg, const FieldT& badField, const double badValue, const bool checkBelow );
};

/**
 *  \class  TemperatureFromE0
 *  \author Nate Yonkee
 *  \date July, 2014
 *
 *  \brief Calculates the temperature of a mixture from
 *  total energy (internal and kinetic energy) and mass fractions.
 *
 *  This class uses Newton's method to solve for the temperature of
 *  a mixture. The internal energy of a mixture is a nonlinear
 *  function of temperature. Typically, this is represented by
 *  a six coefficient polynomial e.g. NASA7,
 *
 * \f[
 * \frac{e(T)}{R} = (a_0 - 1) T + a_1 T^2/2 + a_2 T^3/3 + a_3 T^4/4 + a_4 T^5/5 + a_5
 * \f]
 *
 * Newton's method is used to solve this non-linear equation. For a
 * given \f$ e_0 \f$ and an initial guess \f$ T_i \f$ Newton's method calculates
 * a new guess,
 *
 * \f[
 * T_{i+1} = T_{i} + \frac{e_0-e(T_i)}{c_v(T_i)}
 * \f]
 *
 * This iteration continues until a convergence criteria is met:
 *
 * \f[
 * \mid T_{i+1} - T_{i} \mid < 10^{-4}
 * \f]
 *
 * Sensitivity:
 * \f[
 * dT = de/cv - \sum_{i=1}^{i=nspec-1} (e_i-e_{nspec})/cv * dY_i
 * \f]
 * Here, \f$nspec \f$ is the number of species.
 * In our present codes, e.g., Zodiac, ODT, we choose \f$\rho \f$, \f$T \f$, and \f$Y_i \f$ as the primitive variables. As a result, we want \f$dT/dY_i=0 \f$ and \f$dT/d\rho=0 \f$,
 * and we do not want to use chain rule based on the tree to get \f$dT/dY_i = - \sum_{i=1}^{i=nspec-1} (e_i-e_{nspec})/cv \f$.
 * Therefore, override_sensitivity is used here. If the sensitivity variable is \f$T \f$, the sensitivity is set to be one \f$dT/dT = 1 \f$.
 * If the sensitivity variable is \f$e \f$, the sensitivity is \f$dT/de = 1/cv \f$.
 * Otherwise, the sensitivity is set to be zero. \f$dT/d\rho = 0 \f$. \f$dT/dY_i = 0 \f$.
 * For some other sensitivities, like \f$dT/d(\rho e) \f$, it needs to be given directly by hand, not by this class.
 */

template< typename FieldT >
class TemperatureFromE0
    : public Expr::Expression<FieldT>
{
  typedef std::vector<FieldT*> SpecT;

  DECLARE_FIELDS( FieldT, e0_, ke_, tguess_ )
  DECLARE_VECTOR_OF_FIELDS( FieldT, massFracs_ )

  std::vector< ThermData > specThermVec_;
  std::vector< std::vector<double> > cFracVec_; // vector of polynomial coefficients divided by integers carried from integration
  const size_t nSpec_; // number of species to iterate over

  const double tol_; // tolerance for Newton's method
  const double maxTemp_; //temperatures above this will throw an exception, default is 5000K
  const int maxIterations_; // number of iterations before exception is thrown, default is 20

  std::ostringstream exceptionMsg_; // generic exception to be thrown

  TemperatureFromE0( const Expr::TagList& massFracTags,
                     const Expr::Tag& e0Tag,
                     const Expr::Tag& keTag,
                     const Expr::Tag& temperatureGuessTag,
                     const double tol,
                     const double maxTemp,
                     const int maxIterations );
public:

  class Builder : public Expr::ExpressionBuilder
  {
    const Expr::Tag e0Tag_;
    const Expr::TagList massFracTags_;
    const Expr::Tag keTag_, tempGuessTag_;
    const double tol_, maxTemp_;
    const int maxIterations_;
  public:
    /**
     *  @brief Build a TemperatureFromE0 expression
     *  @param resultTag the tag for the temperature of the mixture
     *  @param massFracTags tag for the mass fraction of each species, ordering is consistent with Cantera input file
     *  @param e0Tag tag for the total energy of the mixture
     *  @param keTag tag for kinetic energy of the mixture
     *  @param temperatureGuessTag if supplied, this field will be used as the
     *         guess for temperature.  If empty, we assume that the result field
     *         is pre-populated with the guess value.
     *  @param tol   tolerance for the non-linear solve
     *  @param maxTemp the maximum allowable temperature
     *  @param maxIterations the maximum iteration count
     *  @param nghost number of ghost cells
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::TagList& massFracTags,
             const Expr::Tag& e0Tag,
             const Expr::Tag& keTag,
             const Expr::Tag& temperatureGuessTag,
             const double tol = 1e-3,
             const double maxTemp = 5000,
             const int maxIterations = 20,
             const SpatialOps::GhostData nghost = DEFAULT_NUMBER_OF_GHOSTS )
    : ExpressionBuilder( resultTag, nghost ),
      massFracTags_( massFracTags ),
      e0Tag_( e0Tag ),
      keTag_( keTag ),
      tempGuessTag_( temperatureGuessTag ),
      tol_( tol ),
      maxTemp_( maxTemp ),
      maxIterations_( maxIterations )
    {
      if( resultTag == temperatureGuessTag ){
        throw std::runtime_error( "Error! Do not provide temperature's tag as the guess tag! Use an empty tag to use T as the guess." );
      }
    }

    Expr::ExpressionBase* build() const{
      return new TemperatureFromE0<FieldT>( massFracTags_, e0Tag_, keTag_, tempGuessTag_, tol_, maxTemp_, maxIterations_ );
    }
  };

  ~TemperatureFromE0(){}
  void evaluate();
  bool override_sensitivity() const{return true;}
  void sensitivity( const Expr::Tag& );
  void find_bad_points( std::ostringstream& msg, const FieldT& badField, const double badValue, const bool checkBelow );
};

// ###################################################################
//
//                          Implementation
//
// ###################################################################

template< typename FieldT >
Temperature<FieldT>::
Temperature( const Expr::TagList& massFracTags,
             const Expr::Tag& enthTag,
             const Expr::Tag& temperatureGuessTag,
             const double tol,
             const double maxTemp,
             const int maxIterations )
  : Expr::Expression<FieldT>(),
    tol_( tol ),
    maxTemp_( maxTemp ),
    maxIterations_( maxIterations ),
    nSpec_( CanteraObjects::number_species() )
{
  this->set_gpu_runnable( true );

  exceptionMsg_ << "\nUnidentified polynomial type somehow evaded detection\n"
                << "This should be have been caught in Cantera Objects\n";

  enth_ = this->template create_field_request<FieldT>( enthTag );
  this->template create_field_vector_request<FieldT>( massFracTags, massFracs_ );

  if( temperatureGuessTag != Expr::Tag() ){
    tguess_ = this->template create_field_request<FieldT>( temperatureGuessTag );
  }

  const std::vector<double>& molecularWeights = CanteraObjects::molecular_weights();
  const double gasConstant = CanteraObjects::gas_constant();

  for( size_t n=0; n<nSpec_; ++n ){
    ThermData tData = CanteraObjects::species_thermo( n );
    std::vector<double>& c = tData.coefficients;
    ThermoPoly type = tData.type;
    std::vector<double> cFrac = c;
    switch ( type ) {
    case CONST_POLY:
      c[1] /= molecularWeights[n]; // convert to mass basis
      c[3] /= molecularWeights[n]; // convert to mass basis
      break;
    case NASA_POLY:
      for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic){
        *ic *= gasConstant / molecularWeights[n]; // dimensionalize the coefficients to mass basis
      }
      cFrac[ 2] = c[ 2] / 2.0;
      cFrac[ 3] = c[ 3] / 3.0;
      cFrac[ 4] = c[ 4] / 4.0;
      cFrac[ 5] = c[ 5] / 5.0;
      cFrac[ 9] = c[ 9] / 2.0;
      cFrac[10] = c[10] / 3.0;
      cFrac[11] = c[11] / 4.0;
      cFrac[12] = c[12] / 5.0;
      break;
    default: {
      std::ostringstream msg;
      msg << __FILE__ << " : " << __LINE__ << "\n Error for spec n = " << n << exceptionMsg_.str();
      throw std::runtime_error( msg.str() );
      }
    }
    cFracVec_.push_back(cFrac);
    specThermVec_.push_back( tData );
  }
}

//--------------------------------------------------------------------

template< typename FieldT >
void
Temperature<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  FieldT& temp = this->value();

  if( tguess_ ){
    temp <<= tguess_->field_ref();
  }

  const FieldT& enth = enth_->field_ref();

  SpatFldPtr<FieldT> delHPtr = SpatialFieldStore::get<FieldT>(temp); // difference between enthalpy field value and enthalpy evaluated at current temperature
  SpatFldPtr<FieldT> dhdTPtr = SpatialFieldStore::get<FieldT>(temp); // dhdT for Newton's method
  SpatFldPtr<FieldT> resPtr  = SpatialFieldStore::get<FieldT>(temp); // change in temperature for new Newton's iteration

  FieldT& delH = *delHPtr;
  FieldT& dhdT = *dhdTPtr;
  FieldT& res  = *resPtr;
# ifdef ENABLE_CUDA
  res.add_device(CPU_INDEX);
# endif
  bool isConverged = false;
  int iterations = 0;
  while( !isConverged ){
    delH <<= enth;
    dhdT <<= 0.0;
#   ifndef ENABLE_CUDA
    const double maxTval = field_max_interior(temp);
    const double minTval = field_min_interior(temp);
    if( maxTval >= maxTemp_ ){
      std::ostringstream msg;
        msg << std::endl
            << "Error in pokitt::Temperature::evaluate()." << std::endl
            << "Temperature is too high" << std::endl;
        find_bad_points( msg, temp, maxTemp_, false );
        msg << "Total   iterations = " << iterations << std::endl
            << "Maximum iterations = " << maxIterations_ << std::endl
            << "Set tolerance      = " << tol_ << std::endl
            << __FILE__ << " : " << __LINE__ << std::endl;
        throw std::runtime_error( msg.str() );
    }
    if( minTval <= 0.0 ){
      std::ostringstream msg;
        msg << std::endl
            << "Error in pokitt::Temperature::evaluate()." << std::endl
            << "Temperature is below 0" << std::endl;
        find_bad_points( msg, temp, 0, true );
        msg << "Total   iterations = " << iterations << std::endl
            << "Maximum iterations = " << maxIterations_ << std::endl
            << "Set tolerance      = " << tol_ << std::endl
            << __FILE__ << " : " << __LINE__ << std::endl;
        throw std::runtime_error( msg.str() );
    }
    if( std::isnan( nebo_sum_interior( temp ) ) ){
      std::ostringstream msg;
        msg << std::endl
            << "Error in pokitt::Temperature::evaluate()." << std::endl
            << "Temperature is NaN" << std::endl;
        find_bad_points( msg, temp, NAN, false );
        msg << "Total   iterations = " << iterations << std::endl
            << "Maximum iterations = " << maxIterations_ << std::endl
            << "Set tolerance      = " << tol_ << std::endl
            << __FILE__ << " : " << __LINE__ << std::endl;
        throw std::runtime_error( msg.str() );
    }
#   endif
    for( size_t n=0; n<nSpec_; ++n ){
      const FieldT& yi = massFracs_[n]->field_ref();
      const ThermData& specTherm = specThermVec_[n];
      const ThermoPoly type = specTherm.type;
      const std::vector<double>& c = specTherm.coefficients;
      const std::vector<double>& cFrac = cFracVec_[n];
      const double minT = specTherm.minTemp;
      const double maxT = specTherm.maxTemp;
#     ifndef ENABLE_CUDA // optimization benefits only the CPU - cond performs betters with if/else than with if/elif/elif/else
      if( maxTval <= maxT && minTval >= minT){ // if true, temperature can only be either high or low
        switch ( type ) {
        case CONST_POLY: // constant cp
          delH <<= delH - yi * ( c[1] + c[3]*(temp-c[0]) );
          dhdT <<= dhdT + yi * c[3];
          break;
        case NASA_POLY:
          delH <<= delH - yi
                 * cond( temp <= c[0], c[13] + temp * ( c[8] + temp * ( cFrac[9] + temp * ( cFrac[10] + temp * ( cFrac[11] + temp * cFrac[12] )))) )  // if low temp
                       (               c[ 6] + temp * ( c[1] + temp * ( cFrac[2] + temp * ( cFrac[ 3] + temp * ( cFrac[ 4] + temp * cFrac[ 5] )))) );  // else if high temp

          dhdT <<= dhdT + yi
                 * cond( temp <= c[0], c[8] + temp * ( c[9] + temp * ( c[10] + temp * ( c[11] + temp * c[12] ))) )  // if low temp
                       (               c[1] + temp * ( c[2] + temp * ( c[ 3] + temp * ( c[ 4] + temp * c[ 5] ))) );  // else if high temp
          break;
        default: {
          std::ostringstream msg;
          msg << __FILE__ << " : " << __LINE__ << "\n Error for spec n = " << n << exceptionMsg_.str();
          throw std::runtime_error( msg.str() );
          }
        } // switch( type )
      }
      else
#     endif
      {
        /* else temperature can be out of bounds low, low temp, high temp, or out of bounds high
         * if out of bounds, properties are interpolated from min or max temp using a constant cp
         */
        switch ( type ){
        case CONST_POLY: // constant cp
          delH<<=delH - yi * ( c[1] + c[3]*(temp-c[0]) );
          dhdT<<=dhdT + yi * c[3];
          break;
        case NASA_POLY:
          delH <<= delH - yi
                 * cond( temp <= c[0] && temp >= minT, c[13] + temp * ( c[8] + temp * ( cFrac[9] + temp * ( cFrac[10] + temp * ( cFrac[11] + temp * cFrac[12] )))) )  // if low temp
                       ( temp >  c[0] && temp <= maxT, c[ 6] + temp * ( c[1] + temp * ( cFrac[2] + temp * ( cFrac[ 3] + temp * ( cFrac[ 4] + temp * cFrac[ 5] )))) )  // else if high temp
                       ( temp < minT,                  c[13] + c[8] * temp + minT * ( c[9] * temp + minT * ( c[10] * temp - cFrac[9] + minT * ( c[11] * temp - 2*cFrac[10] + minT * ( c[12] * temp - 3*cFrac[11] + minT * -4*cFrac[12] )))) )  // else if out of bounds - low
                       (                               c[ 6] + c[1] * temp + maxT * ( c[2] * temp + maxT * ( c[ 3] * temp - cFrac[2] + maxT * ( c[ 4] * temp - 2*cFrac[ 3] + maxT * ( c[ 5] * temp - 3*cFrac[ 4] + maxT * -4*cFrac[ 5] )))) ); // else out of bounds - high

          dhdT <<= dhdT + yi
                 * cond( temp <= c[0] && temp >= minT, c[8] + temp * ( c[9] + temp * ( c[10] + temp * ( c[11] + temp * c[12] ))) )  // if low temp
                       ( temp >  c[0] && temp <= maxT, c[1] + temp * ( c[2] + temp * ( c[ 3] + temp * ( c[ 4] + temp * c[ 5] ))) )  // else if high temp
                       ( temp < minT,                  c[8] + minT * ( c[9] + minT * ( c[10] + minT * ( c[11] + minT * c[12] ))) )  // else if out of bounds - low
                       (                               c[1] + maxT * ( c[2] + maxT * ( c[ 3] + maxT * ( c[ 4] + maxT * c[ 5] ))) ); // else out of bounds - high

          break;
        default: {
          std::ostringstream msg;
          msg << __FILE__ << " : " << __LINE__ << "\n Error for spec n = " << n << exceptionMsg_.str();
          throw std::runtime_error( msg.str() );
          }
        } // switch( type )
      }
    } // species loop

    // Newton's method to find root
    res  <<= delH/dhdT;
    temp <<= temp + res;
#   ifdef ENABLE_CUDA
    res.set_device_as_active(CPU_INDEX);
#   endif
    const double err = field_max_interior( abs(res) );
#   ifdef ENABLE_CUDA
    res.set_device_as_active(GPU_INDEX);
#   endif
    isConverged = err < tol_; // Converged when the temperature has changed by less than specified tolerance
    ++iterations;

    if( !isConverged && iterations == maxIterations_ ){
      std::ostringstream msg;
      msg << std::endl
          << "Error in pokitt::Temperature::evaluate()." << std::endl
          << "Iteration count exceeded" << std::endl
          << "Total   iterations = " << iterations << std::endl
          << "Maximum iterations = " << maxIterations_ << std::endl
          << "Set tolerance              = " << tol_ << std::endl
          << "Largest pointwise residual = " << err << std::endl
          << __FILE__ << " : " << __LINE__ << std::endl;
      throw std::runtime_error( msg.str() );
    }
  }
}
    
//--------------------------------------------------------------------
template< typename FieldT >
void
Temperature<FieldT>::
sensitivity( const Expr::Tag& var)
{
  using namespace SpatialOps;
  FieldT& dfdv = this->sensitivity_result( var );
  if( var == enth_->tag() ){
    FieldT& temp = this->value();
    SpatFldPtr <FieldT> cpPtr = SpatialFieldStore::get<FieldT>( temp ); //heat capacity of mixture
    FieldT& cp = *cpPtr;
    cp <<= 0.0;
    const double maxTval = field_max_interior( temp );
    const double minTval = field_min_interior( temp );
    for( size_t n = 0;n < nSpec_;++n ){
      const FieldT& yi = massFracs_[n]->field_ref();
      const ThermData& specTherm = specThermVec_[n];
      const ThermoPoly type = specTherm.type;
      const std::vector<double>& c = specTherm.coefficients;
      const double minT = specTherm.minTemp;
      const double maxT = specTherm.maxTemp;
#     ifndef ENABLE_CUDA // optimization benefits only the CPU - cond performs betters with if/else than with if/elif/elif/else
      if( maxTval <= maxT && minTval >= minT ){ // if true, temperature can only be either high or low
        switch( type ){
          case CONST_POLY: // constant cp
            cp <<= cp + yi * c[3];
            break;
          case NASA_POLY:
            cp <<= cp + yi
                        * cond( temp <= c[0], c[8] + temp * ( c[9] + temp * ( c[10] + temp * ( c[11] + temp * c[12] ))))   // if low temp
                              ( c[1] + temp * ( c[2] + temp * ( c[3] + temp * ( c[4] + temp * c[5] ))));  // else if high temp
            break;
          default:{
            std::ostringstream msg;
            msg << __FILE__ << " : " << __LINE__ << "\n Error for spec n = " << n << exceptionMsg_.str();
            throw std::runtime_error( msg.str());
          }
        } // switch( type )
      }
      else
#     endif
      {
        /* else temperature can be out of bounds low, low temp, high temp, or out of bounds high
         * if out of bounds, properties are keeped to the same as that for min or max temp
         */
        switch( type ){
          case CONST_POLY: // constant cp
            cp <<= cp + yi * c[3];
            break;
          case NASA_POLY:
            cp <<= cp + yi
                        * cond( temp <= c[0] && temp >= minT, c[8] + temp * ( c[9] + temp * ( c[10] + temp * ( c[11] + temp * c[12] ))))  // if low temp
                              ( temp > c[0] && temp <= maxT, c[1] + temp * ( c[2] + temp * ( c[3] + temp * ( c[4] + temp * c[5] ))))  // else if high temp
                              ( temp < minT, c[8] + minT * ( c[9] + minT * ( c[10] + minT * ( c[11] + minT * c[12] ))))  // else if out of bounds - low
                              ( c[1] + maxT * ( c[2] + maxT * ( c[3] + maxT * ( c[4] + maxT * c[5] )))); // else out of bounds - high

            break;
          default:{
            std::ostringstream msg;
            msg << __FILE__ << " : " << __LINE__ << "\n Error for spec n = " << n << exceptionMsg_.str();
            throw std::runtime_error( msg.str());
          }
        } // switch( type )
      }
    } // species loop
    dfdv <<= 1 / cp;
  }
  else{
    dfdv <<= 0.0;
  }
}

//--------------------------------------------------------------------

template< typename FieldT >
void
Temperature<FieldT>::
find_bad_points( std::ostringstream& msg, const FieldT& badField, const double badValue, const bool checkBelow )
{
  bool checkNaN = false;
  if( std::isnan( badValue ) ) checkNaN = true;
  typename FieldT::const_iterator iField = badField.begin();
  const SpatialOps::MemoryWindow mw = badField.window_with_ghost();
  int badPoints = 0;
  double worstValue = 0;
  int xWorst, yWorst, zWorst;

  for( int z = 1; z <= mw.extent(2); ++z ){
    for( int y = 1; y <= mw.extent(1); ++y ){
      for( int x = 1; x <= mw.extent(0); ++x, ++iField ){
        if( checkNaN ){
          if( std::isnan( *iField ) ){
            ++badPoints;
            if( badPoints < 6 )
              msg << "Value = " << *iField << " at [x,y,z] = [" << x << "," <<  y << "," << z << "]" << std::endl;
            worstValue = *iField;
            xWorst = x; yWorst = y; zWorst = z;
          }
        }
        else if( checkBelow ){
          if( *iField <= badValue  ){
            ++badPoints;
            if( badPoints < 6 )
              msg << "Value = " << *iField << " at [x,y,z] = [" << x << "," <<  y << "," << z << "]" << std::endl;
          }
          if( *iField < worstValue ){
            worstValue = *iField;
            xWorst = x; yWorst = y; zWorst = z;
          }

        }
        else{
          if( *iField >= badValue ){
            ++badPoints;
            if( badPoints < 6 )
              msg << "Value = " << *iField << " at [x,y,z] = [" << x << "," <<  y << "," << z << "]" << std::endl;
          }
          if( *iField > worstValue ){
            worstValue = *iField;
            xWorst = x; yWorst = y; zWorst = z;
          }
        }
      }
    }
  }
  msg << std::endl
      << "Worst value = " << worstValue << " at [x,y,z] = [" << xWorst << "," <<  yWorst << "," << zWorst << "]" << std::endl
      << "A total of " << badPoints << " bad point(s) detected." << std::endl;
}

//--------------------------------------------------------------------

template< typename FieldT >
TemperatureFromE0<FieldT>::
TemperatureFromE0( const Expr::TagList& massFracTags,
                   const Expr::Tag& e0Tag,
                   const Expr::Tag& keTag,
                   const Expr::Tag& temperatureGuessTag,
                   const double tol,
                   const double maxTemp,
                   const int maxIterations )
  : Expr::Expression<FieldT>(),
    tol_( tol ),
    maxTemp_( maxTemp ),
    maxIterations_( maxIterations ),
    nSpec_( CanteraObjects::number_species() )
{
  this->set_gpu_runnable( true );

  exceptionMsg_ << "\nUnidentified polynomial type somehow evaded detection\n"
                << "This should be have been caught in Cantera Objects\n";

  e0_ = this->template create_field_request<FieldT>( e0Tag );
  ke_ = this->template create_field_request<FieldT>( keTag );
  this->template create_field_vector_request<FieldT>( massFracTags, massFracs_ );

  if( temperatureGuessTag != Expr::Tag() ){
    tguess_ = this->template create_field_request<FieldT>( temperatureGuessTag );
  }

  const std::vector<double>& molecularWeights = CanteraObjects::molecular_weights();
  const double gasConstant = CanteraObjects::gas_constant();

  for( size_t n=0; n<nSpec_; ++n ){
    ThermData tData = CanteraObjects::species_thermo( n );
    std::vector<double>& c = tData.coefficients;
    ThermoPoly type = tData.type;
    std::vector<double> cFrac = c;
    switch ( type ) {
    case CONST_POLY:
      c[3] -= gasConstant; // change coefficients from isobaric to isometric
      c[1] -= gasConstant * c[0];
      c[3] /= molecularWeights[n]; // convert to mass basis
      c[1] /= molecularWeights[n];
      break;
    case NASA_POLY:
      c[1] -= 1.0; //change coefficients from isobaric to isometric
      c[8] -= 1.0;
      for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic){
        *ic *= gasConstant / molecularWeights[n]; // dimensionalize the coefficients to mass basis
      }
      cFrac[ 2] = c[ 2] / 2;
      cFrac[ 3] = c[ 3] / 3;
      cFrac[ 4] = c[ 4] / 4;
      cFrac[ 5] = c[ 5] / 5;
      cFrac[ 9] = c[ 9] / 2;
      cFrac[10] = c[10] / 3;
      cFrac[11] = c[11] / 4;
      cFrac[12] = c[12] / 5;
      break;
    default: {
      std::ostringstream msg;
      msg << __FILE__ << " : " << __LINE__ << "\n Error for spec n = " << n << exceptionMsg_.str();
      throw std::runtime_error( msg.str() );
      }
    }
    cFracVec_.push_back( cFrac );
    specThermVec_.push_back( tData );
  }

}

//--------------------------------------------------------------------

template< typename FieldT >
void
TemperatureFromE0<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  FieldT& temp = this->value();

  if( tguess_ ){
    temp <<= tguess_->field_ref();
  }
  const FieldT& e0 = e0_->field_ref();
  const FieldT& ke = ke_->field_ref();

  SpatFldPtr<FieldT> delE0Ptr = SpatialFieldStore::get<FieldT>(temp); // difference between internal energy field value and internal energy at current temperature
  SpatFldPtr<FieldT> dE0dTPtr = SpatialFieldStore::get<FieldT>(temp); // dE0dT for Newton's method
  SpatFldPtr<FieldT> resPtr   = SpatialFieldStore::get<FieldT>(temp); // change in temperature for new Newton's iteration

  FieldT& delE0 = *delE0Ptr;
  FieldT& dE0dT = *dE0dTPtr;
  FieldT& res   = *resPtr;
# ifdef ENABLE_CUDA
  res.add_device(CPU_INDEX);
# endif

  bool isConverged = false;
  int iterations = 0;
  delE0 <<= e0 - ke;
  while( !isConverged ){
    delE0 <<= e0 - ke;
    dE0dT <<= 0.0;
#   ifndef ENABLE_CUDA
    const double maxTval = field_max_interior(temp);
    const double minTval = field_min_interior(temp);
    if( maxTval >= maxTemp_ ){
      std::ostringstream msg;
      msg << std::endl
          << "Error in pokitt::TemperatureFromE0::evaluate()." << std::endl
          << "Temperature is too high" << std::endl;
      find_bad_points( msg, temp, maxTemp_, false );
      msg << "Total   iterations = " << iterations << std::endl
          << "Maximum iterations = " << maxIterations_ << std::endl
          << "Set tolerance      = " << tol_ << std::endl
          << __FILE__ << " : " << __LINE__ << std::endl;
      throw std::runtime_error( msg.str() );
    }
    if( minTval <= 0.0 ){
      std::ostringstream msg;
      msg << std::endl
          << "Error in pokitt::TemperatureFromE0::evaluate()." << std::endl
          << "Temperature is below 0" << std::endl;
      find_bad_points( msg, temp, 0, true );
      msg << "Total   iterations = " << iterations << std::endl
          << "Maximum iterations = " << maxIterations_ << std::endl
          << "Set tolerance      = " << tol_ << std::endl
          << __FILE__ << " : " << __LINE__ << std::endl;
      throw std::runtime_error( msg.str() );
    }
    if( std::isnan( nebo_sum( temp ) ) ){
      std::ostringstream msg;
      msg << std::endl
          << "Error in pokitt::TemperatureFromE0::evaluate()." << std::endl
          << "Temperature is NaN" << std::endl;
      find_bad_points( msg, temp, NAN, false );
      msg << "Total   iterations = " << iterations << std::endl
          << "Maximum iterations = " << maxIterations_ << std::endl
          << "Set tolerance      = " << tol_ << std::endl
          << __FILE__ << " : " << __LINE__ << std::endl;
      throw std::runtime_error( msg.str() );
    }
#   endif
    for( size_t n=0; n<nSpec_; ++n){
      const FieldT& yi = massFracs_[n]->field_ref();
      const ThermData& specTherm = specThermVec_[n];
      const ThermoPoly type = specTherm.type;
      const std::vector<double>& c = specTherm.coefficients;
      const std::vector<double>& cFrac = cFracVec_[n];
      const double minT = specTherm.minTemp;
      const double maxT = specTherm.maxTemp;
#     ifndef ENABLE_CUDA // optimization benefits only the CPU - cond performs betters with if/else than with if/elif/elif/else
      if( maxTval <= maxT && minTval >= minT){ // if true, temperature can only be either high or low
        switch ( type ) {
        case CONST_POLY: // constant cv
          delE0 <<= delE0 - yi * ( c[1] + c[3]*(temp-c[0]) );
          dE0dT <<= dE0dT + yi * c[3];
          break;
        case NASA_POLY:
          delE0 <<= delE0 - yi
                  * cond( temp <= c[0], c[13] + temp * ( c[8] + temp * ( cFrac[9] + temp * ( cFrac[10] + temp * ( cFrac[11] + temp * cFrac[12] )))) )  // if low temp
                        (               c[ 6] + temp * ( c[1] + temp * ( cFrac[2] + temp * ( cFrac[ 3] + temp * ( cFrac[ 4] + temp * cFrac[ 5] )))) ); // else if high temp

          dE0dT <<= dE0dT + yi
                  * cond( temp <= c[0], c[8] + temp * ( c[9] + temp * ( c[10] + temp * ( c[11] + temp * c[12] ))) )  // if low temp
                        (               c[1] + temp * ( c[2] + temp * ( c[ 3] + temp * ( c[ 4] + temp * c[ 5] ))) ); // else if high temp
          break;
        default: {
          std::ostringstream msg;
          msg << __FILE__ << " : " << __LINE__ << "\n Error for spec n = " << n << exceptionMsg_.str();
          throw std::runtime_error( msg.str() );
          }
        } // switch( polyType )
      }
      else
#     endif
      {
        /* else temperature can be out of bounds low, low temp, high temp, or out of bounds high
         * if out of bounds, properties are interpolated from min or max temp using a constant cv
         */
        switch ( type ){
        case CONST_POLY: // constant cv
          delE0 <<= delE0 - yi * ( c[1] + c[3] * (temp-c[0]) );
          dE0dT <<= dE0dT + yi * c[3];
          break;
        case NASA_POLY:
          delE0 <<= delE0 - yi
                  * cond( temp <= c[0] && temp >= minT, c[13] + temp * ( c[8] + temp * ( cFrac[9] + temp * ( cFrac[10] + temp * ( cFrac[11] + temp * cFrac[12] )))) )  // if low temp
                        ( temp >  c[0] && temp <= maxT, c[ 6] + temp * ( c[1] + temp * ( cFrac[2] + temp * ( cFrac[ 3] + temp * ( cFrac[ 4] + temp * cFrac[ 5] )))) )  // else if high temp
                        ( temp < minT,                  c[13] + c[8] * temp + minT * ( c[9] * temp + minT * ( c[10] * temp - cFrac[9] + minT * ( c[11] * temp - 2*cFrac[10] + minT * ( c[12] * temp - 3*cFrac[11] + minT * -4*cFrac[12] )))) )  // else if out of bounds - low
                        (                               c[ 6] + c[1] * temp + maxT * ( c[2] * temp + maxT * ( c[ 3] * temp - cFrac[2] + maxT * ( c[ 4] * temp - 2*cFrac[ 3] + maxT * ( c[ 5] * temp - 3*cFrac[ 4] + maxT * -4*cFrac[ 5] )))) ); // else out of bounds - high


          dE0dT <<= dE0dT + yi
                  * cond( temp <= c[0] && temp >= minT, c[8] + temp * ( c[9] + temp * ( c[10] + temp * ( c[11] + temp * c[12] ))) )  // if low temp
                        ( temp >  c[0] && temp <= maxT, c[1] + temp * ( c[2] + temp * ( c[ 3] + temp * ( c[ 4] + temp * c[ 5] ))) )  // else if high temp
                        ( temp < minT,                  c[8] + minT * ( c[9] + minT * ( c[10] + minT * ( c[11] + minT * c[12] ))) )  // else if out of bounds - low
                        (                               c[1] + maxT * ( c[2] + maxT * ( c[ 3] + maxT * ( c[ 4] + maxT * c[ 5] ))) ); // else out of bounds - high
          break;
        default: {
          std::ostringstream msg;
          msg << __FILE__ << " : " << __LINE__ << "\n Error for spec n = " << n << exceptionMsg_.str();
          throw std::runtime_error( msg.str() );
          }
        } // switch( polyType )
      }
    } // species loop
    // Newton's method to find root
    res  <<= delE0/dE0dT;
    temp <<= temp + res;
#   ifdef ENABLE_CUDA
    res.set_device_as_active(CPU_INDEX);
#   endif
    const double err = field_max_interior( abs(res) );
#   ifdef ENABLE_CUDA
    res.set_device_as_active(GPU_INDEX);
#   endif
    isConverged = err < tol_; // Converged when the temperature has changed by less than set tolerance
    ++iterations;

    if( !isConverged && iterations == maxIterations_ ){
      std::ostringstream msg;
      msg << std::endl
          << "Error in pokitt::TemperatureFromE0::evaluate()." << std::endl
          << "Iteration count exceeded" << std::endl
          << "Total   iterations = " << iterations << std::endl
          << "Maximum iterations = " << maxIterations_ << std::endl
          << "Set tolerance              = " << tol_ << std::endl
          << "Largest pointwise residual = " << err << std::endl
          << __FILE__ << " : " << __LINE__ << std::endl;
      throw std::runtime_error( msg.str() );
    }
  }
}
    
//--------------------------------------------------------------------
template< typename FieldT >
void
TemperatureFromE0<FieldT>::
sensitivity( const Expr::Tag& var)
{
    using namespace SpatialOps;
    FieldT& dfdv = this->sensitivity_result( var );
    if( var == this->get_tag() ){
      dfdv <<= 1.0;
    }
    else if( var == e0_->tag()){
      FieldT& temp = this->value();
      SpatFldPtr<FieldT> cvPtr  = SpatialFieldStore::get<FieldT>(temp); //heat capacity of mixture
      FieldT& cv = *cvPtr;
      cv <<= 0.0;
      const double maxTval = field_max_interior(temp);
      const double minTval = field_min_interior(temp);
      for( size_t n=0; n<nSpec_; ++n ){
          const FieldT& yi = massFracs_[n]->field_ref();
          const ThermData& specTherm = specThermVec_[n];
          const ThermoPoly type = specTherm.type;
          const std::vector<double>& c = specTherm.coefficients;
          const double minT = specTherm.minTemp;
          const double maxT = specTherm.maxTemp;
          #     ifndef ENABLE_CUDA // optimization benefits only the CPU - cond performs betters with if/else than with if/elif/elif/else
          if( maxTval <= maxT && minTval >= minT){ // if true, temperature can only be either high or low
              switch ( type ) {
                  case CONST_POLY: // constant cp
                      cv <<= cv + yi * c[3];
                      break;
                  case NASA_POLY:
                      cv <<= cv + yi
                      * cond( temp <= c[0], c[8] + temp * ( c[9] + temp * ( c[10] + temp * ( c[11] + temp * c[12] ))))   // if low temp
                            (               c[1] + temp * ( c[2] + temp * ( c[ 3] + temp * ( c[ 4] + temp * c[ 5] ))));  // else if high temp;
                      break;
                  default: {
                      std::ostringstream msg;
                      msg << __FILE__ << " : " << __LINE__ << "\n Error for spec n = " << n << exceptionMsg_.str();
                      throw std::runtime_error( msg.str() );
                  }
              } // switch( type )
          }
          else
          #     endif
          {
              /* else temperature can be out of bounds low, low temp, high temp, or out of bounds high
               * if out of bounds, properties are keeped to the same as that for min or max temp
               */
              switch ( type ){
                  case CONST_POLY: // constant cp
                      cv <<= cv + yi * c[3] ;
                      break;
                  case NASA_POLY:
                      cv <<= cv + yi
                      * cond( temp <= c[0] && temp >= minT, c[8] + temp * ( c[9] + temp * ( c[10] + temp * ( c[11] + temp * c[12] ))) )  // if low temp
                            ( temp >  c[0] && temp <= maxT, c[1] + temp * ( c[2] + temp * ( c[ 3] + temp * ( c[ 4] + temp * c[ 5] ))) )  // else if high temp
                            ( temp < minT,                  c[8] + minT * ( c[9] + minT * ( c[10] + minT * ( c[11] + minT * c[12] ))) )  // else if out of bounds - low
                            (                               c[1] + maxT * ( c[2] + maxT * ( c[ 3] + maxT * ( c[ 4] + maxT * c[ 5] ))) ); // else out of bounds - high

                      break;
                  default: {
                      std::ostringstream msg;
                      msg << __FILE__ << " : " << __LINE__ << "\n Error for spec n = " << n << exceptionMsg_.str();
                      throw std::runtime_error( msg.str() );
                  }
              } // switch( type )
          }
      } // species loop
      dfdv <<= 1/cv;
    }
    else{
      dfdv <<= 0.0;
    }
}

//--------------------------------------------------------------------


template< typename FieldT >
void
TemperatureFromE0<FieldT>::
find_bad_points( std::ostringstream& msg, const FieldT& badField, const double badValue, const bool checkBelow )
{
  bool checkNaN = false;
  if( std::isnan( badValue ) ) checkNaN = true;
  typename FieldT::const_iterator iField = badField.begin();
  const SpatialOps::MemoryWindow mw = badField.window_with_ghost();
  int badPoints = 0;
  double worstValue = 0;
  int xWorst, yWorst, zWorst;

  for( size_t z = 1; z <= mw.extent(2); ++z ){
    for( size_t y = 1; y <= mw.extent(1); ++y ){
      for( size_t x = 1; x <= mw.extent(0); ++x, ++iField ){
        if( checkNaN ){
          if( std::isnan( *iField ) ){
            ++badPoints;
            if( badPoints < 6 )
              msg << "Value = " << *iField << " at [x,y,z] = [" << x << "," <<  y << "," << z << "]" << std::endl;
            worstValue = *iField;
            xWorst = x; yWorst = y; zWorst = z;
          }
        }
        else if( checkBelow ){
          if( *iField <= badValue  ){
            ++badPoints;
            if( badPoints < 6 )
              msg << "Value = " << *iField << " at [x,y,z] = [" << x << "," <<  y << "," << z << "]" << std::endl;
          }
          if( *iField < worstValue ){
            worstValue = *iField;
            xWorst = x; yWorst = y; zWorst = z;
          }

        }
        else{
          if( *iField >= badValue ){
            ++badPoints;
            if( badPoints < 6 )
              msg << "Value = " << *iField << " at [x,y,z] = [" << x << "," <<  y << "," << z << "]" << std::endl;
          }
          if( *iField > worstValue ){
            worstValue = *iField;
            xWorst = x; yWorst = y; zWorst = z;
          }
        }
      }
    }
  }
  msg << std::endl
      << "Worst value = " << worstValue << " at [x,y,z] = [" << xWorst << "," <<  yWorst << "," << zWorst << "]" << std::endl
      << "A total of " << badPoints << " bad point(s) detected." << std::endl;
}

//--------------------------------------------------------------------

} // namespace pokitt

#endif // Temperature_Expr_h


