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

#ifndef Temperature_Expr_h
#define Temperature_Expr_h

#include <expression/Expression.h>

#include <cantera/kernel/ct_defs.h> // contains value of gas constant
#include <cantera/kernel/speciesThermoTypes.h> // contains definitions for which polynomial is being used

#include <pokitt/CanteraObjects.h> // include Cantera wrapper

#include <math.h> // isnan error checking

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
 * given h_0 and an initial guess \f$ T_i \f$ Newton's method calculates
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
 */

template< typename FieldT >
class Temperature
    : public Expr::Expression<FieldT>
{
  typedef std::vector<FieldT*> SpecT;
  typedef std::vector<double> PolyVals; // values used for polynomial
  const Expr::Tag enthTag_;
  const Expr::TagList massFracTags_;

  DECLARE_FIELD( FieldT, enth_ )
  DECLARE_VECTOR_OF_FIELDS( FieldT, massFracs_ )

  PolyVals minTVec_; // vector of minimum temperatures for polynomial evaluations
  PolyVals maxTVec_; // vector of maximum temperatures for polynomial evaluations
  std::vector< PolyVals > cVec_; // vector of polynomial coefficients
  std::vector< PolyVals > cFracVec_; // vector of polynomial coefficients divided by integers
  std::vector<int> polyTypeVec_; // vector of polynomial types
  int nSpec_; // number of species to iterate over
  bool shomateFlag_; // flag if shomate polynomial is present
  const double tol_; // tolerance for Newton's method


  Temperature( const Expr::TagList& massFracTags,
               const Expr::Tag& enthTag,
               const double tol );
public:

  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a Temperature expression
     *  @param resultTag the tag for the temperature of the mixture
     *  @param massFracTags tag for the mass fraction of each species, ordering is consistent with Cantera input file
     *  @param enthTag tag for the enthalpy of the mixture
     *  @param tol     tolerance for the non-linear solve
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::TagList& massFracTags,
             const Expr::Tag& enthTag,
             const double tol = 1e-3,
             const int nghost = DEFAULT_NUMBER_OF_GHOSTS );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag enthTag_;
    const Expr::TagList massFracTags_;
    const double tol_;
  };

  ~Temperature(){}
  void evaluate();
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
 * given e_0 and an initial guess \f$ T_i \f$ Newton's method calculates
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
 */

template< typename FieldT >
class TemperatureFromE0
    : public Expr::Expression<FieldT>
{
  typedef std::vector<FieldT*> SpecT;
  typedef std::vector<double> PolyVals; // values used for polynomial

  DECLARE_FIELDS( FieldT, e0_, ke_ )
  DECLARE_VECTOR_OF_FIELDS( FieldT, massFracs_ )

  PolyVals minTVec_; // vector of minimum temperatures for polynomial evaluations
  PolyVals maxTVec_; // vector of maximum temperatures for polynomial evaluations
  std::vector< PolyVals > cVec_; // vector of polynomial coefficients
  std::vector< PolyVals > cFracVec_; // vector of polynomial coefficients divided by integers
  std::vector<int> polyTypeVec_; // vector of polynomial types
  int nSpec_; // number of species to iterate over
  bool shomateFlag_; // flag if Shomate polynomial is present
  const double tol_; // tolerance for Newton's method

  TemperatureFromE0( const Expr::TagList& massFracTags,
                     const Expr::Tag& e0Tag,
                     const Expr::Tag& keTag,
                     const double tol );
public:

  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a TemperatureFromE0 expression
     *  @param resultTag the tag for the temperature of the mixture
     *  @param massFracTags tag for the mass fraction of each species, ordering is consistent with Cantera input file
     *  @param e0Tag tag for the total energy of the mixture
     *  @param keTag tag for kinetic energy of the mixture
     *  @param tol   tolerance for the non-linear solve
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::TagList& massFracTags,
             const Expr::Tag& e0Tag,
             const Expr::Tag& keTag,
             const double tol = 1e-3,
             const int nghost = DEFAULT_NUMBER_OF_GHOSTS );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag e0Tag_;
    const Expr::TagList massFracTags_;
    const Expr::Tag keTag_;
    const double tol_;
  };

  ~TemperatureFromE0(){}
  void evaluate();
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
             const double tol )
  : Expr::Expression<FieldT>(),
    shomateFlag_ ( false ),
    tol_( tol )
{
  this->set_gpu_runnable( true );

  enth_ = this->template create_field_request<FieldT>( enthTag );
  this->template create_field_vector_request<FieldT>( massFracTags, massFracs_ );

  Cantera_CXX::IdealGasMix* const gasMix = CanteraObjects::get_gasmix();
  const Cantera::SpeciesThermo& spThermo = gasMix->speciesThermo();

  nSpec_ = gasMix->nSpecies();
  const std::vector<double> molecularWeights = gasMix->molecularWeights();
  std::vector<double> c(15,0); //vector of Cantera's coefficients
  std::vector<double> cFrac(15,0); //vector of Cantera's coefficients divided by integers
  int polyType; // type of polynomial - const_cp, shomate, or NASA
  double minT; // minimum temperature polynomial is valid
  double maxT; // maximum temperature polynomial is valid
  double refPressure;

  for( size_t n=0; n<nSpec_; ++n ){
    spThermo.reportParams(n, polyType, &c[0], minT, maxT, refPressure);
    cFrac = c;
    switch( polyType ){ // check to ensure that we're using a supported polynomial
    case NASA2   :                      break;
    case SHOMATE2: shomateFlag_ = true; break;
    case SIMPLE  :                      break;
    default:{
      std::ostringstream msg;
      msg << __FILE__ << " : " << __LINE__
          << "\nThermo type not supported,\n Type = " << polyType
          << ", species # " << n << std::endl;
      throw std::invalid_argument( msg.str() );
    }
    }
    polyTypeVec_.push_back(polyType); // vector of polynomial types
    minTVec_.push_back(minT); // vector of minimum temperatures for polynomial evaluations
    maxTVec_.push_back(maxT); // vector of maximum temperatures for polynomial evaluations
    switch (polyType) {
    case SIMPLE:
      c[1] /= molecularWeights[n]; // convert to mass basis
      c[3] /= molecularWeights[n]; // convert to mass basis
      break;
    case NASA2:
      for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic){
        *ic *= Cantera::GasConstant / molecularWeights[n]; // dimensionalize the coefficients to mass basis
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
    case SHOMATE2:
      for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic ){
        *ic *= 1e6 / molecularWeights[n]; // scale the coefficients to keep units consistent on mass basis
      }
      cFrac[ 2] = c[ 2] / 2;
      cFrac[ 3] = c[ 3] / 3;
      cFrac[ 4] = c[ 4] / 4;
      cFrac[ 9] = c[ 9] / 2;
      cFrac[10] = c[10] / 3;
      cFrac[11] = c[11] / 4;
      break;
    }
    cVec_.push_back(c); // vector of polynomial coefficients
    cFracVec_.push_back(cFrac);
  }
  CanteraObjects::restore_gasmix(gasMix);
}

//--------------------------------------------------------------------

template< typename FieldT >
void
Temperature<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  using namespace Cantera;
  FieldT& temp = this->value();

  const FieldT& enth = enth_->field_ref();

  SpatFldPtr<FieldT> recipT;      // may be used for Shomate polynomial
  SpatFldPtr<FieldT> recipRecipT; // may be used for Shomate polynomial
  if( shomateFlag_ ) {
    recipT      = SpatialFieldStore::get<FieldT>(temp);
    recipRecipT = SpatialFieldStore::get<FieldT>(temp);
  }

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
  const int maxIts = 20;
  const double exceptionTemp = 5000; // maximum temperature before exception is thrown
  while( !isConverged ){
    delH <<= enth;
    dhdT <<= 0.0;
    // pre-compute powers of temperature used in polynomial evaluations
    if( shomateFlag_ == true ){
      *recipT      <<= 1/ temp;
      *recipRecipT <<= *recipT * *recipT;
    }
#   ifndef ENABLE_CUDA
    const double maxTval = nebo_max(temp);
    const double minTval = nebo_min(temp);
    if( maxTval >= exceptionTemp ){
      std::ostringstream msg;
        msg << std::endl
            << "Error in pokitt::Temperature::evaluate()." << std::endl
            << "Temperature is too high" << std::endl;
        find_bad_points( msg, temp, exceptionTemp, false );
        msg << "Total   iterations = " << iterations << std::endl
            << "Maximum iterations = " << maxIts << std::endl
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
            << "Maximum iterations = " << maxIts << std::endl
            << "Set tolerance      = " << tol_ << std::endl
            << __FILE__ << " : " << __LINE__ << std::endl;
        throw std::runtime_error( msg.str() );
    }
    if( isnan( nebo_sum_interior( temp ) ) ){
      std::ostringstream msg;
        msg << std::endl
            << "Error in pokitt::Temperature::evaluate()." << std::endl
            << "Temperature is NaN" << std::endl;
        find_bad_points( msg, temp, NAN, false );
        msg << "Total   iterations = " << iterations << std::endl
            << "Maximum iterations = " << maxIts << std::endl
            << "Set tolerance      = " << tol_ << std::endl
            << __FILE__ << " : " << __LINE__ << std::endl;
        throw std::runtime_error( msg.str() );
    }
#   endif
    for( size_t n=0; n<nSpec_; ++n ){
      const FieldT& yi = massFracs_[n]->field_ref();
      const int polyType = polyTypeVec_[n];
      const std::vector<double>& c = cVec_[n];
      const std::vector<double>& cFrac = cFracVec_[n];
      const double minT = minTVec_[n];
      const double maxT = maxTVec_[n];
#     ifndef ENABLE_CUDA // optimization benefits only the CPU - cond performs betters with if/else than with if/elif/elif/else
      if( maxTval <= maxT && minTval >= minT){ // if true, temperature can only be either high or low
        switch (polyType) {
        case SIMPLE: // constant cp
          delH <<= delH - yi * ( c[1] + c[3]*(temp-c[0]) );
          dhdT <<= dhdT + yi * c[3];
          break;
        case NASA2:
          delH <<= delH - yi
                 * cond( temp <= c[0], c[ 6] + temp * ( c[1] + temp * ( cFrac[2] + temp * ( cFrac[ 3] + temp * ( cFrac[ 4] + temp * cFrac[ 5] )))) )  // if low temp
                       (               c[13] + temp * ( c[8] + temp * ( cFrac[9] + temp * ( cFrac[10] + temp * ( cFrac[11] + temp * cFrac[12] )))) );  // else if high temp

          dhdT <<= dhdT + yi
                 * cond( temp <= c[0], c[1] + temp * ( c[2] + temp * ( c[ 3] + temp * ( c[ 4] + temp * c[ 5] ))) )  // if low temp
                       (               c[8] + temp * ( c[9] + temp * ( c[10] + temp * ( c[11] + temp * c[12] ))) );  // else if high temp
          break;
        case SHOMATE2:
          delH <<= delH - yi
                 * cond( temp <= c[0], c[ 6] + temp*1e-3 * ( c[1] + temp*1e-3 * ( cFrac[2] + temp*1e-3 * ( cFrac[ 3] + temp*1e-3 * cFrac[ 4] ))) - c[ 5] * *recipT*1e3 ) // if low temp
                       (               c[13] + temp*1e-3 * ( c[8] + temp*1e-3 * ( cFrac[9] + temp*1e-3 * ( cFrac[10] + temp*1e-3 * cFrac[11] ))) - c[12] * *recipT*1e3 );  // else if high range

          dhdT <<= dhdT + yi * 1e-3 // 1e-3 is a factor for units conversion
                 * cond( temp <= c[0], c[1] + temp*1e-3 * ( c[2] + temp*1e-3 * ( c[ 3] + temp*1e-3 * c[ 4] )) + c[ 5] * *recipRecipT*1e6 )  // if low temp
                       (               c[8] + temp*1e-3 * ( c[9] + temp*1e-3 * ( c[10] + temp*1e-3 * c[11] )) + c[12] * *recipRecipT*1e6 );  // else if high temp
          break;
        } // switch( polyType )
      }
    else
#   endif
      {
        /* else temperature can be out of bounds low, low temp, high temp, or out of bounds high
         * if out of bounds, properties are interpolated from min or max temp using a constant cp
         */
        switch (polyType){
        case SIMPLE: // constant cp
          delH<<=delH - yi * ( c[1] + c[3]*(temp-c[0]) );
          dhdT<<=dhdT + yi * c[3];
          break;
        case NASA2:
          delH <<= delH - yi
                 * cond( temp <= c[0] && temp >= minT, c[ 6] + temp * ( c[1] + temp * ( cFrac[2] + temp * ( cFrac[ 3] + temp * ( cFrac[ 4] + temp * cFrac[ 5] )))) )  // if low temp
                       ( temp >  c[0] && temp <= maxT, c[13] + temp * ( c[8] + temp * ( cFrac[9] + temp * ( cFrac[10] + temp * ( cFrac[11] + temp * cFrac[12] )))) )  // else if high temp
                       ( temp < minT,                  c[ 6] + c[1] * temp + minT * ( c[2] * temp + minT * ( c[ 3] * temp - cFrac[2] + minT * ( c[ 4] * temp - 2*cFrac[ 3] + minT * ( c[ 5] * temp - 3*cFrac[ 4] + minT * -4*cFrac[ 5] )))) )  // else if out of bounds - low
                       (                               c[13] + c[8] * temp + maxT * ( c[9] * temp + maxT * ( c[10] * temp - cFrac[9] + maxT * ( c[11] * temp - 2*cFrac[10] + maxT * ( c[12] * temp - 3*cFrac[11] + maxT * -4*cFrac[12] )))) ); // else out of bounds - high

          dhdT <<= dhdT + yi
                 * cond( temp <= c[0] && temp >= minT, c[1] + temp * ( c[2] + temp * ( c[ 3] + temp * ( c[ 4] + temp * c[ 5] ))) )  // if low temp
                       ( temp >  c[0] && temp <= maxT, c[8] + temp * ( c[9] + temp * ( c[10] + temp * ( c[11] + temp * c[12] ))) )  // else if high temp
                       ( temp < minT,                  c[1] + minT * ( c[2] + minT * ( c[ 3] + minT * ( c[ 4] + minT * c[ 5] ))) )  // else if out of bounds - low
                       (                               c[8] + maxT * ( c[9] + maxT * ( c[10] + maxT * ( c[11] + maxT * c[12] ))) ); // else out of bounds - high

          break;
        case SHOMATE2:
          delH <<= delH - yi
                 * cond( temp <= c[0] && temp >= minT, c[ 6] + temp*1e-3 * ( c[1] + temp*1e-3 * ( cFrac[2] + temp*1e-3 * ( cFrac[ 3] + temp*1e-3 * cFrac[ 4] ))) - c[ 5] * *recipT*1e3 ) // if low temp
                       ( temp >  c[0] && temp <= maxT, c[13] + temp*1e-3 * ( c[8] + temp*1e-3 * ( cFrac[9] + temp*1e-3 * ( cFrac[10] + temp*1e-3 * cFrac[11] ))) - c[12] * *recipT*1e3 )  // else if high range
                       ( temp <  minT,                 c[ 6] + c[1] * temp*1e-3 + minT*1e-3 * ( c[2] * temp*1e-3 + minT*1e-3 * ( c[ 3] * temp*1e-3 - cFrac[2] + minT*1e-3 * ( c[ 4] * temp*1e-3 - 2*cFrac[ 3] + minT*1e-3 * -3*cFrac[ 4] ))) + ( c[ 5] * temp / minT - 2*c[ 5] ) / (minT*1e-3) ) // else if out of bounds - low
                       (                               c[13] + c[8] * temp*1e-3 + maxT*1e-3 * ( c[9] * temp*1e-3 + maxT*1e-3 * ( c[10] * temp*1e-3 - cFrac[9] + maxT*1e-3 * ( c[11] * temp*1e-3 - 2*cFrac[10] + maxT*1e-3 * -3*cFrac[11] ))) + ( c[12] * temp / maxT - 2*c[12] ) / (maxT*1e-3) ); // else out of bounds - high

          dhdT <<= dhdT + yi * 1e-3 // 1e-3 is a factor for units conversion
                 * cond( temp <= c[0] && temp >= minT, c[1] + temp*1e-3 * ( c[2] + temp*1e-3 * ( c[ 3] + temp*1e-3 * c[ 4] )) + c[ 5] * *recipRecipT*1e6  )  // if low temp
                       ( temp >  c[0] && temp <= maxT, c[8] + temp*1e-3 * ( c[9] + temp*1e-3 * ( c[10] + temp*1e-3 * c[11] )) + c[12] * *recipRecipT*1e6  )  // else if high temp
                       ( temp < minT,                  c[1] + minT*1e-3 * ( c[2] + minT*1e-3 * ( c[ 3] + minT*1e-3 * c[ 4] )) + c[ 5] / (minT*minT*1e-6) )  // else if out of bounds - low
                       (                               c[8] + maxT*1e-3 * ( c[9] + maxT*1e-3 * ( c[10] + maxT*1e-3 * c[11] )) + c[12] / (maxT*maxT*1e-6) ); // else out of bounds - high

        } // switch( polyType )
      }
    } // species loop
    // Newton's method to find root
    res  <<= delH/dhdT;
    temp <<= temp + res;
#   ifdef ENABLE_CUDA
    res.set_device_as_active(CPU_INDEX);
#   endif
    const double err = nebo_max( abs(res) );
#   ifdef ENABLE_CUDA
    res.set_device_as_active(GPU_INDEX);
#   endif
    isConverged = err < tol_; // Converged when the temperature has changed by less than specified tolerance
    ++iterations;

    if( !isConverged && iterations == maxIts ){
      std::ostringstream msg;
        msg << std::endl
            << "Error in pokitt::Temperature::evaluate()." << std::endl
            << "Iteration count exceeded" << std::endl
            << "Total   iterations = " << iterations << std::endl
            << "Maximum iterations = " << maxIts << std::endl
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
find_bad_points( std::ostringstream& msg, const FieldT& badField, const double badValue, const bool checkBelow )
{
  bool checkNaN = false;
  if( isnan( badValue ) ) checkNaN = true;
  typename FieldT::const_iterator iField = badField.begin();
  const SpatialOps::MemoryWindow mw = badField.window_with_ghost();
  int badPoints = 0;
  double worstValue = 0;
  int xWorst, yWorst, zWorst;

  for( int z = 1; z <= mw.extent(2); ++z ){
    for( int y = 1; y <= mw.extent(1); ++y ){
      for( int x = 1; x <= mw.extent(0); ++x, ++iField ){
        if( checkNaN ){
          if( isnan( *iField ) ){
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
Temperature<FieldT>::
Builder::Builder( const Expr::Tag& resultTag,
                  const Expr::TagList& massFracTags,
                  const Expr::Tag& enthTag,
                  const double tol,
                  const int nghost )
: ExpressionBuilder( resultTag, nghost ),
  massFracTags_( massFracTags ),
  enthTag_( enthTag ),
  tol_( tol )
{}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
Temperature<FieldT>::
Builder::build() const
{
  return new Temperature<FieldT>( massFracTags_, enthTag_, tol_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
TemperatureFromE0<FieldT>::
TemperatureFromE0( const Expr::TagList& massFracTags,
                   const Expr::Tag& e0Tag,
                   const Expr::Tag& keTag,
                   const double tol )
  : Expr::Expression<FieldT>(),
    tol_( tol ),
    shomateFlag_ ( false )
{
  this->set_gpu_runnable( true );

  e0_ = this->template create_field_request<FieldT>( e0Tag );
  ke_ = this->template create_field_request<FieldT>( keTag );
  this->template create_field_vector_request<FieldT>( massFracTags, massFracs_ );

  Cantera_CXX::IdealGasMix* const gasMix = CanteraObjects::get_gasmix();
  const Cantera::SpeciesThermo& spThermo = gasMix->speciesThermo();

  nSpec_ = gasMix->nSpecies();

  const std::vector<double> molecularWeights = gasMix->molecularWeights();
  std::vector<double> c(15,0); //vector of Cantera's coefficients
  std::vector<double> cFrac(15,0); //vector of Cantera's coefficients divided by integers
  int polyType; // type of polynomial - const_cp, shomate, or NASA
  double minT; // minimum temperature polynomial is valid
  double maxT; // maximum temperature polynomial is valid
  double refPressure;
  for( size_t n=0; n<nSpec_; ++n ){
    spThermo.reportParams(n, polyType, &c[0], minT, maxT, refPressure);
    cFrac = c;
    switch( polyType ){ // check to ensure that we are using a supported polynomial type
    case NASA2   :                      break;
    case SHOMATE2: shomateFlag_ = true; break;
    case SIMPLE  :                      break;
    default:{
      std::ostringstream msg;
      msg << __FILE__ << " : " << __LINE__
          << "\nThermo type not supported,\n Type = " << polyType
          << ", species # " << n << std::endl;
      throw std::invalid_argument( msg.str() );
    }
    }
    polyTypeVec_.push_back(polyType); // vector of polynomial types
    minTVec_.push_back(minT); // vector of minimum temperatures for polynomial evaluations
    maxTVec_.push_back(maxT); // vector of maximum temperatures for polynomial evaluations
    switch (polyType) {
    case SIMPLE:
      c[3] -= Cantera::GasConstant; // change coefficients from isobaric to isometric
      c[1] -= Cantera::GasConstant * c[0];
      c[3] /= molecularWeights[n]; // convert to mass basis
      c[1] /= molecularWeights[n];
      break;
    case NASA2:
      c[1] -= 1.0; //change coefficients from isobaric to isometric
      c[8] -= 1.0;
      for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic){
        *ic *= Cantera::GasConstant / molecularWeights[n]; // dimensionalize the coefficients to mass basis
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
    case SHOMATE2:
      c[1] -= Cantera::GasConstant * 1e-3; //change coefficients from isobaric to isometric
      c[8] -= Cantera::GasConstant * 1e-3;
      for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic ){
        *ic *= 1e6 / molecularWeights[n]; // scale the coefficients to keep units consistent on mass basis
      }
      cFrac[ 2] = c[ 2] / 2;
      cFrac[ 3] = c[ 3] / 3;
      cFrac[ 4] = c[ 4] / 4;
      cFrac[ 9] = c[ 9] / 2;
      cFrac[10] = c[10] / 3;
      cFrac[11] = c[11] / 4;
      break;
    }
    cVec_.push_back(c); // vector of polynomial coefficients
    cFracVec_.push_back(cFrac);
  }
  CanteraObjects::restore_gasmix(gasMix);
}

//--------------------------------------------------------------------

template< typename FieldT >
void
TemperatureFromE0<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  using namespace Cantera;
  FieldT& temp = this->value();

  const FieldT& e0 = e0_->field_ref();
  const FieldT& ke = ke_->field_ref();

  SpatFldPtr<FieldT> recipT;      // may be used for Shomate polynomial
  SpatFldPtr<FieldT> recipRecipT; // may be used for Shomate polynomial
  if( shomateFlag_ ) {
    recipT      = SpatialFieldStore::get<FieldT>(temp);
    recipRecipT = SpatialFieldStore::get<FieldT>(temp);
  }

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
  const int maxIts = 20;
  const double exceptionTemp = 5000; // maximum temperature before exception is thrown
    while( !isConverged ){
      delE0 <<= e0 - ke;
      dE0dT <<= 0.0;
      // pre-compute powers of temperature used in polynomial evaluations
      if( shomateFlag_ == true ){
        *recipT      <<= 1/ temp;
        *recipRecipT <<= *recipT * *recipT;
      }
#     ifndef ENABLE_CUDA
      const double maxTval = nebo_max(temp);
      const double minTval = nebo_min(temp);
      if( maxTval >= exceptionTemp ){
        std::ostringstream msg;
          msg << std::endl
              << "Error in pokitt::TemperatureFromE0::evaluate()." << std::endl
              << "Temperature is too high" << std::endl;
          find_bad_points( msg, temp, exceptionTemp, false );
          msg << "Total   iterations = " << iterations << std::endl
              << "Maximum iterations = " << maxIts << std::endl
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
              << "Maximum iterations = " << maxIts << std::endl
              << "Set tolerance      = " << tol_ << std::endl
              << __FILE__ << " : " << __LINE__ << std::endl;
          throw std::runtime_error( msg.str() );
      }
      if( isnan( nebo_sum( temp ) ) ){
        std::ostringstream msg;
          msg << std::endl
              << "Error in pokitt::TemperatureFromE0::evaluate()." << std::endl
              << "Temperature is NaN" << std::endl;
          find_bad_points( msg, temp, NAN, false );
          msg << "Total   iterations = " << iterations << std::endl
              << "Maximum iterations = " << maxIts << std::endl
              << "Set tolerance      = " << tol_ << std::endl
              << __FILE__ << " : " << __LINE__ << std::endl;
          throw std::runtime_error( msg.str() );
      }
#     endif
      for( size_t n=0; n<nSpec_; ++n){
        const FieldT& yi = massFracs_[n]->field_ref();
        const int polyType = polyTypeVec_[n];
        const std::vector<double>& c = cVec_[n];
        const std::vector<double>& cFrac = cFracVec_[n];
        const double minT = minTVec_[n];
        const double maxT = maxTVec_[n];
#       ifndef ENABLE_CUDA // optimization benefits only the CPU - cond performs betters with if/else than with if/elif/elif/else
        if( maxTval <= maxT && minTval >= minT){ // if true, temperature can only be either high or low
          switch (polyType) {
          case SIMPLE: // constant cv
            delE0 <<= delE0 - yi * ( c[1] + c[3]*(temp-c[0]) );
            dE0dT <<= dE0dT + yi * c[3];
            break;
          case NASA2:
            delE0 <<= delE0 - yi
                    * cond( temp <= c[0], c[ 6] + temp * ( c[1] + temp * ( cFrac[2] + temp * ( cFrac[ 3] + temp * ( cFrac[ 4] + temp * cFrac[ 5] )))) )  // if low temp
                          (               c[13] + temp * ( c[8] + temp * ( cFrac[9] + temp * ( cFrac[10] + temp * ( cFrac[11] + temp * cFrac[12] )))) );  // else if high temp

            dE0dT <<= dE0dT + yi
                    * cond( temp <= c[0], c[1] + temp * ( c[2] + temp * ( c[ 3] + temp * ( c[ 4] + temp * c[ 5] ))) )  // if low temp
                          (               c[8] + temp * ( c[9] + temp * ( c[10] + temp * ( c[11] + temp * c[12] ))) );  // else if high temp
            break;
          case SHOMATE2:
            delE0 <<= delE0 - yi
                    * cond( temp <= c[0], c[ 6] + temp*1e-3 * ( c[1] + temp*1e-3 * ( cFrac[2] + temp*1e-3 * ( cFrac[ 3] + temp*1e-3 * cFrac[ 4] ))) - c[ 5] * *recipT*1e3 ) // if low temp
                          (               c[13] + temp*1e-3 * ( c[8] + temp*1e-3 * ( cFrac[9] + temp*1e-3 * ( cFrac[10] + temp*1e-3 * cFrac[11] ))) - c[12] * *recipT*1e3 );  // else if high range

            dE0dT <<= dE0dT + yi * 1e-3 // 1e-3 is a factor for units conversion
                    * cond( temp <= c[0], c[1] + temp*1e-3 * ( c[2] + temp*1e-3 * ( c[ 3] + temp*1e-3 * c[ 4])) + c[ 5] * *recipRecipT*1e6  )  // if low temp
                          (               c[8] + temp*1e-3 * ( c[9] + temp*1e-3 * ( c[10] + temp*1e-3 * c[11])) + c[12] * *recipRecipT*1e6  );  // else if high temp
            break;
          } // switch( polyType )
        }
        else
#       endif
        {
          /* else temperature can be out of bounds low, low temp, high temp, or out of bounds high
           * if out of bounds, properties are interpolated from min or max temp using a constant cv
           */
          switch (polyType) {
          case SIMPLE: // constant cv
            delE0 <<= delE0 - yi * ( c[1] + c[3] * (temp-c[0]) );
            dE0dT <<= dE0dT + yi * c[3];
            break;
          case NASA2:
            delE0 <<= delE0 - yi
                    * cond( temp <= c[0] && temp >= minT, c[ 6] + temp * ( c[1] + temp * ( cFrac[2] + temp * ( cFrac[ 3] + temp * ( cFrac[ 4] + temp * cFrac[ 5] )))) )  // if low temp
                          ( temp >  c[0] && temp <= maxT, c[13] + temp * ( c[8] + temp * ( cFrac[9] + temp * ( cFrac[10] + temp * ( cFrac[11] + temp * cFrac[12] )))) )  // else if high temp
                          ( temp < minT,                  c[ 6] + c[1] * temp + minT * ( c[2] * temp + minT * ( c[ 3] * temp - cFrac[2] + minT * ( c[ 4] * temp - 2*cFrac[ 3] + minT * ( c[ 5] * temp - 3*cFrac[ 4] + minT * -4*cFrac[ 5] )))) )  // else if out of bounds - low
                          (                               c[13] + c[8] * temp + maxT * ( c[9] * temp + maxT * ( c[10] * temp - cFrac[9] + maxT * ( c[11] * temp - 2*cFrac[10] + maxT * ( c[12] * temp - 3*cFrac[11] + maxT * -4*cFrac[12] )))) ); // else out of bounds - high


            dE0dT <<= dE0dT + yi
                    * cond( temp <= c[0] && temp >= minT, c[1] + temp * ( c[2] + temp * ( c[ 3] + temp * ( c[ 4] + temp * c[ 5] ))) )  // if low temp
                          ( temp >  c[0] && temp <= maxT, c[8] + temp * ( c[9] + temp * ( c[10] + temp * ( c[11] + temp * c[12] ))) )  // else if high temp
                          ( temp < minT,                  c[1] + minT * ( c[2] + minT * ( c[ 3] + minT * ( c[ 4] + minT * c[ 5] ))) )  // else if out of bounds - low
                          (                               c[8] + maxT * ( c[9] + maxT * ( c[10] + maxT * ( c[11] + maxT * c[12] ))) ); // else out of bounds - high


            break;
          case SHOMATE2:
            delE0 <<= delE0 - yi
                    * cond( temp <= c[0] && temp >= minT, c[ 6] + temp*1e-3 * ( c[1] + temp*1e-3 * ( cFrac[2] + temp*1e-3 * ( cFrac[ 3] + temp*1e-3 * cFrac[ 4] ))) - c[ 5] * *recipT*1e3 ) // if low temp
                          ( temp >  c[0] && temp <= maxT, c[13] + temp*1e-3 * ( c[8] + temp*1e-3 * ( cFrac[9] + temp*1e-3 * ( cFrac[10] + temp*1e-3 * cFrac[11] ))) - c[12] * *recipT*1e3 )  // else if high range
                          ( temp <  minT,                 c[ 6] + c[1] * temp*1e-3 + minT*1e-3 * ( c[2] * temp*1e-3 + minT*1e-3 * ( c[ 3] * temp*1e-3 - cFrac[2] + minT*1e-3 * ( c[ 4] * temp*1e-3 - 2*cFrac[ 3] + minT*1e-3 * -3*cFrac[ 4] ))) + ( c[ 5] * temp / minT - 2*c[ 5] ) / (minT*1e-3) ) // else if out of bounds - low
                          (                               c[13] + c[8] * temp*1e-3 + maxT*1e-3 * ( c[9] * temp*1e-3 + maxT*1e-3 * ( c[10] * temp*1e-3 - cFrac[9] + maxT*1e-3 * ( c[11] * temp*1e-3 - 2*cFrac[10] + maxT*1e-3 * -3*cFrac[11] ))) + ( c[12] * temp / maxT - 2*c[12] ) / (maxT*1e-3) ); // else out of bounds - high

            dE0dT <<= dE0dT + yi * 1e-3 // 1e-3 is a factor for units conversion
                    * cond( temp <= c[0] && temp >= minT, c[1] + temp*1e-3 * ( c[2] + temp*1e-3 * ( c[ 3] + temp*1e-3 * c[ 4])) + c[ 5] * *recipRecipT*1e6  )  // if low temp
                          ( temp >  c[0] && temp <= maxT, c[8] + temp*1e-3 * ( c[9] + temp*1e-3 * ( c[10] + temp*1e-3 * c[11])) + c[12] * *recipRecipT*1e6  )  // else if high temp
                          ( temp < minT,                  c[1] + minT*1e-3 * ( c[2] + minT*1e-3 * ( c[ 3] + minT*1e-3 * c[ 4])) + c[ 5] / (minT*minT*1e-6) )  // else if out of bounds - low
                          (                               c[8] + maxT*1e-3 * ( c[9] + maxT*1e-3 * ( c[10] + maxT*1e-3 * c[11])) + c[12] / (maxT*maxT*1e-6) ); // else out of bounds - high
            break;
          } // switch( polyType )
        }
      } // species loop
      // Newton's method to find root
      res  <<= delE0/dE0dT;
      temp <<= temp + res;
#     ifdef ENABLE_CUDA
      res.set_device_as_active(CPU_INDEX);
#     endif
      const double err = nebo_max( abs(res) );
#     ifdef ENABLE_CUDA
      res.set_device_as_active(GPU_INDEX);
#     endif
      isConverged = ( err < tol_ ); // Converged when the temperature has changed by less than set tolerance
      ++iterations;

      if( !isConverged && iterations == maxIts ){
        std::ostringstream msg;
          msg << std::endl
              << "Error in pokitt::TemperatureFromE0::evaluate()." << std::endl
              << "Iteration count exceeded" << std::endl
              << "Total   iterations = " << iterations << std::endl
              << "Maximum iterations = " << maxIts << std::endl
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
find_bad_points( std::ostringstream& msg, const FieldT& badField, const double badValue, const bool checkBelow )
{
  bool checkNaN = false;
  if( isnan( badValue ) ) checkNaN = true;
  typename FieldT::const_iterator iField = badField.begin();
  const SpatialOps::MemoryWindow mw = badField.window_with_ghost();
  int badPoints = 0;
  double worstValue = 0;
  int xWorst, yWorst, zWorst;

  for( int z = 1; z <= mw.extent(2); ++z ){
    for( int y = 1; y <= mw.extent(1); ++y ){
      for( int x = 1; x <= mw.extent(0); ++x, ++iField ){
        if( checkNaN ){
          if( isnan( *iField ) ){
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
Builder::Builder( const Expr::Tag& resultTag,
                  const Expr::TagList& massFracTags,
                  const Expr::Tag& e0Tag,
                  const Expr::Tag& keTag,
                  const double tol,
                  const int nghost )
: ExpressionBuilder( resultTag, nghost ),
  massFracTags_( massFracTags ),
  e0Tag_( e0Tag ),
  keTag_( keTag ),
  tol_( tol )
{}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
TemperatureFromE0<FieldT>::
Builder::build() const
{
  return new TemperatureFromE0<FieldT>( massFracTags_, e0Tag_, keTag_, tol_ );
}

} // namespace pokitt

#endif // Temperature_Expr_h

