#ifndef Temperature_Expr_h
#define Temperature_Expr_h

#include <expression/Expression.h>

#include <cantera/kernel/ct_defs.h> // contains value of gas constant
#include <cantera/kernel/speciesThermoTypes.h> // contains definitions for which polynomial is being used
#include <cantera/IdealGasMix.h>

#include <pokitt/CanteraObjects.h> // include Cantera wrapper

namespace Cantera_CXX{ class IdealGasMix; } // location of polynomial coefficients

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
  Expr::TagList massFracTags_;
  const FieldT* enth_;
  std::vector<const FieldT*> massFracs_;

  PolyVals minTVec_; // vector of minimum temperatures for polynomial evaluations
  PolyVals maxTVec_; // vector of maximum temperatures for polynomial evaluations
  std::vector< PolyVals > cVec_; // vector of polynomial coefficients
  std::vector<int> polyTypeVec_; // vector of polynomial types
  int nSpec_; // number of species to iterate over
  bool nasaFlag_; // flag if NASA polynomial is present
  bool shomateFlag_; // flag if shomate polynomial is present


  Temperature( const Expr::Tag& massFracTag,
               const Expr::Tag& enthTag );
public:

  static const Expr::TagList& temperature_powers_tags();

  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a Temperature expression
     *  @param resultTag the tag for the temperature of the mixture
     *  @param massFracTag tag for the mass fraction of each species, ordering is consistent with Cantera input file
     *  @param enthTag tag for the enthalpy of the mixture
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::Tag& massFracTag,
             const Expr::Tag& enthTag  );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag enthTag_;
    const Expr::Tag massFracTag_;
  };

  ~Temperature();
  void advertise_dependents( Expr::ExprDeps& exprDeps );
  void bind_fields( const Expr::FieldManagerList& fml );
  void evaluate();
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
  const Expr::Tag e0Tag_;
  const Expr::Tag keTag_;
  Expr::TagList massFracTags_;
  const FieldT* e0_;
  const FieldT* ke_;
  std::vector<const FieldT*> massFracs_;

  PolyVals minTVec_; // vector of minimum temperatures for polynomial evaluations
  PolyVals maxTVec_; // vector of maximum temperatures for polynomial evaluations
  std::vector< PolyVals > cVec_; // vector of polynomial coefficients
  std::vector<int> polyTypeVec_; // vector of polynomial types
  int nSpec_; // number of species to iterate over
  bool nasaFlag_; // flag if NASA polynomial is present
  bool shomateFlag_; // flag if Shomate polynomial is present

  TemperatureFromE0( const Expr::Tag& massFracTag,
                     const Expr::Tag& e0Tag,
                     const Expr::Tag& keTag );
public:

  static const Expr::TagList& temperature_powers_tags();

  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a TemperatureFromE0 expression
     *  @param resultTag the tag for the temperature of the mixture
     *  @param massFracTag tag for the mass fraction of each species, ordering is consistent with Cantera input file
     *  @param e0Tag tag for the total energy of the mixture
     *  @param keTag tag for kinetic energy of the mixture
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::Tag& massFracTag,
             const Expr::Tag& e0Tag,
             const Expr::Tag& keTag );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag e0Tag_;
    const Expr::Tag massFracTag_;
    const Expr::Tag keTag_;
  };

  ~TemperatureFromE0();
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
const Expr::TagList&
Temperature<FieldT>::temperature_powers_tags()
{
  using namespace Expr;
  static TagList tags = tag_list( Tag("T^2"  ,STATE_NONE),
                                  Tag("T^3"  ,STATE_NONE),
                                  Tag("T^4"  ,STATE_NONE),
                                  Tag("T^5"  ,STATE_NONE),
                                  Tag("1/T"  ,STATE_NONE),
                                  Tag("1/T^2",STATE_NONE) );
  return tags;
}


template< typename FieldT >
Temperature<FieldT>::
Temperature( const Expr::Tag& massFracTag,
             const Expr::Tag& enthTag )
  : Expr::Expression<FieldT>(),
    enthTag_( enthTag ),
    nasaFlag_ ( false ),
    shomateFlag_ ( false )
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
  int polyType; // type of polynomial - const_cp, shomate, or NASA
  double minT; // minimum temperature polynomial is valid
  double maxT; // maximum temperature polynomial is valid
  double refPressure;

  for( size_t n=0; n<nSpec_; ++n ){
    spThermo.reportParams(n, polyType, &c[0], minT, maxT, refPressure);
    switch( polyType ){ // check to ensure that we're using a supported polynomial
    case NASA2   : nasaFlag_    = true; break;
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
Temperature<FieldT>::
~Temperature()
{

}

//--------------------------------------------------------------------

template< typename FieldT >
void
Temperature<FieldT>::
advertise_dependents( Expr::ExprDeps& exprDeps )
{
  exprDeps.requires_expression( massFracTags_ );
  exprDeps.requires_expression( enthTag_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
Temperature<FieldT>::
bind_fields( const Expr::FieldManagerList& fml )
{
  const typename Expr::FieldMgrSelector<FieldT>::type& fm = fml.field_manager<FieldT>();
  massFracs_.clear();
  BOOST_FOREACH( const Expr::Tag& tag, massFracTags_ ){
    massFracs_.push_back( &fm.field_ref(tag) );
  }
  enth_ = &fm.field_ref( enthTag_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
Temperature<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  using namespace Cantera;
  SpecT& temp_vec = this->get_value_vec();
  FieldT& temp        = *temp_vec[0];
  FieldT& t2          = *temp_vec[1]; // t^2
  FieldT& t3          = *temp_vec[2]; // t^3
  FieldT& t4          = *temp_vec[3]; // t^4
  FieldT& t5          = *temp_vec[4]; // t^5
  FieldT& recipT      = *temp_vec[5]; // t^-1
  FieldT& recipRecipT = *temp_vec[6]; // t^-2

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
  while( !isConverged ){
    delH <<= *enth_;
    dhdT <<= 0.0;
    // pre-compute powers of temperature used in polynomial evaluations
    if( nasaFlag_ == true || shomateFlag_ == true ){
      t2 <<= temp * temp;
      t3 <<= t2 * temp;
      t4 <<= t3 * temp;
      if( nasaFlag_ == true )
        t5 <<= t4 * temp;
      if( shomateFlag_ == true ){
        recipT <<= 1/ temp;
        recipRecipT <<= recipT * recipT;
      }
    }
    for( size_t n=0; n<nSpec_; ++n ){
      const int polyType = polyTypeVec_[n];
      const std::vector<double>& c = cVec_[n];
      const double minT = minTVec_[n];
      const double maxT = maxTVec_[n];
#     ifndef ENABLE_CUDA // optimization benefits only the CPU - cond performs betters with if/else than with if/elif/elif/else
      if( nebo_max(temp) <= maxT && nebo_min(temp) >= minT){ // if true, temperature can only be either high or low
        switch (polyType) {
        case SIMPLE: // constant cp
          delH <<= delH - *massFracs_[n] * ( c[1] + c[3]*(temp-c[0]) );
          dhdT <<= dhdT + *massFracs_[n] * c[3];
          break;
        case NASA2:
          delH <<= delH - *massFracs_[n] * cond( temp <= c[0], c[ 6] + c[1] * temp + c[2]/2 * t2 + c[ 3]/3 * t3 + c[ 4]/4 * t4 + c[ 5]/5 * t5 ) // if low temp
                                               (               c[13] + c[8] * temp + c[9]/2 * t2 + c[10]/3 * t3 + c[11]/4 * t4 + c[12]/5 * t5 );// else if high temp

          dhdT <<= dhdT + *massFracs_[n] * cond( temp <= c[0], c[1] + c[2] * temp + c[ 3] * t2 + c[ 4] * t3 + c[ 5] * t4) // if low temp
                                               (               c[8] + c[9] * temp + c[10] * t2 + c[11] * t3 + c[12] * t4);// else if high temp
          break;
        case SHOMATE2:
          delH <<= delH - *massFracs_[n] * cond( temp <= c[0], c[ 6] + c[1] * temp*1e-3 + c[2]/2 * t2*1e-6 + c[ 3]/3 * t3*1e-9 + c[ 4]/4 * t4*1e-12 - c[ 5] * recipT*1e3 ) // if low temp
                                               (               c[13] + c[8] * temp*1e-3 + c[9]/2 * t2*1e-6 + c[10]/3 * t3*1e-9 + c[11]/4 * t4*1e-12 - c[12] * recipT*1e3 ); // else if high range
                              // factor of 1e-3 is for units conversion
          dhdT <<= dhdT + *massFracs_[n] * 1e-3 * cond( temp <= c[0], c[1] + c[2] * temp*1e-3 + c[ 3] * t2*1e-6 + c[ 4] * t3*1e-9 + c[ 5] * recipRecipT*1e6 ) // if low temp
                                                      (               c[8] + c[9] * temp*1e-3 + c[10] * t2*1e-6 + c[11] * t3*1e-9 + c[12] * recipRecipT*1e6 ); // else if high range
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
          delH<<=delH - *massFracs_[n] * ( c[1] + c[3]*(temp-c[0]) );
          dhdT<<=dhdT + *massFracs_[n] * c[3];
          break;
        case NASA2:
          delH <<= delH - *massFracs_[n] * cond( temp <= c[0] && temp >= minT, c[ 6] + c[1] * temp + c[2]/2 * t2 + c[ 3]/3 * t3 + c[ 4]/4 * t4 + c[ 5]/5 * t5 ) // if low temp
                                               ( temp >  c[0] && temp <= maxT, c[13] + c[8] * temp + c[9]/2 * t2 + c[10]/3 * t3 + c[11]/4 * t4 + c[12]/5 * t5 )  // else if high temp
                                               ( temp < minT, c[ 6] + c[1] * temp + c[2] * minT * (temp - minT/2) + c[ 3] * minT * minT * (temp - 2*minT/3) + c[ 4]*pow(minT,3) * (temp - 3*minT/4) + c[ 5]*pow(minT,4) * (temp - 4*minT/5) ) // else if out of bounds - low
                                               (              c[13] + c[8] * temp + c[9] * maxT * (temp - maxT/2) + c[10] * maxT * maxT * (temp - 2*maxT/3) + c[11]*pow(maxT,3) * (temp - 3*maxT/4) + c[12]*pow(maxT,4) * (temp - 4*maxT/5) ); // else out of bounds - high

          dhdT <<= dhdT + *massFracs_[n] * cond( temp <= c[0] && temp >= minT, c[1] + c[2] * temp + c[ 3] * t2 + c[ 4] * t3 + c[ 5] * t4) // if low temp
                                               ( temp >  c[0] && temp <= maxT, c[8] + c[9] * temp + c[10] * t2 + c[11] * t3 + c[12] * t4)  // else if high temp
                                               ( temp < minT, c[1] + c[2] * minT + c[ 3] * minT * minT + c[ 4] * pow(minT,3) + c[ 5] * pow(minT,4))  // else if out of bounds - low
                                               (              c[8] + c[9] * maxT + c[10] * maxT * maxT + c[11] * pow(maxT,3) + c[12] * pow(maxT,4)); // else out of bounds - high

          break;
        case SHOMATE2:
          delH <<= delH - *massFracs_[n] * cond( temp <= c[0] && temp >= minT, c[ 6] + c[1] * temp*1e-3 + c[2]/2 * t2*1e-6 + c[ 3]/3 * t3*1e-9 + c[ 4]/4 * t4*1e-12 - c[ 5] * recipT*1e3 ) // if low temp
                                               ( temp >  c[0] && temp <= maxT, c[13] + c[8] * temp*1e-3 + c[9]/2 * t2*1e-6 + c[10]/3 * t3*1e-9 + c[11]/4 * t4*1e-12 - c[12] * recipT*1e3 )  // else if high temp
                                               ( temp < minT, c[1] * temp*1e-3 + c[2] * minT*1e-3 * (temp*1e-3 - minT*1e-3/2) + c[ 3] * minT*1e-3 * minT*1e-3 * (temp*1e-3 - 2*minT*1e-3/3) + c[ 4]*pow(minT*1e-3,3) * (temp*1e-3 - 3*minT*1e-3/4) - c[ 5]*pow(minT*1e-3,-1) * (-temp*1e-3 / minT*1e-3 + 2) + c[ 6] ) // else if out of bounds - low
                                               (              c[8] * temp*1e-3 + c[9] * maxT*1e-3 * (temp*1e-3 - maxT*1e-3/2) + c[10] * maxT*1e-3 * maxT*1e-3 * (temp*1e-3 - 2*maxT*1e-3/3) + c[11]*pow(maxT*1e-3,3) * (temp*1e-3 - 3*maxT*1e-3/4) - c[12]*pow(maxT*1e-3,-1) * (-temp*1e-3 / maxT*1e-3 + 2) + c[13] ); // else out of bounds - high
                              // factor of 1e-3 is for units conversion
          dhdT <<= dhdT + *massFracs_[n] * 1e-3 * cond( temp <= c[0] && temp >= minT, c[1] + c[2] * temp*1e-3 + c[ 3] * t2*1e-6 + c[ 4] * t3*1e-9 + c[ 5] * recipRecipT*1e6 ) // if low temp
                                                      ( temp >  c[0] && temp <= maxT, c[8] + c[9] * temp*1e-3 + c[10] * t2*1e-6 + c[11] * t3*1e-9 + c[12] * recipRecipT*1e6 )  // else if high temp
                                                      ( temp < minT, c[1] + c[2] * minT*1e-3 + c[ 3] * minT*1e-3 * minT*1e-3 + c[ 4] * pow(minT*1e-3,3) + c[ 5] * pow(minT*1e-3,-2) )  // else if out of bounds - low
                                                      (              c[8] + c[9] * maxT*1e-3 + c[10] * maxT*1e-3 * maxT*1e-3 + c[11] * pow(maxT*1e-3,3) + c[12] * pow(maxT*1e-3,-2) ); // else out of bounds - high

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
    isConverged = err < 1e-4; // Converged when the temperature has changed by less than 1e-4
  }
  t2          <<= temp * temp;
  t3          <<= t2 * temp;
  t4          <<= t3 * temp;
  t5          <<= t4 * temp;
  recipT      <<= 1/ temp;
  recipRecipT <<= recipT * recipT;
}


//--------------------------------------------------------------------

template< typename FieldT >
Temperature<FieldT>::
Builder::Builder( const Expr::Tag& resultTag,
                  const Expr::Tag& massFracTag,
                  const Expr::Tag& enthTag )
: ExpressionBuilder( tag_list( resultTag, Temperature<FieldT>::temperature_powers_tags() ) ),
  massFracTag_ ( massFracTag ),
  enthTag_( enthTag )
{}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
Temperature<FieldT>::
Builder::build() const
{
  return new Temperature<FieldT>( massFracTag_, enthTag_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
const Expr::TagList&
TemperatureFromE0<FieldT>::temperature_powers_tags()
{
  using namespace Expr;
  static const TagList tags = tag_list( Tag("T^2"  ,STATE_NONE),
                                        Tag("T^3"  ,STATE_NONE),
                                        Tag("T^4"  ,STATE_NONE),
                                        Tag("T^5"  ,STATE_NONE),
                                        Tag("1/T"  ,STATE_NONE),
                                        Tag("1/T^2",STATE_NONE) );
  return tags;
}

template< typename FieldT >
TemperatureFromE0<FieldT>::
TemperatureFromE0( const Expr::Tag& massFracTag,
                   const Expr::Tag& e0Tag,
                   const Expr::Tag& keTag )
  : Expr::Expression<FieldT>(),
    e0Tag_( e0Tag ),
    keTag_( keTag ),
    nasaFlag_ ( false ),
    shomateFlag_ ( false )
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
  int polyType; // type of polynomial - const_cp, shomate, or NASA
  double minT; // minimum temperature polynomial is valid
  double maxT; // maximum temperature polynomial is valid
  double refPressure;
  for( size_t n=0; n<nSpec_; ++n ){
    spThermo.reportParams(n, polyType, &c[0], minT, maxT, refPressure);
    switch( polyType ){ // check to ensure that we are using a supported polynomial type
    case NASA2   : nasaFlag_    = true; break;
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
      for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic)
        *ic *= Cantera::GasConstant / molecularWeights[n]; // dimensionalize the coefficients to mass basis
      break;
    case SHOMATE2:
      c[1] -= Cantera::GasConstant * 1e-3; //change coefficients from isobaric to isometric
      c[8] -= Cantera::GasConstant * 1e-3;
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
TemperatureFromE0<FieldT>::
~TemperatureFromE0()
{

}

//--------------------------------------------------------------------

template< typename FieldT >
void
TemperatureFromE0<FieldT>::
advertise_dependents( Expr::ExprDeps& exprDeps )
{
  exprDeps.requires_expression( massFracTags_ );
  exprDeps.requires_expression( e0Tag_ );
  exprDeps.requires_expression( keTag_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
TemperatureFromE0<FieldT>::
bind_fields( const Expr::FieldManagerList& fml )
{
  const typename Expr::FieldMgrSelector<FieldT>::type& fm = fml.field_manager<FieldT>();
  massFracs_.clear();
  BOOST_FOREACH( const Expr::Tag& tag, massFracTags_ ){
    massFracs_.push_back( &fm.field_ref(tag) );
  }
  e0_ = &fm.field_ref( e0Tag_ );
  ke_ = &fm.field_ref( keTag_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
TemperatureFromE0<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  using namespace Cantera;
  SpecT& temp_vec = this->get_value_vec();
  FieldT& temp        = *temp_vec[0];
  FieldT& t2          = *temp_vec[1]; // t^2
  FieldT& t3          = *temp_vec[2]; // t^3
  FieldT& t4          = *temp_vec[3]; // t^4
  FieldT& t5          = *temp_vec[4]; // t^5
  FieldT& recipT      = *temp_vec[5]; // t^-1
  FieldT& recipRecipT = *temp_vec[6]; // t^-2

  SpatFldPtr<FieldT> delE0Ptr = SpatialFieldStore::get<FieldT>(temp); // difference between internal energy field value and internal energy at current temperature
  SpatFldPtr<FieldT> dE0dTPtr = SpatialFieldStore::get<FieldT>(temp); // dE0dT for Newton's method
  SpatFldPtr<FieldT> resPtr   = SpatialFieldStore::get<FieldT>(temp); // change in temperature for new Newton's iteration

  FieldT& delE0 = *delE0Ptr;
  FieldT& dE0dT = *dE0dTPtr;
  FieldT& res   = *resPtr;

  res.add_device(CPU_INDEX);
  bool isConverged = false;
  while( !isConverged ){
    delE0 <<= *e0_ - *ke_;
    dE0dT <<= 0.0;
    // pre-compute powers of temperature used in polynomial evaluations
    if( nasaFlag_ == true | shomateFlag_ == true ){
      t2 <<= temp * temp;
      t3 <<= t2 * temp;
      t4 <<= t3 * temp;
      if( nasaFlag_ == true )
        t5<<=t4 * temp;
      if( shomateFlag_ == true ){
        recipT <<= 1/ temp;
        recipRecipT <<= recipT * recipT;
      }
    }
    for( size_t n=0; n<nSpec_; ++n){
      const int polyType = polyTypeVec_[n];
      const std::vector<double>& c = cVec_[n];
      const double minT = minTVec_[n];
      const double maxT = maxTVec_[n];
#     ifndef ENABLE_CUDA // optimization benefits only the CPU - cond performs betters with if/else than with if/elif/elif/else
      if( nebo_max(temp) <= maxT && nebo_min(temp) >= minT){ // if true, temperature can only be either high or low
        switch (polyType) {
          case SIMPLE: // constant cv
            delE0 <<= delE0 - *massFracs_[n] * ( c[1] + c[3]*(temp-c[0]) );
            dE0dT <<= dE0dT + *massFracs_[n] * c[3];
            break;
          case NASA2:
            delE0 <<= delE0 - *massFracs_[n] * cond( temp <= c[0], c[ 6] + c[1] * temp + c[2]/2 * t2 + c[ 3]/3 * t3 + c[ 4]/4 * t4 + c[ 5]/5 * t5 ) // if low temp
                                                   (               c[13] + c[8] * temp + c[9]/2 * t2 + c[10]/3 * t3 + c[11]/4 * t4 + c[12]/5 * t5 );// else if high temp

            dE0dT <<= dE0dT + *massFracs_[n] * cond( temp <= c[0], c[1] + c[2] * temp + c[ 3] * t2 + c[ 4] * t3 + c[ 5] * t4) // if low temp
                                                   (               c[8] + c[9] * temp + c[10] * t2 + c[11] * t3 + c[12] * t4);// else if high temp
            break;
          case SHOMATE2:
            delE0 <<= delE0 - *massFracs_[n] * cond( temp <= c[0], c[ 6] + c[1] * temp*1e-3 + c[2]/2 * t2*1e-6 + c[ 3]/3 * t3*1e-9 + c[ 4]/4 * t4*1e-12 - c[ 5] * recipT*1e3 ) // if low temp
                                                   (               c[13] + c[8] * temp*1e-3 + c[9]/2 * t2*1e-6 + c[10]/3 * t3*1e-9 + c[11]/4 * t4*1e-12 - c[12] * recipT*1e3 ); // else if high range
                                            // 1e-3 is a factor for units conversion
            dE0dT <<= dE0dT + *massFracs_[n] * 1e-3 * cond( temp <= c[0], c[1] + c[2] * temp*1e-3 + c[ 3] * t2*1e-6 + c[ 4] * t3*1e-9 + c[ 5] * recipRecipT*1e6 ) // if low temp
                                                          (               c[8] + c[9] * temp*1e-3 + c[10] * t2*1e-6 + c[11] * t3*1e-9 + c[12] * recipRecipT*1e6 ); // else if high range
            break;
        } // switch( polyType )
      }
    else
#   endif
      {
      /* else temperature can be out of bounds low, low temp, high temp, or out of bounds high
       * if out of bounds, properties are interpolated from min or max temp using a constant cv
       */
        switch (polyType) {
          case SIMPLE: // constant cv
            delE0<<=delE0 - *massFracs_[n] * ( c[1] + c[3] * (temp-c[0]) );
            dE0dT<<=dE0dT + *massFracs_[n] * c[3];
            break;
          case NASA2:
            delE0 <<= delE0 - *massFracs_[n] * cond( temp <= c[0] && temp >= minT, c[ 6] + c[1] * temp + c[2]/2 * t2 + c[ 3]/3 * t3 + c[ 4]/4 * t4 + c[ 5]/5 * t5 ) // if low temp
                                                   ( temp >  c[0] && temp <= maxT, c[13] + c[8] * temp + c[9]/2 * t2 + c[10]/3 * t3 + c[11]/4 * t4 + c[12]/5 * t5 )  // else if high temp
                                                   ( temp < minT, c[1] * temp + c[2] * minT * (temp - minT/2) + c[ 3] * minT * minT * (temp - 2*minT/3) + c[ 4]*pow(minT,3) * (temp - 3*minT/4) + c[ 5]*pow(minT,4) * (temp - 4*minT/5) + c[ 6] ) // else if out of bounds - low
                                                   (              c[8] * temp + c[9] * maxT * (temp - maxT/2) + c[10] * maxT * maxT * (temp - 2*maxT/3) + c[11]*pow(maxT,3) * (temp - 3*maxT/4) + c[12]*pow(maxT,4) * (temp - 4*maxT/5) + c[13] ); // else out of bounds - high


            dE0dT <<= dE0dT + *massFracs_[n] * cond( temp <= c[0] && temp >= minT, c[1] + c[2] * temp + c[ 3] * t2 + c[ 4] * t3 + c[ 5] * t4 ) // if low temp
                                                     ( temp >  c[0] && temp <= maxT, c[8] + c[9] * temp + c[10] * t2 + c[11] * t3 + c[12] * t4 )  // else if high temp
                                                     ( temp < minT, c[1] + c[2] * minT + c[ 3] * minT * minT + c[ 4] * pow(minT,3) + c[ 5] * pow(minT,4) )  // else if out of bounds - low
                                                     (              c[8] + c[9] * maxT + c[10] * maxT * maxT + c[11] * pow(maxT,3) + c[12] * pow(maxT,4) ); // else out of bounds - high

            break;
          case SHOMATE2:
            delE0 <<= delE0 - *massFracs_[n] * cond( temp <= c[0] && temp >= minT, c[ 6] + c[1] * temp*1e-3 + c[2]/2 * t2*1e-6 + c[ 3]/3 * t3*1e-9 + c[ 4]/4 * t4*1e-12 - c[ 5] * recipT*1e3 ) // if low temp
                                                   ( temp >  c[0] && temp <= maxT, c[13] + c[8] * temp*1e-3 + c[9]/2 * t2*1e-6 + c[10]/3 * t3*1e-9 + c[11]/4 * t4*1e-12 - c[12] * recipT*1e3 )  // else if high temp
                                                   ( temp < minT, c[1] * temp*1e-3 + c[2] * minT*1e-3 * (temp*1e-3 - minT*1e-3/2) + c[ 3] * minT*1e-3 * minT*1e-3 * (temp*1e-3 - 2*minT*1e-3/3) + c[ 4]*pow(minT*1e-3,3) * (temp*1e-3 - 3*minT*1e-3/4) - c[ 5]*pow(minT*1e-3,-1) * (-temp*1e-3 / minT*1e-3 + 2) + c[ 6] ) // else if out of bounds - low
                                                   (              c[8] * temp*1e-3 + c[9] * maxT*1e-3 * (temp*1e-3 - maxT*1e-3/2) + c[10] * maxT*1e-3 * maxT*1e-3 * (temp*1e-3 - 2*maxT*1e-3/3) + c[11]*pow(maxT*1e-3,3) * (temp*1e-3 - 3*maxT*1e-3/4) - c[12]*pow(maxT*1e-3,-1) * (-temp*1e-3 / maxT*1e-3 + 2) + c[13] ); // else out of bounds - high
                                            // 1e-3 is a factor for units conversion
            dE0dT <<= dE0dT + *massFracs_[n] * 1e-3 * cond( temp <= c[0] && temp >= minT, c[1] + c[2] * temp*1e-3 + c[ 3] * t2*1e-6 + c[ 4] * t3*1e-9 + c[ 5] * recipRecipT*1e6 ) // if low temp
                                                          ( temp >  c[0] && temp <= maxT, c[8] + c[9] * temp*1e-3 + c[10] * t2*1e-6 + c[11] * t3*1e-9 + c[12] * recipRecipT*1e6 )  // else if high temp
                                                          ( temp < minT, c[1] + c[2] * minT*1e-3 + c[ 3] * minT*1e-3 * minT*1e-3 + c[ 4] * pow(minT*1e-3,3) + c[ 5] * pow(minT*1e-3,-2) )  // else if out of bounds - low
                                                          (              c[8] + c[9] * maxT*1e-3 + c[10] * maxT*1e-3 * maxT*1e-3 + c[11] * pow(maxT*1e-3,3) + c[12] * pow(maxT*1e-3,-2) ); // else out of bounds - high

        } // switch( polyType )
      }
    } // species loop
    // Newton's method to find root
    res  <<= delE0/dE0dT;
    temp <<= temp + res;
    res.set_device_as_active(CPU_INDEX);
    const double err = nebo_max( abs(res) );
#   ifdef ENABLE_CUDA
    res.set_device_as_active(GPU_INDEX);
#   endif
    isConverged = err < 1e-3; // Converged when the temperature has changed by less than 1e-3
  }
  t2          <<= temp * temp;
  t3          <<= t2 * temp;
  t4          <<= t3 * temp;
  t5          <<= t4 * temp;
  recipT      <<= 1/ temp;
  recipRecipT <<= recipT * recipT;
}


//--------------------------------------------------------------------

template< typename FieldT >
TemperatureFromE0<FieldT>::
Builder::Builder( const Expr::Tag& resultTag,
                  const Expr::Tag& massFracTag,
                  const Expr::Tag& e0Tag,
                  const Expr::Tag& keTag )
: ExpressionBuilder( tag_list( resultTag, TemperatureFromE0<FieldT>::temperature_powers_tags()) ),
  massFracTag_ ( massFracTag ),
  e0Tag_( e0Tag ),
  keTag_( keTag )
{}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
TemperatureFromE0<FieldT>::
Builder::build() const
{
  return new TemperatureFromE0<FieldT>( massFracTag_, e0Tag_, keTag_ );
}

#endif // Temperature_Expr_h
