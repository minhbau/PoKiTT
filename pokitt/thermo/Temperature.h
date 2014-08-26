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
  const Expr::Tag enthTag_;
  Expr::TagList massFracTags_;
  const FieldT* enth_;
  std::vector<const FieldT*> massFracs_;
  Cantera_CXX::IdealGasMix* const gasMix_; // gas mixture to be evaluated
  const int nSpec_; // number of species to iterate over
  bool nasaFlag_; // flag if NASA polynomial is present
  bool shomateFlag_; // flag if shomate polynomial is present


  Temperature( const Expr::Tag& massFracTag,
               const Expr::Tag& enthTag );
public:
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
 * \mid T_{i+1} - T_{i} \mid < 10^{-3}
 * \f]
 */

template< typename FieldT >
class TemperatureFromE0
    : public Expr::Expression<FieldT>
{
  const Expr::Tag e0Tag_;
  const Expr::Tag keTag_;
  Expr::TagList massFracTags_;
  const FieldT* e0_;
  const FieldT* ke_;
  std::vector<const FieldT*> massFracs_;
  Cantera_CXX::IdealGasMix* const gasMix_; // gas mixture to be evaluated
  const int nSpec_; // number of species to iterate over
  bool nasaFlag_; // flag if NASA polynomial is present
  bool shomateFlag_; // flag if Shomate polynomial is present

  TemperatureFromE0( const Expr::Tag& massFracTag,
                     const Expr::Tag& e0Tag,
                     const Expr::Tag& keTag );
public:
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
Temperature<FieldT>::
Temperature( const Expr::Tag& massFracTag,
             const Expr::Tag& enthTag )
  : Expr::Expression<FieldT>(),
    enthTag_( enthTag ),
    gasMix_( CanteraObjects::get_gasmix() ),
    nSpec_( gasMix_->nSpecies() ),
    nasaFlag_ ( false ),
    shomateFlag_ ( false )
{
  // not GPU runnable until we get reductions supported on the GPU.
//  this->set_gpu_runnable( true );

  massFracTags_.clear();
  for( size_t n=0; n<nSpec_; ++n ){
    std::ostringstream name;
    name << massFracTag.name() << "_" << n;
    massFracTags_.push_back( Expr::Tag(name.str(),massFracTag.context()) );
  }

  const Cantera::SpeciesThermo& spThermo = gasMix_->speciesThermo();
  const int nSpec = gasMix_->nSpecies();
  for( size_t n=0; n<nSpec; ++n ){
    const int type = spThermo.reportType(n);
    switch( type ){
      case NASA2   : nasaFlag_    = true; break;
      case SHOMATE2: shomateFlag_ = true; break;
      case SIMPLE  :                      break;
      default:{
        std::ostringstream msg;
        msg << __FILE__ << " : " << __LINE__
            << "\nThermo type not supported,\n Type = " << type
            << ", species # " << n << std::endl;
        throw std::invalid_argument( msg.str() );
      }
    }
  }
      }

//--------------------------------------------------------------------

template< typename FieldT >
Temperature<FieldT>::
~Temperature()
{
  CanteraObjects::restore_gasmix(gasMix_);
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
  FieldT& temp = this->value();

  SpatFldPtr<FieldT> delHPtr  = SpatialFieldStore::get<FieldT>(temp); // difference between enthalpy field value and enthalpy evaluated at current temperature
  SpatFldPtr<FieldT> dhdTPtr  = SpatialFieldStore::get<FieldT>(temp); // dhdT for Newton's method
  SpatFldPtr<FieldT> resPtr   = SpatialFieldStore::get<FieldT>(temp); // change in temperature for new Newton's iteration
  FieldT& delH = *delHPtr;
  FieldT& res = *resPtr;
  FieldT& dhdT = *dhdTPtr;

  SpatFldPtr<FieldT> t2; // t^2
  SpatFldPtr<FieldT> t3; // t^3
  SpatFldPtr<FieldT> t4; // t^4
  SpatFldPtr<FieldT> t5; // t^5
  SpatFldPtr<FieldT> recipT; // t^-1
  SpatFldPtr<FieldT> recipRecipT; // t^-2


  if( nasaFlag_ == true || shomateFlag_ == true){

    t2 = SpatialFieldStore::get<FieldT>(temp);
    t3 = SpatialFieldStore::get<FieldT>(temp);
    t4 = SpatialFieldStore::get<FieldT>(temp);

    if( nasaFlag_ == true ){
      t5 = SpatialFieldStore::get<FieldT>(temp);
    }

    if( shomateFlag_ == true ){
      recipT = SpatialFieldStore::get<FieldT>(temp);
      recipRecipT = SpatialFieldStore::get<FieldT>(temp);
    }

  }

  std::vector<double> c(15,0); //vector of Cantera's polynomial coefficients
  int polyType; // type of polynomial
  double minT; // minimum temperature where polynomial is valid
  double maxT; // maximum temperature where polynomial is valid
  double refPressure;
  const Cantera::SpeciesThermo& spThermo = gasMix_->speciesThermo();

  const std::vector<double>& molecularWeights = gasMix_->molecularWeights();

  bool isConverged = false;
  while( !isConverged ){
    delH <<= *enth_;
    dhdT <<= 0.0;

    // pre-compute powers of temperature used in polynomial evaluations
    if( nasaFlag_ == true || shomateFlag_ == true ){
      *t2 <<= temp * temp;
      *t3 <<= *t2 * temp;
      *t4 <<= *t3 * temp;
      if( nasaFlag_ == true )
        *t5 <<= *t4 * temp;
      if( shomateFlag_ == true ){
        *recipT <<= 1/ temp;
        *recipRecipT <<= *recipT * *recipT;
      }
    }

    for( size_t n=0; n<nSpec_; ++n ){
      spThermo.reportParams(n, polyType, &c[0], minT, maxT, refPressure);
      std::vector<double>::iterator ic = c.begin() + 1;
      std::vector<double>::iterator icend = c.end();
      if( nebo_max(temp) <= maxT && nebo_min(temp) >= minT ){
        /* if true, temperature can only be either high or low
         * This reduces the number of checks and reduces evaluation time
         */
        switch (polyType){
        case SIMPLE: // constant cp
          delH <<= delH - *massFracs_[n] * ( c[1] + c[3] * (temp-c[0]) )
                                         / molecularWeights[n];

          dhdT <<= dhdT + *massFracs_[n] * c[3]
                                         / molecularWeights[n];
          break;
        case NASA2:
          for( ; ic != icend; ++ic)
            *ic *= GasConstant / molecularWeights[n]; // dimensionalize the coefficients
          delH <<= delH - *massFracs_[n] * cond( temp <= c[0] , c[ 6] + c[1] * temp + c[2]/2 * *t2 + c[ 3]/3 * *t3 + c[ 4]/4 * *t4 + c[ 5]/5 * *t5) // if low temp
                                               (                c[13] + c[8] * temp + c[9]/2 * *t2 + c[10]/3 * *t3 + c[11]/4 * *t4 + c[12]/5 * *t5); // else if high temp


          dhdT <<= dhdT + *massFracs_[n] * cond( temp <= c[0] , c[1] + c[2] * temp + c[ 3] * *t2 + c[ 4] * *t3 + c[ 5] * *t4) // if low temp
                                               (                c[8] + c[9] * temp + c[10] * *t2 + c[11] * *t3 + c[12] * *t4); // else if high temp

          break;
        case SHOMATE2:
          for( ; ic != icend; ++ic)
            *ic *= 1e6 / molecularWeights[n]; // scale the coefficients to keep units consistent
          delH<<= delH - *massFracs_[n] * cond( temp <= c[0] , c[ 6] + c[1] * temp*1e-3 + c[2]/2 * *t2*1e-6 + c[ 3]/3 * *t3*1e-9 + c[ 4]/4 * *t4*1e-12 - c[ 5] * *recipT*1e3 ) // if low temp
                                              (                c[13] + c[8] * temp*1e-3 + c[9]/2 * *t2*1e-6 + c[10]/3 * *t3*1e-9 + c[11]/4 * *t4*1e-12 - c[12] * *recipT*1e3 ); // else if high temp


          for( ic = c.begin() + 1; ic != icend; ++ic)
            *ic *= 1e-3; // scale the coefficients again to keep units consistent
          dhdT<<= dhdT + *massFracs_[n] * cond( temp <= c[0] , c[1] + c[2] * temp*1e-3 + c[ 3] * *t2*1e-6 + c[ 4] * *t3*1e-9 + c[ 5] * *recipRecipT*1e6) // if low temp
                                              (                c[8] + c[9] * temp*1e-3 + c[10] * *t2*1e-6 + c[11] * *t3*1e-9 + c[12] * *recipRecipT*1e6); // else if high temp

          break;
        }
      }
      else{
        /* else temperature can be out of bounds low, low temp, high temp, or out of bounds high
         * if out of bounds, properties are interpolated from min or max temp using a constant cp
         */
        switch (polyType){
        case SIMPLE: // constant cp
          delH<<=delH - *massFracs_[n] * ( c[1] + c[3]*(temp-c[0]) )
                                       / molecularWeights[n];

          dhdT<<=dhdT + *massFracs_[n] * c[3]
                                       / molecularWeights[n];
          break;
        case NASA2:
          for( ; ic != icend; ++ic)
            *ic *= GasConstant / molecularWeights[n]; // dimensionalize the coefficients
          delH <<= delH - *massFracs_[n] * cond( temp <= c[0] && temp >= minT, c[ 6] + c[1] * temp + c[2]/2 * *t2 + c[ 3]/3 * *t3 + c[ 4]/4 * *t4 + c[ 5]/5 * *t5 ) // if low temp
                                               ( temp >  c[0] && temp <= maxT, c[13] + c[8] * temp + c[9]/2 * *t2 + c[10]/3 * *t3 + c[11]/4 * *t4 + c[12]/5 * *t5 )  // else if high temp
                                               ( temp < minT, c[ 6] + c[1] * temp + c[2] * minT * (temp - minT/2) + c[ 3] * minT * minT * (temp - 2*minT/3) + c[ 4]*pow(minT,3) * (temp - 3*minT/4) + c[ 5]*pow(minT,4) * (temp - 4*minT/5) ) // else if out of bounds - low
                                               (              c[13] + c[8] * temp + c[9] * maxT * (temp - maxT/2) + c[10] * maxT * maxT * (temp - 2*maxT/3) + c[11]*pow(maxT,3) * (temp - 3*maxT/4) + c[12]*pow(maxT,4) * (temp - 4*maxT/5) ); // else out of bounds - high

          dhdT <<= dhdT + *massFracs_[n] * cond( temp <= c[0] && temp >= minT, c[1] + c[2] * temp + c[ 3] * *t2 + c[ 4] * *t3 + c[ 5] * *t4) // if low temp
                                               ( temp >  c[0] && temp <= maxT, c[8] + c[9] * temp + c[10] * *t2 + c[11] * *t3 + c[12] * *t4)  // else if high temp
                                               ( temp < minT, c[1] + c[2] * minT + c[ 3] * minT * minT + c[ 4] * pow(minT,3) + c[ 5] * pow(minT,4))  // else if out of bounds - low
                                               (              c[8] + c[9] * maxT + c[10] * maxT * maxT + c[11] * pow(maxT,3) + c[12] * pow(maxT,4)); // else out of bounds - high

          break;
        case SHOMATE2:
          double minTScaled = minT/1000;
          double maxTScaled = maxT/1000;
          for( ; ic != icend; ++ic)
            *ic *= 1e6 / molecularWeights[n] ; // scale the coefficients to keep units consistent
          delH <<= delH - *massFracs_[n] * cond( temp <= c[0] && temp >= minT, c[ 6] + c[1] * temp*1e-3 + c[2]/2 * *t2*1e-6 + c[ 3]/3 * *t3*1e-9 + c[ 4]/4 * *t4*1e-12 - c[ 5] * *recipT*1e3 ) // if low temp
                                               ( temp >  c[0] && temp <= maxT, c[13] + c[8] * temp*1e-3 + c[9]/2 * *t2*1e-6 + c[10]/3 * *t3*1e-9 + c[11]/4 * *t4*1e-12 - c[12] * *recipT*1e3 )  // else if high temp
                                               ( temp < minT, c[1] * temp*1e-3 + c[2] * minTScaled * (temp*1e-3 - minTScaled/2) + c[ 3] * minTScaled * minTScaled * (temp*1e-3 - 2*minTScaled/3) + c[ 4]*pow(minTScaled,3) * (temp*1e-3 - 3*minTScaled/4) - c[ 5]*pow(minTScaled,-1) * (-temp*1e-3 / minTScaled + 2) + c[ 6] ) // else if out of bounds - low
                                               (              c[8] * temp*1e-3 + c[9] * maxTScaled * (temp*1e-3 - maxTScaled/2) + c[10] * maxTScaled * maxTScaled * (temp*1e-3 - 2*maxTScaled/3) + c[11]*pow(maxTScaled,3) * (temp*1e-3 - 3*maxTScaled/4) - c[12]*pow(maxTScaled,-1) * (-temp*1e-3 / maxTScaled + 2) + c[13] ); // else out of bounds - high


          for( ic = c.begin() + 1; ic != icend; ++ic)
            *ic *= 1e-3; // scale the coefficients again to keep units consistent
          dhdT <<= dhdT + *massFracs_[n] * cond( temp <= c[0] && temp >= minT, c[1] + c[2] * temp*1e-3 + c[ 3] * *t2*1e-6 + c[ 4] * *t3*1e-9 + c[ 5] * *recipRecipT*1e6 ) // if low temp
                                               ( temp >  c[0] && temp <= maxT, c[8] + c[9] * temp*1e-3 + c[10] * *t2*1e-6 + c[11] * *t3*1e-9 + c[12] * *recipRecipT*1e6 )  // else if high temp
                                               ( temp < minT, c[1] + c[2] * minTScaled + c[ 3] * minTScaled * minTScaled + c[ 4] * pow(minTScaled,3) + c[ 5] * pow(minTScaled,-2) )  // else if out of bounds - low
                                               (              c[8] + c[9] * maxTScaled + c[10] * maxTScaled * maxTScaled + c[11] * pow(maxTScaled,3) + c[12] * pow(maxTScaled,-2) ); // else out of bounds - high

        } // switch
      }
    }
    // Newton's method to find root
    res <<= delH/dhdT;
    temp <<= temp + res;
    isConverged = nebo_max( abs(res) ) < 1e-3; // Converged when the temperature has changed by less than 1e-3
  }
}


//--------------------------------------------------------------------

template< typename FieldT >
Temperature<FieldT>::
Builder::Builder( const Expr::Tag& resultTag,
                  const Expr::Tag& massFracTag,
                  const Expr::Tag& enthTag )
: ExpressionBuilder( resultTag ),
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
TemperatureFromE0<FieldT>::
TemperatureFromE0( const Expr::Tag& massFracTag,
                   const Expr::Tag& e0Tag,
                   const Expr::Tag& keTag )
  : Expr::Expression<FieldT>(),
    e0Tag_( e0Tag ),
    keTag_( keTag ),
    gasMix_( CanteraObjects::get_gasmix() ),
    nSpec_( gasMix_->nSpecies() ),
    nasaFlag_ ( false ),
    shomateFlag_ ( false )
{
  this->set_gpu_runnable( true );

  massFracTags_.clear();
  for( size_t n=0; n<nSpec_; ++n ){
    std::ostringstream name;
    name << massFracTag.name() << "_" << n;
    massFracTags_.push_back( Expr::Tag(name.str(),massFracTag.context()) );
  }

  const Cantera::SpeciesThermo& spThermo = gasMix_->speciesThermo();
  const int nSpec = gasMix_->nSpecies();
  for( size_t n=0; n<nSpec; ++n ){
    const int type = spThermo.reportType(n);
    switch( type ){
      case NASA2   : nasaFlag_    = true; break;
      case SHOMATE2: shomateFlag_ = true; break;
      case SIMPLE  :                      break;
      default:{
        std::ostringstream msg;
        msg << __FILE__ << " : " << __LINE__
            << "\nThermo type not supported,\n Type = " << type
            << ", species # " << n << std::endl;
        throw std::invalid_argument( msg.str() );
      }
    }
  }
}

//--------------------------------------------------------------------

template< typename FieldT >
TemperatureFromE0<FieldT>::
~TemperatureFromE0()
{
  CanteraObjects::restore_gasmix(gasMix_);
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
  FieldT& temp = this->value();

  SpatFldPtr<FieldT> delE0Ptr = SpatialFieldStore::get<FieldT>(temp); // difference between internal energy field value and internal energy at current temperature
  SpatFldPtr<FieldT> dE0dTPtr = SpatialFieldStore::get<FieldT>(temp); // dE0dT for Newton's method
  SpatFldPtr<FieldT> resPtr   = SpatialFieldStore::get<FieldT>(temp); // change in temperature for new Newton's iteration

  FieldT& delE0 = *delE0Ptr;
  FieldT& res   = *resPtr;
  FieldT& dE0dT = *dE0dTPtr;

  SpatFldPtr<FieldT> t2; // t^2
  SpatFldPtr<FieldT> t3; // t^3
  SpatFldPtr<FieldT> t4; // t^4
  SpatFldPtr<FieldT> t5; // t^5
  SpatFldPtr<FieldT> recipT; // t^-1
  SpatFldPtr<FieldT> recipRecipT; // t^-2

  if( nasaFlag_ || shomateFlag_ ){
    t2 = SpatialFieldStore::get<FieldT>(temp);
    t3 = SpatialFieldStore::get<FieldT>(temp);
    t4 = SpatialFieldStore::get<FieldT>(temp);
    if( nasaFlag_ )
      t5 = SpatialFieldStore::get<FieldT>(temp);
    if( shomateFlag_ ){
      recipT = SpatialFieldStore::get<FieldT>(temp);
      recipRecipT = SpatialFieldStore::get<FieldT>(temp);
    }
  }

  std::vector<double> c(15,0); //vector of Cantera's coefficients
  int polyType; // type of polynomial
  double minT; // minimum temperature where polynomial is valid
  double maxT; // maximum temperature where polynomial is valid
  double refPressure;
  const Cantera::SpeciesThermo& spThermo = gasMix_->speciesThermo();
  const std::vector<double>& molecularWeights = gasMix_->molecularWeights();

  bool isConverged = false;

  while( !isConverged ){
    delE0 <<= *e0_ - *ke_;
    dE0dT <<= 0.0;

    // pre-compute powers of temperature used in polynomial evaluations
    if( nasaFlag_ == true | shomateFlag_ == true ){
      *t2 <<= temp * temp;
      *t3 <<= *t2 * temp;
      *t4 <<= *t3 * temp;
      if( nasaFlag_ == true )
        *t5<<=*t4 * temp;
      if( shomateFlag_ == true ){
        *recipT <<= 1/ temp;
        *recipRecipT <<= *recipT * *recipT;
      }
    }

    for( size_t n=0; n<nSpec_; ++n){
      spThermo.reportParams(n, polyType, &c[0], minT, maxT, refPressure);
      std::vector<double>::iterator ic = c.begin() + 1;
      std::vector<double>::iterator icend = c.end();
      if( nebo_max(temp) <= maxT && nebo_min(temp) >= minT){ // if true, temperature can only be either high or low
        switch (polyType) {

          case SIMPLE: // constant cv
            c[3] -= GasConstant; // change coefficients from enthalpy to internal energy
            c[1] -= GasConstant * c[0];
            delE0<<=delE0 - *massFracs_[n] * ( c[1] + c[3]*(temp-c[0]) ) / molecularWeights[n];
            dE0dT<<=dE0dT + *massFracs_[n] * c[3] / molecularWeights[n];
            break;

          case NASA2:
            c[1] -= 1.0; //change coefficients from enthalpy to internal energy
            c[8] -= 1.0;
            for( ; ic != icend; ++ic)
              *ic *= GasConstant / molecularWeights[n]; // dimensionalize the coefficients
            delE0 <<= delE0 - *massFracs_[n] * cond( temp <= c[0], c[ 6] + c[1] * temp + c[2]/2 * *t2 + c[ 3]/3 * *t3 + c[ 4]/4 * *t4 + c[ 5]/5 * *t5 ) // if low temp
                                                   (               c[13] + c[8] * temp + c[9]/2 * *t2 + c[10]/3 * *t3 + c[11]/4 * *t4 + c[12]/5 * *t5 );// else if high temp

            dE0dT <<= dE0dT + *massFracs_[n] * cond( temp <= c[0], c[1] + c[2] * temp + c[ 3] * *t2 + c[ 4] * *t3 + c[ 5] * *t4) // if low temp
                                                   (               c[8] + c[9] * temp + c[10] * *t2 + c[11] * *t3 + c[12] * *t4);// else if high temp
            break;
          case SHOMATE2:
            c[1] -= GasConstant * 1e-3; //change coefficients from enthalpy to internal energy
            c[8] -= GasConstant * 1e-3;
            for( ; ic != icend; ++ic)
              *ic *= 1e6 / molecularWeights[n]; // scale the coefficients again to keep units consistent
            delE0 <<= delE0 - *massFracs_[n] * cond( temp <= c[0], c[ 6] + c[1] * temp*1e-3 + c[2]/2 * *t2*1e-6 + c[ 3]/3 * *t3*1e-9 + c[ 4]/4 * *t4*1e-12 - c[ 5] * *recipT*1e3 ) // if low temp
                                                     (               c[13] + c[8] * temp*1e-3 + c[9]/2 * *t2*1e-6 + c[10]/3 * *t3*1e-9 + c[11]/4 * *t4*1e-12 - c[12] * *recipT*1e3 ); // else if high range

            for( ic = c.begin() + 1; ic != icend; ++ic)
              *ic *= 1e-3; // scale the coefficients again to keep units consistent
            dE0dT <<= dE0dT + *massFracs_[n] * cond( temp <= c[0], c[1] + c[2] * temp*1e-3 + c[ 3] * *t2*1e-6 + c[ 4] * *t3*1e-9 + c[ 5] * *recipRecipT*1e6 ) // if low temp
                                                   (               c[8] + c[9] * temp*1e-3 + c[10] * *t2*1e-6 + c[11] * *t3*1e-9 + c[12] * *recipRecipT*1e6 ); // else if high range
            break;
        }
      }
    else{
      /* else temperature can be out of bounds low, low temp, high temp, or out of bounds high
       * if out of bounds, properties are interpolated from min or max temp using a constant cv
       */
        switch (polyType) {
          case SIMPLE: // constant cv
            c[3] -= GasConstant; // change coefficients from enthalpy to internal energy
            c[1] -= GasConstant * c[0];
            delE0<<=delE0 - *massFracs_[n] * ( c[1] + c[3] * (temp-c[0]) )
                                             / molecularWeights[n];

            dE0dT<<=dE0dT + *massFracs_[n] * c[3]
                                               / molecularWeights[n];
            break;
          case NASA2:
            c[1] -= 1.0; //change coefficients from enthalpy to internal energy
            c[8] -= 1.0;
            for( ; ic != icend; ++ic)
              *ic *= GasConstant / molecularWeights[n];
            delE0 <<= delE0 - *massFracs_[n] * cond( temp <= c[0] && temp >= minT, c[ 6] + c[1] * temp + c[2]/2 * *t2 + c[ 3]/3 * *t3 + c[ 4]/4 * *t4 + c[ 5]/5 * *t5 ) // if low temp
                                                   ( temp >  c[0] && temp <= maxT, c[13] + c[8] * temp + c[9]/2 * *t2 + c[10]/3 * *t3 + c[11]/4 * *t4 + c[12]/5 * *t5 )  // else if high temp
                                                   ( temp < minT, c[1] * temp + c[2] * minT * (temp - minT/2) + c[ 3] * minT * minT * (temp - 2*minT/3) + c[ 4]*pow(minT,3) * (temp - 3*minT/4) + c[ 5]*pow(minT,4) * (temp - 4*minT/5) + c[ 6] ) // else if out of bounds - low
                                                   (              c[8] * temp + c[9] * maxT * (temp - maxT/2) + c[10] * maxT * maxT * (temp - 2*maxT/3) + c[11]*pow(maxT,3) * (temp - 3*maxT/4) + c[12]*pow(maxT,4) * (temp - 4*maxT/5) + c[13] ); // else out of bounds - high


            dE0dT <<= dE0dT + *massFracs_[n] * cond( temp <= c[0] && temp >= minT, c[1] + c[2] * temp + c[ 3] * *t2 + c[ 4] * *t3 + c[ 5] * *t4 ) // if low temp
                                                     ( temp >  c[0] && temp <= maxT, c[8] + c[9] * temp + c[10] * *t2 + c[11] * *t3 + c[12] * *t4 )  // else if high temp
                                                     ( temp < minT, c[1] + c[2] * minT + c[ 3] * minT * minT + c[ 4] * pow(minT,3) + c[ 5] * pow(minT,4) )  // else if out of bounds - low
                                                     (              c[8] + c[9] * maxT + c[10] * maxT * maxT + c[11] * pow(maxT,3) + c[12] * pow(maxT,4) ); // else out of bounds - high

            break;
          case SHOMATE2:
            double minTScaled = minT/1000;
            double maxTScaled = maxT/1000;
            c[1] -= GasConstant * 1e-3; //change coefficients from enthalpy to internal energy
            c[8] -= GasConstant * 1e-3;
            for( ; ic != icend; ++ic)
              *ic *= 1e6 / molecularWeights[n]; // scale the coefficients again to keep units consistent
            delE0 <<= delE0 - *massFracs_[n] * cond( temp <= c[0] && temp >= minT, c[ 6] + c[1] * temp*1e-3 + c[2]/2 * *t2*1e-6 + c[ 3]/3 * *t3*1e-9 + c[ 4]/4 * *t4*1e-12 - c[ 5] * *recipT*1e3 ) // if low temp
                                                   ( temp >  c[0] && temp <= maxT, c[13] + c[8] * temp*1e-3 + c[9]/2 * *t2*1e-6 + c[10]/3 * *t3*1e-9 + c[11]/4 * *t4*1e-12 - c[12] * *recipT*1e3 )  // else if high temp
                                                   ( temp < minT, c[1] * temp*1e-3 + c[2] * minTScaled * (temp*1e-3 - minTScaled/2) + c[ 3] * minTScaled * minTScaled * (temp*1e-3 - 2*minTScaled/3) + c[ 4]*pow(minTScaled,3) * (temp*1e-3 - 3*minTScaled/4) - c[ 5]*pow(minTScaled,-1) * (-temp*1e-3 / minTScaled + 2) + c[ 6] ) // else if out of bounds - low
                                                   (              c[8] * temp*1e-3 + c[9] * maxTScaled * (temp*1e-3 - maxTScaled/2) + c[10] * maxTScaled * maxTScaled * (temp*1e-3 - 2*maxTScaled/3) + c[11]*pow(maxTScaled,3) * (temp*1e-3 - 3*maxTScaled/4) - c[12]*pow(maxTScaled,-1) * (-temp*1e-3 / maxTScaled + 2) + c[13] ); // else out of bounds - high


            for( ic = c.begin() + 1; ic != icend; ++ic)
              *ic *= 1e-3; // scale the coefficients again to keep units consistent
            dE0dT <<= dE0dT + *massFracs_[n] * cond( temp <= c[0] && temp >= minT, c[1] + c[2] * temp*1e-3 + c[ 3] * *t2*1e-6 + c[ 4] * *t3*1e-9 + c[ 5] * *recipRecipT*1e6 ) // if low temp
                                                   ( temp >  c[0] && temp <= maxT, c[8] + c[9] * temp*1e-3 + c[10] * *t2*1e-6 + c[11] * *t3*1e-9 + c[12] * *recipRecipT*1e6 )  // else if high temp
                                                   ( temp < minT, c[1] + c[2] * minTScaled + c[ 3] * minTScaled * minTScaled + c[ 4] * pow(minTScaled,3) + c[ 5] * pow(minTScaled,-2) )  // else if out of bounds - low
                                                   (              c[8] + c[9] * maxTScaled + c[10] * maxTScaled * maxTScaled + c[11] * pow(maxTScaled,3) + c[12] * pow(maxTScaled,-2) ); // else out of bounds - high

        }
      }
    }
    // Newton's method to find root
    res <<= delE0/dE0dT;
    temp <<= temp + res;
    isConverged = nebo_max( abs(res) ) < 1e-3; // Converged when the temperature has changed by less than 1e-3
  }
}


//--------------------------------------------------------------------

template< typename FieldT >
TemperatureFromE0<FieldT>::
Builder::Builder( const Expr::Tag& resultTag,
                  const Expr::Tag& massFracTag,
                  const Expr::Tag& e0Tag,
                  const Expr::Tag& keTag )
: ExpressionBuilder( resultTag ),
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
