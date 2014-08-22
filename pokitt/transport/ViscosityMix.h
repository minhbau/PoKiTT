#ifndef Viscosity_Expr_h
#define Viscosity_Expr_h

//#define TIMINGS

#include <expression/Expression.h>
#ifdef TIMINGS
#include <boost/timer.hpp>
#endif

#include <pokitt/CanteraObjects.h> //include cantera wrapper

#include <cantera/transport.h>

//#define NFIELDS

/**
 *  \class  Viscosity
 *  \author Nate Yonkee
 *  \date July, 2014
 *
 *  \brief Calculates the viscosity using a mixing rule [kg/m/s].
 *
 *  This class calculates the viscosity using the Wilke
 *  mixing rule,
 *
 * \f[
 * \mu = \sum_k \frac{\mu_k x_k}{\sum_j \Phi_{k,j} x_j}.
 * \f]
 *
 * \f[
 * \Phi_{k,j} = \frac{\left[1
 * + \sqrt{\left(\frac{\mu_k}{\mu_j}\sqrt{\frac{M_j}{M_k}}\right)}\right]^2}
 * {\sqrt{8}\sqrt{1 + M_k/M_j}}
 * \f]
 *
 * Where \f$ \mu_k \f$ is the viscosity of pure species k and
 * \f$ M_k \f$ is the molecular weight of species k.
 *
 * Units are kg/m/s.
 *
 */

template< typename FieldT >
class Viscosity
    : public Expr::Expression<FieldT>
{
  const Expr::Tag temperatureTag_;
  Expr::TagList massFracTags_;
  const FieldT* temperature_;
  std::vector<const FieldT*> massFracs_;
  Cantera::MixTransport* const trans_; // transport for mixture to be evaluated
  const int nSpec_; //number of species to iterate over
  std::vector< std::vector<double > > molecularWeightRatios_; // pre-compute ratios of molecular weights to reduce evaluation time
  std::vector< std::vector<double > > denominator_; // pre-compute the denominator of the mixing rule to reduce evaluation time

  Viscosity ( const Expr::Tag& temperatureTag,
              const Expr::Tag& massFracTag );
public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a Viscosity expression
     *  @param resultTag tag for the mixture averaged viscosity
     *  @param temperatureTag temperature
     *  @param massFracTag tag for mass fraction of each species, ordering must be consistent with Cantera input
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::Tag& temperatureTag,
             const Expr::Tag& massFracTag );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag temperatureTag_;
    const Expr::Tag massFracTag_;
  };

  ~Viscosity();
  void advertise_dependents( Expr::ExprDeps& exprDeps );
  void bind_fields( const Expr::FieldManagerList& fml );
  void evaluate();
};

/**
 *  @class SutherlandViscosity
 *  @date  October, 2008
 *  @author James C. Sutherland
 *
 *  Calculates the viscosity using Sutherland's law,
 *  \f[
 *     \mu = \mu_0 \frac{T_0 + C}{T+C} \left( \frac{T}{T_0} \right)^{3/2}
 *  \f]
 *  where:
 *  <ul>
 *   <li> \f$C\f$ is a constant with units Kelvin
 *   <li> \f$\mu_0\f$ is the reference viscosity with units \f$\mathrm{Pa\cdot s}\f$,
 *   <li> \f$T_0\f$ is the reference temperature with units Kelvin.
 *  </ul>
 */
template<typename FieldT>
class SutherlandViscosity
  : public Expr::Expression<FieldT>
{
  const double c_, refVisc_, refTemp_;
  const Expr::Tag tTag_;
  const FieldT* temp_;

  SutherlandViscosity( const double c,
                       const double refVisc,
                       const double refTemp,
                       const Expr::Tag& temperatureTag );

  void evaluate();
  void advertise_dependents( Expr::ExprDeps& exprDeps );
  void bind_fields( const Expr::FieldManagerList& fml );

public:

  class Builder : public Expr::ExpressionBuilder
  {
    const double c_, mu0_, t0_;
    const Expr::Tag tTag_;
  public:
    /**
     *  @brief Build a SutherlandViscosity expression
     *  @param result tag for the viscosity
     *  @param temperatureTag temperature
     *  @param suthConstant Constant for use in the viscosity correlation (K)
     *  @param refVisc reference viscosity (Pa s)
     *  @param refTemp reference temperature (K)
     */
    Builder( const Expr::Tag& result,
             const Expr::Tag& temperatureTag,
             const double suthConstant, ///< Constant for use in the viscosity correlation (K)
             const double refVisc,      ///< reference viscosity (Pa s)
             const double refTemp );    ///< reference temperature (K)
    Builder( const Expr::Tag& result,
             const Expr::Tag& temperatureTag );

    Expr::ExpressionBase* build() const;
  };

};

// ###################################################################
//
//                          Implementation
//
// ###################################################################



template< typename FieldT >
Viscosity<FieldT>::
Viscosity( const Expr::Tag& temperatureTag,
    const Expr::Tag& massFracTag )
    : Expr::Expression<FieldT>(),
      temperatureTag_( temperatureTag ),
      trans_( dynamic_cast<Cantera::MixTransport*>( CanteraObjects::get_transport() )), // cast gas transport object as mix transport
      nSpec_( trans_->thermo().nSpecies() )
      {
  this->set_gpu_runnable( true );

  massFracTags_.clear();
  for( size_t n=0; n<nSpec_; ++n ){
    std::ostringstream name;
    name << massFracTag.name() << "_" << n;
    massFracTags_.push_back( Expr::Tag(name.str(),massFracTag.context()) );
  }

  for( size_t n=0; n!=nSpec_; ++n){
    molecularWeightRatios_.push_back( std::vector<double>(nSpec_,0.0) ); // nSpec_ by nSpec_ vector of vectors
    denominator_.push_back( std::vector<double>(nSpec_,0.0) );
  }
  const std::vector<double>& molecularWeights = trans_->thermo().molecularWeights();

#ifdef NFIELDS // evaluation requires nSpec_ temporary fields
  for(  size_t k=0; k!=nSpec_; ++k){
    for( size_t j=0; j!=nSpec_; ++j){
      molecularWeightRatios_[k][j] = pow( molecularWeights[k]/molecularWeights[j], 0.25); // pre-compute value to reduce evaluation time
      denominator_[j][k] = sqrt(8.0) * sqrt( 1 + molecularWeights[k]/molecularWeights[j] ) * molecularWeights[j] / molecularWeights[k]; // pre-compute value to reduce evaluation time
      denominator_[j][k] = 1/ denominator_[j][k];
    }
  }
#else // evaluation requires 2*nSpec_ temporary fields
  for(  size_t k=0; k!=nSpec_; ++k){
    for( size_t j=0; j!=nSpec_; ++j){
      denominator_[j][k] = sqrt(8.0) * sqrt( 1 + molecularWeights[k] / molecularWeights[j] ) * molecularWeights[j];
      denominator_[j][k] = 1/ denominator_[j][k];
    }
  }
  for( size_t k=0; k!=nSpec_; ++k){
    for( size_t j=0; j!=k; ++j){
      molecularWeightRatios_[k][j] =      molecularWeights[j] / molecularWeights[k]       ;
      molecularWeightRatios_[j][k] = pow( molecularWeights[j] / molecularWeights[k], 0.25);
    }
  }
#endif

      }

//--------------------------------------------------------------------

template< typename FieldT >
Viscosity<FieldT>::
~Viscosity()
{
  CanteraObjects::restore_transport( trans_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
Viscosity<FieldT>::
advertise_dependents( Expr::ExprDeps& exprDeps )
{
  exprDeps.requires_expression( temperatureTag_ );
  exprDeps.requires_expression( massFracTags_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
Viscosity<FieldT>::
bind_fields( const Expr::FieldManagerList& fml )
{
  const typename Expr::FieldMgrSelector<FieldT>::type& fm = fml.field_manager<FieldT>();
  temperature_ = &fm.field_ref( temperatureTag_ );
  massFracs_.clear();
  BOOST_FOREACH( const Expr::Tag& tag, massFracTags_ ){
    massFracs_.push_back( &fm.field_ref(tag) );
  }
}

//--------------------------------------------------------------------

template< typename FieldT >
void
Viscosity<FieldT>::
evaluate()
{
#ifdef TIMINGS
boost::timer timer;
#endif
  using namespace SpatialOps;

  FieldT& result = this->value();

  std::vector< SpatFldPtr<FieldT> > sqrtSpeciesVis;
  for( size_t n=0; n<nSpec_; ++n)
    sqrtSpeciesVis.push_back(SpatialFieldStore::get<FieldT>(*temperature_));

  // pre-compute power of log(t) for the species viscosity polynomial
  SpatFldPtr<FieldT> logtPtr  = SpatialFieldStore::get<FieldT>(*temperature_);
  SpatFldPtr<FieldT> logttPtr  = SpatialFieldStore::get<FieldT>(*temperature_);
  SpatFldPtr<FieldT> logtttPtr  = SpatialFieldStore::get<FieldT>(*temperature_);
  SpatFldPtr<FieldT> logt4Ptr;
  SpatFldPtr<FieldT> tOneFourthPtr;

  if( trans_->model() == Cantera::cMixtureAveraged ) { // as opposed to CK mode
    logt4Ptr  = SpatialFieldStore::get<FieldT>(*temperature_);
    tOneFourthPtr  = SpatialFieldStore::get<FieldT>(*temperature_);
  }

  FieldT& logt = *logtPtr;
  FieldT& logtt = *logttPtr;
  FieldT& logttt = *logtttPtr;
  FieldT& logt4 = *logt4Ptr;
  FieldT& tOneFourth = *tOneFourthPtr;

  logt <<= log( *temperature_ );
  logtt <<= logt * logt;
  logttt <<= logtt * logt;

  if( trans_->model() == Cantera::cMixtureAveraged ){
    tOneFourth <<= pow( *temperature_, 0.25 );
    logt4 <<= logttt * logt;
  }

  const std::vector<std::vector<double> >& viscosityCoefs = trans_->getViscosityCoefficients();

  const std::vector<double>& molecularWeights = trans_->thermo().molecularWeights();

  std::vector<double> molecularWeightsInv(nSpec_);
  for( size_t n=0; n<nSpec_; ++n)
    molecularWeightsInv[n] = 1 / molecularWeights[n];

  for( size_t n = 0; n<nSpec_; ++n ){
    if( trans_->model() == Cantera::cMixtureAveraged ) // Cantera uses a 5 coefficient polynomial in temperature
      *sqrtSpeciesVis[n] <<= tOneFourth * ( viscosityCoefs[n][0] + viscosityCoefs[n][1] * logt + viscosityCoefs[n][2] * logtt + viscosityCoefs[n][3] * logttt + viscosityCoefs[n][4] * logt4 );
    else // CK mode
      *sqrtSpeciesVis[n] <<= exp ( 0.5 * (viscosityCoefs[n][0] + viscosityCoefs[n][1] * logt + viscosityCoefs[n][2] * logtt + viscosityCoefs[n][3] * logttt) );
  }

  result <<= 0.0; // set result to 0 before summing species contributions

#ifdef NFIELDS // requires nSpec_ temporary fields
  SpatFldPtr<FieldT> phiPtr  = SpatialFieldStore::get<FieldT>(*temperature_);
  FieldT& phi = *phiPtr;

  for( size_t k=0; k!=nSpec_; ++k){ // start looping over species contributions
    phi <<= *massFracs_[k]; // begin sum with the k=j case when phi = y_k
    for( size_t j=0; j!=nSpec_; ++j){
      if(j!=k)
        phi <<= phi + *massFracs_[j] * ( 1 + *sqrtSpeciesVis[k] / *sqrtSpeciesVis[j] * molecularWeightRatios_[j][k] ) * ( 1 + *sqrtSpeciesVis[k] / *sqrtSpeciesVis[j] * molecularWeightRatios_[j][k] ) * denominator_[j][k]; // mixing rule
    }
    result <<= result + *sqrtSpeciesVis[k] * *sqrtSpeciesVis[k] * *massFracs_[k] / phi; // mixing rule
  }
#else // requires 2*nSpec_ temporary fields
  SpatFldPtr<FieldT> temporary = SpatialFieldStore::get<FieldT>(*temperature_);
  std::vector< SpatFldPtr<FieldT> > phivec;

  for(size_t k=0; k<nSpec_; ++k){
    phivec.push_back(SpatialFieldStore::get<FieldT>(*temperature_));
    *phivec[k] <<= *massFracs_[k] * molecularWeightsInv[k]; // begin sum with the k=j case, phi = y_k / M_k
  }

  for( size_t k=0; k!=nSpec_; ++k){
    for( size_t j=0; j!=k; ++j){
      *temporary <<= ( 1 + *sqrtSpeciesVis[k] / *sqrtSpeciesVis[j] * molecularWeightRatios_[j][k] ) * ( 1 + *sqrtSpeciesVis[k] / *sqrtSpeciesVis[j] * molecularWeightRatios_[j][k] ) * denominator_[j][k]; // mixing rule
      *phivec[k] <<= *phivec[k] + *massFracs_[j] * *temporary;
      /* phi[j][k] is proportional to phi[k][j]
       * This saves some evaluation time at the cost of doubling memory requirements
       */
      *phivec[j] <<= *phivec[j] + *massFracs_[k] * *temporary * *sqrtSpeciesVis[j] * *sqrtSpeciesVis[j] / ( *sqrtSpeciesVis[k] * *sqrtSpeciesVis[k] );
    }
  }

  for( size_t k=0; k!=nSpec_; ++k){
    result<<=result + *massFracs_[k] * ( *sqrtSpeciesVis[k] * *sqrtSpeciesVis[k] ) / ( *phivec[k] * molecularWeights[k] ); // mixing rule
  }
#endif
#ifdef TIMINGS
    std::cout<<"visc time "<<timer.elapsed()<<std::endl;
#endif
}

//--------------------------------------------------------------------

template< typename FieldT >
Viscosity<FieldT>::
Builder::Builder( const Expr::Tag& resultTag,
                  const Expr::Tag& temperatureTag,
                  const Expr::Tag& massFracTag )
: ExpressionBuilder( resultTag ),
  temperatureTag_( temperatureTag ),
  massFracTag_( massFracTag )
{}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
Viscosity<FieldT>::
Builder::build() const
{
  return new Viscosity<FieldT>( temperatureTag_, massFracTag_ );
}

//====================================================================



//--------------------------------------------------------------------

template< typename FieldT >
SutherlandViscosity<FieldT>::
SutherlandViscosity( const double c,
                     const double refVisc,
                     const double refTemp,
                     const Expr::Tag& temperatureTag )
  : Expr::Expression<FieldT>(),
    c_( c ),
    refVisc_( refVisc ),
    refTemp_( refTemp ),
    tTag_( temperatureTag )
{
  this->set_gpu_runnable(true);
}

//--------------------------------------------------------------------

template< typename FieldT >
void
SutherlandViscosity<FieldT>::
advertise_dependents( Expr::ExprDeps& exprDeps )
{
  exprDeps.requires_expression( tTag_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
SutherlandViscosity<FieldT>::
bind_fields( const Expr::FieldManagerList& fml )
{
  temp_ = &fml.field_manager<FieldT>().field_ref( tTag_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
SutherlandViscosity<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  FieldT& visc = this->value();
  visc <<= refVisc_ * ( refTemp_ + c_ ) / ( *temp_ + c_ ) * pow( *temp_ / refTemp_, 1.5 );
}

//--------------------------------------------------------------------

template< typename FieldT >
SutherlandViscosity<FieldT>::Builder::
Builder( const Expr::Tag& result,
         const Expr::Tag& temperatureTag )
  : ExpressionBuilder(result),
    c_  ( 120      ), // constant for air (K)
    mu0_( 18.27e-6 ), // ref visc for air (Pa s)
    t0_( 291.15    ), // ref temp for air (K)
    tTag_( temperatureTag )
{}

//--------------------------------------------------------------------

template< typename FieldT >
SutherlandViscosity<FieldT>::Builder::
Builder( const Expr::Tag& result,
         const Expr::Tag& temperatureTag,
         const double suthConstant,
         const double refVisc,
         const double refTemp )
  : ExpressionBuilder(result),
    c_( suthConstant ),
    mu0_( refVisc ), // ref visc for air (Pa s)
    t0_( refTemp ), // ref temp for air (K)
    tTag_( temperatureTag )
{}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
SutherlandViscosity<FieldT>::Builder::
build() const
{
  return new SutherlandViscosity<FieldT>(c_,mu0_,t0_,tTag_);
}

#endif // Viscosity_Expr_h
