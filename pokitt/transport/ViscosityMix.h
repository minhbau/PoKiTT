#ifndef Viscosity_Expr_h
#define Viscosity_Expr_h

#include <expression/Expression.h>

#include <pokitt/CanteraObjects.h> //include cantera wrapper

namespace pokitt{

// viscosity can be calculated using either n temporary fields or 2n temporary fields (default)
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
  const Expr::TagList massFracTags_;
  const FieldT* temperature_;
  std::vector<const FieldT*> massFracs_;

  int nSpec_; //number of species to iterate over
  std::vector< std::vector<double> > viscosityCoefs_; // Cantera's vector of coefficients for the pure viscosity polynomials
  int modelType_; // type of model used by Cantera to estimate pure viscosity
  std::vector<double> molecularWeights_; // molecular weights
  std::vector<double> molecularWeightsInv_; // inverse of molecular weights (diving by MW is expensive)
  std::vector< std::vector<double> > molecularWeightRatios_; // pre-compute ratios of molecular weights to reduce evaluation time
  std::vector< std::vector<double> > denominator_; // pre-compute the denominator of the mixing rule to reduce evaluation time

  Viscosity( const Expr::Tag& temperatureTag,
             const Expr::TagList& massFracTags );
public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a Viscosity expression
     *  @param resultTag tag for the mixture averaged viscosity
     *  @param temperatureTag temperature
     *  @param massFracTags tag for mass fraction of each species, ordering must be consistent with Cantera input
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::Tag& temperatureTag,
             const Expr::TagList& massFracTags );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag temperatureTag_;
    const Expr::TagList massFracTags_;
  };

  ~Viscosity(){}
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
             const double suthConstant,
             const double refVisc,
             const double refTemp );
    /**
     *  @brief Build a SutherlandViscosity expression
     *  @param result tag for the viscosity
     *  @param temperatureTag temperature
     */
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
           const Expr::TagList& massFracTags )
  : Expr::Expression<FieldT>(),
    temperatureTag_( temperatureTag ),
    massFracTags_( massFracTags )
{
  this->set_gpu_runnable( true );

  Cantera::MixTransport* trans = dynamic_cast<Cantera::MixTransport*>( CanteraObjects::get_transport() ); // cast gas transport object as mix transport
  nSpec_ = trans->thermo().nSpecies();

  viscosityCoefs_ = trans->getViscosityCoefficients();
  modelType_ = trans->model();

  molecularWeights_ = trans->thermo().molecularWeights();
  molecularWeightsInv_.resize(nSpec_);
  for( size_t n=0; n<nSpec_; ++n)
    molecularWeightsInv_[n] = 1 / molecularWeights_[n];

  for( size_t n=0; n!=nSpec_; ++n){
    molecularWeightRatios_.push_back( std::vector<double>(nSpec_,0.0) ); // nSpec_ by nSpec_ vector of vectors
    denominator_.push_back( std::vector<double>(nSpec_,0.0) ); // nSpec_ by nSpec_ vector of vectors
  }

# ifdef NFIELDS // evaluation requires nSpec_ temporary fields
  for(  size_t k=0; k!=nSpec_; ++k){
    for( size_t j=0; j!=nSpec_; ++j){
      molecularWeightRatios_[k][j] = pow( molecularWeights_[k]/molecularWeights_[j], 0.25); // pre-compute value to reduce evaluation time
      denominator_[j][k] = sqrt(8.0) * sqrt( 1 + molecularWeights_[k]/molecularWeights_[j] ) * molecularWeights_[j] / molecularWeights_[k]; // pre-compute value to reduce evaluation time
      denominator_[j][k] = 1/ denominator_[j][k];
    }
  }
# else // evaluation requires 2*nSpec_ temporary fields
  for(  size_t k=0; k!=nSpec_; ++k){
    for( size_t j=0; j!=nSpec_; ++j){
      denominator_[j][k] = sqrt(8.0) * sqrt( 1 + molecularWeights_[k] / molecularWeights_[j] ) * molecularWeights_[j];
      denominator_[j][k] = 1/ denominator_[j][k];
    }
  }
  for( size_t k=0; k!=nSpec_; ++k){
    for( size_t j=0; j!=k; ++j){
      molecularWeightRatios_[k][j] =      molecularWeights_[j] / molecularWeights_[k]       ;
      molecularWeightRatios_[j][k] = pow( molecularWeights_[j] / molecularWeights_[k], 0.25);
    }
  }
# endif
  CanteraObjects::restore_transport( trans );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
Viscosity<FieldT>::
advertise_dependents( Expr::ExprDeps& exprDeps )
{
  exprDeps.requires_expression( temperatureTag_ );
  exprDeps.requires_expression( massFracTags_   );
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
  using namespace SpatialOps;

  FieldT& mixVis = this->value();

  std::vector< SpatFldPtr<FieldT> > sqrtSpeciesVis;
  for( size_t n=0; n<nSpec_; ++n)
    sqrtSpeciesVis.push_back(SpatialFieldStore::get<FieldT>(*temperature_));

  // pre-compute power of log(t) for the species viscosity polynomial
  SpatFldPtr<FieldT> tOneFourthPtr; // t^(1/4) -- may be used later on
  SpatFldPtr<FieldT> logtPtr = SpatialFieldStore::get<FieldT>(*temperature_);

  FieldT& logt = *logtPtr;
  logt <<= log( *temperature_ );

  if( modelType_ == Cantera::cMixtureAveraged ) { // as opposed to CK mode
    tOneFourthPtr = SpatialFieldStore::get<FieldT>(*temperature_);
    *tOneFourthPtr <<= pow( *temperature_, 0.25 );
  }

  for( size_t n = 0; n<nSpec_; ++n ){
    const std::vector<double>& viscCoefs = viscosityCoefs_[n];
    if( modelType_ == Cantera::cMixtureAveraged )
      *sqrtSpeciesVis[n] <<= *tOneFourthPtr * ( viscCoefs[0] + logt * ( viscCoefs[1] + logt * ( viscCoefs[2] + logt * ( viscCoefs[3] + logt * viscCoefs[4] ))) );
    else
      *sqrtSpeciesVis[n] <<=      exp ( 0.5 * ( viscCoefs[0] + logt * ( viscCoefs[1] + logt * ( viscCoefs[2] + logt *   viscCoefs[3]                        )) ));
  }

  mixVis <<= 0.0; // set result to 0 before summing species contributions
# ifdef NFIELDS // requires nSpec_ temporary fields
  SpatFldPtr<FieldT> phiPtr  = SpatialFieldStore::get<FieldT>(*temperature_);
  FieldT& phi = *phiPtr;

  for( size_t k=0; k!=nSpec_; ++k){ // start looping over species contributions
    phi <<= *massFracs_[k]; // begin sum with the k=j case when phi = y_k
    for( size_t j=0; j!=nSpec_; ++j){
      if( j!=k )
        phi <<= phi + *massFracs_[j] * ( 1 + *sqrtSpeciesVis[k] / *sqrtSpeciesVis[j] * molecularWeightRatios_[j][k] ) * ( 1 + *sqrtSpeciesVis[k] / *sqrtSpeciesVis[j] * molecularWeightRatios_[j][k] ) * denominator_[j][k]; // mixing rule
    }
    mixVis <<= mixVis + *sqrtSpeciesVis[k] * *sqrtSpeciesVis[k] * *massFracs_[k] / phi; // mixing rule
  }
# else // requires 2*nSpec_ temporary fields
  SpatFldPtr<FieldT> temporary = SpatialFieldStore::get<FieldT>(*temperature_);
  std::vector< SpatFldPtr<FieldT> > phivec;

  for(size_t k=0; k<nSpec_; ++k){
    phivec.push_back(SpatialFieldStore::get<FieldT>(*temperature_));
    *phivec[k] <<= *massFracs_[k] * molecularWeightsInv_[k]; // begin sum with the k=j case, phi = y_k / M_k
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
    mixVis <<= mixVis + *massFracs_[k] * ( *sqrtSpeciesVis[k] * *sqrtSpeciesVis[k] ) / ( *phivec[k] * molecularWeights_[k] ); // mixing rule
  }
# endif

}

//--------------------------------------------------------------------

template< typename FieldT >
Viscosity<FieldT>::
Builder::Builder( const Expr::Tag& resultTag,
                  const Expr::Tag& temperatureTag,
                  const Expr::TagList& massFracTags )
: ExpressionBuilder( resultTag ),
  temperatureTag_( temperatureTag ),
  massFracTags_( massFracTags )
{}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
Viscosity<FieldT>::
Builder::build() const
{
  return new Viscosity<FieldT>( temperatureTag_, massFracTags_ );
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
    t0_( refTemp ),  // ref temp for air (K)
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

} // namespace pokitt

#endif // Viscosity_Expr_h
