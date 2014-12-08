#ifndef ThermalConductivity_Expr_h
#define ThermalConductivity_Expr_h

#include <expression/Expression.h>

#include <pokitt/CanteraObjects.h> //include cantera wrapper

namespace pokitt{

/**
 *  \class  ThermalConductivity
 *  \author Nate Yonkee
 *  \date July, 2014
 *
 *  \brief Calculates the thermal conductivity using a mixing rule [W/m/K].
 *
 *  This class calculates the thermal conductivity using the following
 *  mixing rule,
 *
 * \f[
 * \lambda = 0.5 \left( \sum_n x_n \lambda_n
 * + \frac{1}{\sum_n x_n/\lambda_n}\right)
 * \f]
 *
 * Where \f$ \lambda_n \f$ is the thermal conductivity of pure n and
 * \f$ x_n \f$ is the mole fraction of species n.
 *
 * Units are W/m/K.
 *
 */

template< typename FieldT >
class ThermalConductivity
    : public Expr::Expression<FieldT>
{
  const Expr::Tag temperatureTag_;
  Expr::TagList massFracTags_;
  const Expr::Tag mmwTag_;
  const FieldT* temperature_;
  const FieldT* mmw_; // mixture molecular weight, needed to convert from mass to mole fractions
  std::vector<const FieldT*> massFracs_;

  int nSpec_; //number of species to iterate over
  std::vector< std::vector<double> > tCondCoefs_;
  int modelType_; // type of model used by Cantera to estimate pure viscosity
  std::vector<double> molecularWeights_; // molecular weights
  std::vector<double> molecularWeightsInv_; // inverse of molecular weights (diving by MW is expensive)

  ThermalConductivity( const Expr::Tag& temperatureTag,
                       const Expr::Tag& massFracTag,
                       const Expr::Tag& mmwTag );
public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a ThermalConductivity expression
     *  @param resultTag the tag for the mixture averaged thermal conductivity
     *  @param temperatureTag temperature
     *  @param massFracTag tag for mass fraction of each species, ordering must be consistent with Cantera input
     *  @param mmwTag tag for mixture molecular weight
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::Tag& temperatureTag,
             const Expr::Tag& massFracTag,
             const Expr::Tag& mmwTag);

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag temperatureTag_;
    const Expr::Tag massFracTag_;
    const Expr::Tag mmwTag_;
  };

  ~ThermalConductivity();
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
ThermalConductivity<FieldT>::
ThermalConductivity( const Expr::Tag& temperatureTag,
                     const Expr::Tag& massFracTag,
                     const Expr::Tag& mmwTag )
  : Expr::Expression<FieldT>(),
    temperatureTag_( temperatureTag ),
    mmwTag_( mmwTag )
{
  this->set_gpu_runnable( true );
  Cantera::MixTransport* trans = dynamic_cast<Cantera::MixTransport*>( CanteraObjects::get_transport() ); // cast gas transport object as mix transport
  nSpec_ = trans->thermo().nSpecies();

  massFracTags_.clear();
  for( size_t n=0; n<nSpec_; ++n ){
    std::ostringstream name;
    name << massFracTag.name() << "_" << n;
    massFracTags_.push_back( Expr::Tag(name.str(),massFracTag.context()) );
  }

  tCondCoefs_ = trans->getConductivityCoefficients();
  modelType_ = trans->model();

  molecularWeights_ = trans->thermo().molecularWeights();
  molecularWeightsInv_.resize(nSpec_);
  for( size_t n=0; n<nSpec_; ++n)
    molecularWeightsInv_[n] = 1 / molecularWeights_[n];

  CanteraObjects::restore_transport( trans );
}

//--------------------------------------------------------------------

template< typename FieldT >
ThermalConductivity<FieldT>::
~ThermalConductivity()
{

}

//--------------------------------------------------------------------

template< typename FieldT >
void
ThermalConductivity<FieldT>::
advertise_dependents( Expr::ExprDeps& exprDeps )
{
  exprDeps.requires_expression( temperatureTag_ );
  exprDeps.requires_expression( massFracTags_ );
  exprDeps.requires_expression( mmwTag_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
ThermalConductivity<FieldT>::
bind_fields( const Expr::FieldManagerList& fml )
{
  const typename Expr::FieldMgrSelector<FieldT>::type& fm = fml.field_manager<FieldT>();
  temperature_ = &fm.field_ref( temperatureTag_ );
  mmw_ = &fm.field_ref( mmwTag_ );
  massFracs_.clear();
  BOOST_FOREACH( const Expr::Tag& tag, massFracTags_ ){
    massFracs_.push_back( &fm.field_ref(tag) );
  }
}

//--------------------------------------------------------------------

template< typename FieldT >
void
ThermalConductivity<FieldT>::
evaluate()
{
  using namespace SpatialOps;

  FieldT& mixTCond = this->value();

  // pre-compute powers of temperature used in polynomial evaluations
  SpatFldPtr<FieldT> sqrtTPtr;
  SpatFldPtr<FieldT> logtPtr   = SpatialFieldStore::get<FieldT>(*temperature_); // log(t)
  SpatFldPtr<FieldT> logttPtr  = SpatialFieldStore::get<FieldT>(*temperature_); // log(t)*log(t)
  SpatFldPtr<FieldT> logtttPtr = SpatialFieldStore::get<FieldT>(*temperature_); // log(t)*log(t)*log(t)
  SpatFldPtr<FieldT> logt4Ptr; // log(t)*log(t)*log(t)*log(t)

  FieldT& logt   = *logtPtr;
  FieldT& logtt  = *logttPtr;
  FieldT& logttt = *logtttPtr;

  logt   <<= log( *temperature_ );
  logtt  <<= logt * logt;
  logttt <<= logtt * logt;

  SpatFldPtr<FieldT> speciesTCondPtr = SpatialFieldStore::get<FieldT>(*temperature_); // temporary to store the thermal conductivity for an individual species
  SpatFldPtr<FieldT> sumPtr          = SpatialFieldStore::get<FieldT>(*temperature_); // for mixing rule
  SpatFldPtr<FieldT> inverseSumPtr   = SpatialFieldStore::get<FieldT>(*temperature_); // 1/sum()

  FieldT& speciesTCond = *speciesTCondPtr;
  FieldT& sum          = *sumPtr;
  FieldT& inverseSum   = *inverseSumPtr;

  sum        <<= 0.0; // set sum to 0 before loop
  inverseSum <<= 0.0; // set inverse sum to 0 before loop

  if( modelType_ == Cantera::cMixtureAveraged ) { // as opposed to CK mode
    logt4Ptr = SpatialFieldStore::get<FieldT>(*temperature_);
    sqrtTPtr = SpatialFieldStore::get<FieldT>(*temperature_);

    *logt4Ptr <<= logttt * logt;
    *sqrtTPtr <<= sqrt( *temperature_ );
  }

  for( size_t n = 0; n < nSpec_; ++n){
    const std::vector<double>& tCondCoefs = tCondCoefs_[n];
    if( modelType_ == Cantera::cMixtureAveraged )
      speciesTCond <<= *sqrtTPtr * ( tCondCoefs[0] + tCondCoefs[1] * logt + tCondCoefs[2] * logtt + tCondCoefs[3] * logttt + tCondCoefs[4] * *logt4Ptr );
    else
      speciesTCond <<=         exp ( tCondCoefs[0] + tCondCoefs[1] * logt + tCondCoefs[2] * logtt + tCondCoefs[3] * logttt );
    sum        <<= sum        + *massFracs_[n] * speciesTCond * molecularWeightsInv_[n];
    inverseSum <<= inverseSum + *massFracs_[n] / speciesTCond * molecularWeightsInv_[n];
  }

  mixTCond <<= 0.5 * ( sum * *mmw_ + 1 / (inverseSum * *mmw_) ); // mixing rule
}

//--------------------------------------------------------------------

template< typename FieldT >
ThermalConductivity<FieldT>::
Builder::Builder( const Expr::Tag& resultTag,
                  const Expr::Tag& temperatureTag,
                  const Expr::Tag& massFracTag,
                  const Expr::Tag& mmwTag )
: ExpressionBuilder( resultTag ),
  temperatureTag_( temperatureTag ),
  massFracTag_( massFracTag ),
  mmwTag_( mmwTag )
{}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
ThermalConductivity<FieldT>::
Builder::build() const
{
  return new ThermalConductivity<FieldT>( temperatureTag_, massFracTag_, mmwTag_ );
}

} // namespace pokitt

#endif // ThermalConductivity_Expr_h
