#ifndef DiffusionCoeff_Expr_h
#define DiffusionCoeff_Expr_h

#include <expression/Expression.h>

#include <pokitt/CanteraObjects.h> //include cantera wrapper

#include <cantera/transport.h>

/**
 *  \class  DiffusionCoeff
 *  \author Nate Yonkee
 *  \date July, 2014
 *
 *  \brief Calculates the mixture averaged diffusion coefficients [m^2/s].
 *
 *  This class calculates the diffusion coefficient
 *  appropriate for calculating the mass averaged diffusive flux with respect
 *  to the mass averaged velocity using gradients of the mass fraction.
 *  The diffusion coefficients are calculated using the following mixing rule,
 *
 * \f[
 * \frac{1}{D_{k,mix}'} = \sum^K_{j \neq k} \frac{x_j}{D_{kj}} + \frac{x_k}{1-y_{k}} \sum^K_{j \neq k} \frac{y_j}{D_{kj}}
 * \f]
 *
 * Where \f$ D_kj \f$ is the binary diffusion coefficient of species k and j,
 * \f$ y_j \f$ is the mass fraction of species j, and \f$ x_j \f$ is the mole
 * fraction of species j.
 *
 * Units are m^2/s.
 *
 */
template< typename FieldT >
class DiffusionCoeff
    : public Expr::Expression<FieldT>
{
  const Expr::Tag temperatureTag_;
  const Expr::Tag pTag_;
  Expr::TagList massFracTags_;
  const Expr::Tag mmwTag_;
  const FieldT* temperature_;
  const FieldT* p_;
  const FieldT* mmw_; // mixture molecular weight
  std::vector<const FieldT*> massFracs_;

  int nSpec_; //number of species to iterate over
  int modelType_; // type of model used by Cantera to estimate pure viscosity
  std::vector<double> molecularWeights_; // molecular weights
  std::vector<double> molecularWeightsInv_; // inverse of molecular weights (diving by MW is expensive)
  std::vector< std::vector<double> > binaryDCoefs_; // coefficients used by Cantera to calculate binary diffusino coefficients
  /* Cantera uses a polynomial in temperature to evaluate the binary diffusion coefficient of each pair [i][j] = [j][i]
   * indicies_[i][j] stores the index of the set of polynomial coefficients for the pair [i][j]
   * This is to simplify bookkeeping and does not affect the evaluation itself
   */
  std::vector< std::vector< int > > indices_; // keeps track of which polynomial coefficients correspond to which binary pair of species

  DiffusionCoeff( const Expr::Tag& temperatureTag,
                  const Expr::Tag& pTag,
                  const Expr::Tag& massFracTag,
                  const Expr::Tag& mmwTag );
public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a DiffusionCoeff expression
     *  @param resultTags the tags for the diffusion coefficient of each species, ordering is consistent with Cantera input file
     *  @param temperatureTag temperature
     *  @param pTag pressure
     *  @param massFracTag tag for mass fractions of each species, ordering is consistent with Cantera input
     *  @param mmwTag tag for mixture molecular weight
     */
    Builder( const Expr::TagList& resultTags,
             const Expr::Tag& temperatureTag,
             const Expr::Tag& pTag,
             const Expr::Tag& massFracTag,
             const Expr::Tag& mmwTag );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag temperatureTag_;
    const Expr::Tag pTag_;
    const Expr::Tag massFracTag_;
    const Expr::Tag mmwTag_;
  };

  ~DiffusionCoeff();
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
DiffusionCoeff<FieldT>::
DiffusionCoeff( const Expr::Tag& temperatureTag,
                const Expr::Tag& pTag,
                const Expr::Tag& massFracTag,
                const Expr::Tag& mmwTag )
  : Expr::Expression<FieldT>(),
    temperatureTag_( temperatureTag ),
    pTag_( pTag ),
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

  binaryDCoefs_ = trans->getDiffusionPolyCoefficients();
  modelType_ = trans->model();

  molecularWeights_ = trans->thermo().molecularWeights();
  molecularWeightsInv_.resize(nSpec_);
  for( size_t n=0; n<nSpec_; ++n)
    molecularWeightsInv_[n] = 1 / molecularWeights_[n];

  for( size_t n=0; n<nSpec_; ++n) // nSpec_ by nSpec_ vector of vectors
    indices_.push_back(std::vector<int>(nSpec_));

  size_t ij=0;
  for( size_t i=0; i<nSpec_; ++i ){
    for( size_t j=i; j<nSpec_; ++j, ++ij ){
      indices_[i][j]=ij;
      indices_[j][i]=ij;
    }
  }
  CanteraObjects::restore_transport( trans );
}

//--------------------------------------------------------------------

template< typename FieldT >
DiffusionCoeff<FieldT>::
~DiffusionCoeff()
{

}

//--------------------------------------------------------------------

template< typename FieldT >
void
DiffusionCoeff<FieldT>::
advertise_dependents( Expr::ExprDeps& exprDeps )
{
  exprDeps.requires_expression( temperatureTag_ );
  exprDeps.requires_expression( pTag_ );
  exprDeps.requires_expression( massFracTags_ );
  exprDeps.requires_expression( mmwTag_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
DiffusionCoeff<FieldT>::
bind_fields( const Expr::FieldManagerList& fml )
{
  const typename Expr::FieldMgrSelector<FieldT>::type& fm = fml.field_manager<FieldT>();
  temperature_ = &fm.field_ref( temperatureTag_ );
  p_ = &fm.field_ref( pTag_ );
  mmw_ = &fm.field_ref( mmwTag_ );
  massFracs_.clear();
  BOOST_FOREACH( const Expr::Tag& tag, massFracTags_ ){
    massFracs_.push_back( &fm.field_ref(tag) );
  }
}

//--------------------------------------------------------------------

template< typename FieldT >
void
DiffusionCoeff<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  std::vector< FieldT* >& mixD = this->get_value_vec();

  const FieldT& p = *p_;

  // pre-compute power of log(t) for the species viscosity polynomial
  SpatFldPtr<FieldT> tThreeHalvesPtr; // t^(3/2)
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

  SpatFldPtr<FieldT> dPtr    = SpatialFieldStore::get<FieldT>(*temperature_);
  SpatFldPtr<FieldT> sum1Ptr = SpatialFieldStore::get<FieldT>(*temperature_);
  SpatFldPtr<FieldT> sum2Ptr = SpatialFieldStore::get<FieldT>(*temperature_);

  FieldT& d    = *dPtr;
  FieldT& sum1 = *sum1Ptr;
  FieldT& sum2 = *sum2Ptr;

  if( modelType_ == Cantera::cMixtureAveraged ) { // as opposed to CK mode
    logt4Ptr        = SpatialFieldStore::get<FieldT>(*temperature_);
    tThreeHalvesPtr = SpatialFieldStore::get<FieldT>(*temperature_);

    *logt4Ptr        <<= logttt * logt;
    *tThreeHalvesPtr <<= pow( *temperature_, 1.5 );
  }

  for( size_t i=0; i<nSpec_; ++i){
    sum1 <<= 0.0;
    sum2 <<= 0.0;
    for( size_t j=0; j<nSpec_; ++j){
      if( j != i){
        const std::vector<double>& coefs = binaryDCoefs_[indices_[i][j]]; // coefficients for pair [i][j]
        if( modelType_ == Cantera::cMixtureAveraged )
          d <<= *massFracs_[j] / ( *tThreeHalvesPtr * ( coefs[0] + coefs[1] * logt + coefs[2] * logtt + coefs[3] * logttt + coefs[4] * *logt4Ptr ) ); // polynomial in t for binary diffusion coefficients
        else
          d <<= *massFracs_[j] / (                exp ( coefs[0] + coefs[1] * logt + coefs[2] * logtt + coefs[3] * logttt ) );
        sum1 <<= sum1 + d * molecularWeightsInv_[j];
        sum2 <<= sum2 + d;
      }
    }
    *mixD[i] <<= 1 / ( p * *mmw_ * ( sum1 + sum2 * *massFracs_[i] / ( molecularWeights_[i] - molecularWeights_[i] * *massFracs_[i] ) ) ); // mixing rule
  }

}

//--------------------------------------------------------------------

template< typename FieldT >
DiffusionCoeff<FieldT>::
Builder::Builder( const Expr::TagList& resultTags,
                  const Expr::Tag& temperatureTag,
                  const Expr::Tag& pTag,
                  const Expr::Tag& massFracTag,
                  const Expr::Tag& mmwTag )
: ExpressionBuilder( resultTags ),
  temperatureTag_( temperatureTag ),
  pTag_( pTag ),
  massFracTag_( massFracTag ),
  mmwTag_( mmwTag )
{}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
DiffusionCoeff<FieldT>::
Builder::build() const
{
  return new DiffusionCoeff<FieldT>( temperatureTag_, pTag_, massFracTag_, mmwTag_ );
}


#endif // DiffusionCoeff_Expr_h
