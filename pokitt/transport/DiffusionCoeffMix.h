#ifndef DiffusionCoeff_Expr_h
#define DiffusionCoeff_Expr_h

#include <expression/Expression.h>

#include <pokitt/CanteraObjects.h> //include cantera wrapper

namespace pokitt{

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
 * \frac{1}{D_{k,mass}'} = \sum^K_{j \neq k} \frac{x_j}{D_{kj}} + \frac{x_k}{1-y_{k}} \sum^K_{j \neq k} \frac{y_j}{D_{kj}}
 * \f]
 *
 * Where \f$ D_{kj} \f$ is the binary diffusion coefficient of species k and j,
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
  DECLARE_FIELDS( FieldT, temperature_, p_, mmw_ )
  DECLARE_VECTOR_OF_FIELDS( FieldT, massFracs_ )

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
                  const Expr::TagList& massFracTags,
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
     *  @param massFracTags tag for mass fractions of each species, ordering is consistent with Cantera input
     *  @param mmwTag tag for mixture molecular weight
     *  @param nghost the number of ghost cells to compute in
     */
    Builder( const Expr::TagList& resultTags,
             const Expr::Tag& temperatureTag,
             const Expr::Tag& pTag,
             const Expr::TagList& massFracTags,
             const Expr::Tag& mmwTag,
             const int nghost = DEFAULT_NUMBER_OF_GHOSTS );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag temperatureTag_;
    const Expr::Tag pTag_;
    const Expr::TagList massFracTags_;
    const Expr::Tag mmwTag_;
  };

  ~DiffusionCoeff(){}
  void evaluate();
};

//--------------------------------------------------------------------
//--------------------------------------------------------------------

/**
*  \class  DiffusionCoeffMol
 *  \author Nate Yonkee
 *  \date December, 2014
 *
 *  \brief Calculates the mixture averaged diffusion coefficients [m^2/s] for a mass flux
 *  relative to a mass averaged velocity using the gradient of the mole fraction.
 *
 *  This class calculates the diffusion coefficient
 *  appropriate for calculating the mass averaged diffusive flux with respect
 *  to the mass averaged velocity using gradients of the mole fraction.
 *  The diffusion coefficients are calculated using the following mixing rule,
 *
 * \f[
 * D_{k,mol} = (1-y_k) \left( \sum_{j\neq k} \frac{y_k \bar{M}}{M_k D_{kj}} \right)^{-1}
 * \f]
 *
 * Where \f$ D_kj \f$ is the binary diffusion coefficient of species k and j,
 * \f$ y_k \f$ is the mass fraction of species k, \f$ \bar{M} \f$ is the mixture molecular weight,
 * and \f$ M_k \f$ is the molecular weight of species k.
 *
 * Units are m^2/s.
 *
 */
template< typename FieldT >
class DiffusionCoeffMol
    : public Expr::Expression<FieldT>
{
  DECLARE_FIELDS( FieldT, temperature_, p_, mmw_ )
  DECLARE_VECTOR_OF_FIELDS( FieldT, massFracs_ )

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

  DiffusionCoeffMol( const Expr::Tag& temperatureTag,
                     const Expr::Tag& pTag,
                     const Expr::TagList& massFracTags,
                     const Expr::Tag& mmwTag );
public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a DiffusionCoeffMol expression
     *  @param resultTags the tags for the diffusion coefficient of each species, ordering is consistent with Cantera input file
     *  @param temperatureTag temperature
     *  @param pTag pressure
     *  @param massFracTags tag for mass fractions of each species, ordering is consistent with Cantera input
     *  @param mmwTag tag for mixture molecular weight
     *  @param nghost the number of ghost cells to compute in
     */
    Builder( const Expr::TagList& resultTags,
             const Expr::Tag& temperatureTag,
             const Expr::Tag& pTag,
             const Expr::TagList& massFracTags,
             const Expr::Tag& mmwTag,
             const int nghost = DEFAULT_NUMBER_OF_GHOSTS );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag temperatureTag_;
    const Expr::Tag pTag_;
    const Expr::TagList massFracTags_;
    const Expr::Tag mmwTag_;
  };

  ~DiffusionCoeffMol(){}
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
                const Expr::TagList& massFracTags,
                const Expr::Tag& mmwTag )
  : Expr::Expression<FieldT>()
{
  this->set_gpu_runnable( true );

  temperature_ = this->template create_field_request<FieldT>( temperatureTag );
  p_           = this->template create_field_request<FieldT>( pTag           );
  mmw_         = this->template create_field_request<FieldT>( mmwTag         );

  this->template create_field_vector_request<FieldT>( massFracTags, massFracs_ );


  Cantera::MixTransport* trans = dynamic_cast<Cantera::MixTransport*>( CanteraObjects::get_transport() ); // cast gas transport object as mix transport
  nSpec_ = trans->thermo().nSpecies();

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
void
DiffusionCoeff<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  std::vector< FieldT* >& mixD = this->get_value_vec();

  const FieldT& p    = p_          ->field_ref();
  const FieldT& temp = temperature_->field_ref();
  const FieldT& mmw  = mmw_        ->field_ref();

  // pre-compute power of log(t) for the species viscosity polynomial
  SpatFldPtr<FieldT> tThreeHalvesPtr; // t^(3/2)
  SpatFldPtr<FieldT> logtPtr   = SpatialFieldStore::get<FieldT>(temp); // log(t)

  FieldT& logt   = *logtPtr;

  logt   <<= log( temp );

  SpatFldPtr<FieldT> dPtr    = SpatialFieldStore::get<FieldT>(temp);
  SpatFldPtr<FieldT> sum1Ptr = SpatialFieldStore::get<FieldT>(temp);
  SpatFldPtr<FieldT> sum2Ptr = SpatialFieldStore::get<FieldT>(temp);

  FieldT& d    = *dPtr;
  FieldT& sum1 = *sum1Ptr;
  FieldT& sum2 = *sum2Ptr;

  if( modelType_ == Cantera::cMixtureAveraged ) { // as opposed to CK mode
    tThreeHalvesPtr = SpatialFieldStore::get<FieldT>(temp);
    *tThreeHalvesPtr <<= pow( temp, 1.5 );
  }

  for( size_t i=0; i<nSpec_; ++i){
    const FieldT& yi = massFracs_[i]->field_ref();
    sum1 <<= 0.0;
    sum2 <<= 0.0;
    for( size_t j=0; j<nSpec_; ++j){
      const FieldT& yj = massFracs_[j]->field_ref();
      if( j != i){
        const std::vector<double>& coefs = binaryDCoefs_[indices_[i][j]]; // coefficients for pair [i][j]
        if( modelType_ == Cantera::cMixtureAveraged )
          d <<= yj / ( *tThreeHalvesPtr * ( coefs[0] + logt * ( coefs[1] + logt * ( coefs[2] + logt * ( coefs[3] + logt * coefs[4] ))) )); // polynomial in t for binary diffusion coefficients
        else
          d <<= yj / (                exp ( coefs[0] + logt * ( coefs[1] + logt * ( coefs[2] + logt *   coefs[3]                    )) ));
        sum1 <<= sum1 + d * molecularWeightsInv_[j];
        sum2 <<= sum2 + d;
      }
    }
    *mixD[i] <<= 1 / ( p * mmw * ( sum1 + sum2 * yi / ( molecularWeights_[i] - molecularWeights_[i] * yi ) ) ); // mixing rule
  }

}

//--------------------------------------------------------------------

template< typename FieldT >
DiffusionCoeff<FieldT>::
Builder::Builder( const Expr::TagList& resultTags,
                  const Expr::Tag& temperatureTag,
                  const Expr::Tag& pTag,
                  const Expr::TagList& massFracTags,
                  const Expr::Tag& mmwTag,
                  const int nghost )
: ExpressionBuilder( resultTags, nghost ),
  temperatureTag_( temperatureTag ),
  pTag_( pTag ),
  massFracTags_( massFracTags ),
  mmwTag_( mmwTag )
{}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
DiffusionCoeff<FieldT>::
Builder::build() const
{
  return new DiffusionCoeff<FieldT>( temperatureTag_, pTag_, massFracTags_, mmwTag_ );
}


// ###################################################################
//
//                          Implementation
//
// ###################################################################



template< typename FieldT >
DiffusionCoeffMol<FieldT>::
DiffusionCoeffMol( const Expr::Tag& temperatureTag,
                   const Expr::Tag& pTag,
                   const Expr::TagList& massFracTags,
                   const Expr::Tag& mmwTag )
  : Expr::Expression<FieldT>()
{
  this->set_gpu_runnable( true );

  temperature_ = this->template create_field_request<FieldT>( temperatureTag );
  p_           = this->template create_field_request<FieldT>( pTag           );
  mmw_         = this->template create_field_request<FieldT>( mmwTag         );

  this->template create_field_vector_request<FieldT>( massFracTags, massFracs_ );


  Cantera::MixTransport* trans = dynamic_cast<Cantera::MixTransport*>( CanteraObjects::get_transport() ); // cast gas transport object as mix transport
  nSpec_ = trans->thermo().nSpecies();

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
void
DiffusionCoeffMol<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  std::vector< FieldT* >& mixD = this->get_value_vec();

  const FieldT& p    = p_          ->field_ref();
  const FieldT& temp = temperature_->field_ref();
  const FieldT& mmw  = mmw_        ->field_ref();

  // pre-compute power of log(t) for the species viscosity polynomial
  SpatFldPtr<FieldT> tThreeHalvesPtr; // t^(3/2)
  SpatFldPtr<FieldT> logtPtr   = SpatialFieldStore::get<FieldT>(temp); // log(t)

  FieldT& logt   = *logtPtr;

  logt   <<= log( temp );

  SpatFldPtr<FieldT> dInvPtr = SpatialFieldStore::get<FieldT>(temp);
  SpatFldPtr<FieldT> sum1Ptr = SpatialFieldStore::get<FieldT>(temp);

  FieldT& dInv    = *dInvPtr;
  FieldT& sum1 = *sum1Ptr;

  if( modelType_ == Cantera::cMixtureAveraged ) { // as opposed to CK mode
    tThreeHalvesPtr = SpatialFieldStore::get<FieldT>(temp);
    *tThreeHalvesPtr <<= pow( temp, 1.5 );
  }

  for( size_t i=0; i<nSpec_; ++i){
    const FieldT& yi = massFracs_[i]->field_ref();
    sum1 <<= 0.0;
    for( size_t j=0; j<nSpec_; ++j){
    const FieldT& yj = massFracs_[j]->field_ref();
      if( j != i){
        const std::vector<double>& coefs = binaryDCoefs_[indices_[i][j]]; // coefficients for pair [i][j]
        if( modelType_ == Cantera::cMixtureAveraged )
          dInv <<= yj / ( *tThreeHalvesPtr * ( coefs[0] + logt * ( coefs[1] + logt * ( coefs[2] + logt * ( coefs[3] + logt * coefs[4] ))) )); // polynomial in t for binary diffusion coefficients
        else
          dInv <<= yj / (                exp ( coefs[0] + logt * ( coefs[1] + logt * ( coefs[2] + logt *   coefs[3]                    )) ));
        sum1 <<= sum1 + dInv * molecularWeightsInv_[j];
      }
    }
    *mixD[i] <<= ( 1 - yi ) / ( p * sum1 * mmw ); // mixing rule
  }

}

//--------------------------------------------------------------------

template< typename FieldT >
DiffusionCoeffMol<FieldT>::
Builder::Builder( const Expr::TagList& resultTags,
                  const Expr::Tag& temperatureTag,
                  const Expr::Tag& pTag,
                  const Expr::TagList& massFracTags,
                  const Expr::Tag& mmwTag,
                  const int nghost )
: ExpressionBuilder( resultTags, nghost ),
  temperatureTag_( temperatureTag ),
  pTag_( pTag ),
  massFracTags_( massFracTags ),
  mmwTag_( mmwTag )
{}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
DiffusionCoeffMol<FieldT>::
Builder::build() const
{
  return new DiffusionCoeffMol<FieldT>( temperatureTag_, pTag_, massFracTags_, mmwTag_ );
}

} // namespace pokitt

#endif // DiffusionCoeff_Expr_h
