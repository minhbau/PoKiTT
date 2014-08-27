#ifndef DiffusionCoeff_Expr_h
#define DiffusionCoeff_Expr_h

//#define TIMINGS

#include <expression/Expression.h>
#ifdef TIMINGS
#include <boost/timer.hpp>
#endif

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
  Cantera::MixTransport* const transport_; //gas mixture to be evaluated
  const int nSpec_; //number of species to iterate over
  const double pressure_;
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
    mmwTag_( mmwTag ),
    transport_( dynamic_cast<Cantera::MixTransport*>( CanteraObjects::get_transport() )),
    nSpec_( transport_->thermo().nSpecies() ),
    pressure_( transport_->thermo().pressure() )
{
  this->set_gpu_runnable( true );

  massFracTags_.clear();
  for( size_t n=0; n<nSpec_; ++n ){
    std::ostringstream name;
    name << massFracTag.name() << "_" << n;
    massFracTags_.push_back( Expr::Tag(name.str(),massFracTag.context()) );
  }

  for( size_t n=0; n<nSpec_; ++n)
    indices_.push_back(std::vector<int>(nSpec_));

  size_t ij=0;
  for( size_t i=0; i<nSpec_; ++i ){
    for( size_t j=i; j<nSpec_; ++j, ++ij ){
      indices_[i][j]=ij;
      indices_[j][i]=ij;
    }
  }
}

//--------------------------------------------------------------------

template< typename FieldT >
DiffusionCoeff<FieldT>::
~DiffusionCoeff()
{
  CanteraObjects::restore_transport( transport_ );
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
#ifdef TIMINGS
boost::timer timer;
#endif

  using namespace SpatialOps;
  std::vector< FieldT* >& results = this->get_value_vec();

  const FieldT& p = *p_;

  //pre-compute powers of log(t) for polynomial evaluation of binary diffusion coefficients
  SpatFldPtr<FieldT> logtPtr   = SpatialFieldStore::get<FieldT>(*temperature_);
  SpatFldPtr<FieldT> logttPtr  = SpatialFieldStore::get<FieldT>(*temperature_);
  SpatFldPtr<FieldT> logtttPtr = SpatialFieldStore::get<FieldT>(*temperature_);
  SpatFldPtr<FieldT> tThreeHalvesPtr;
  SpatFldPtr<FieldT> logt4Ptr;

  if( transport_->model() == Cantera::cMixtureAveraged ) { // as opposed to CK mode
    logt4Ptr        = SpatialFieldStore::get<FieldT>(*temperature_);
    tThreeHalvesPtr = SpatialFieldStore::get<FieldT>(*temperature_);
  }

  FieldT& logt   = *logtPtr; // log(t)
  FieldT& logtt  = *logttPtr; // log(t)*log(t)
  FieldT& logttt = *logtttPtr; // log(t)*log(t)*log(t)
  FieldT& logt4  = *logt4Ptr; // pow( log(t), 4)
  FieldT& tThreeHalves = *tThreeHalvesPtr; // t^(3/2)

  logt   <<= log( *temperature_ );
  logtt  <<= logt * logt;
  logttt <<= logtt * logt;
  if( transport_->model() == Cantera::cMixtureAveraged ) {
    logt4        <<= logttt * logt;
    tThreeHalves <<= pow( *temperature_, 1.5 );
  }

  SpatFldPtr<FieldT> dPtr = SpatialFieldStore::get<FieldT>(*temperature_);
  FieldT& d = *dPtr;
  SpatFldPtr<FieldT> sum1Ptr = SpatialFieldStore::get<FieldT>(*temperature_);
  FieldT& sum1 = *sum1Ptr;
  SpatFldPtr<FieldT> sum2Ptr = SpatialFieldStore::get<FieldT>(*temperature_);
  FieldT& sum2 = *sum2Ptr;

  const std::vector< std::vector<double> >& diffusionPolyCoeffs = transport_->getDiffusionPolyCoefficients();

  const std::vector<double>& molecularWeights = transport_->thermo().molecularWeights();

  std::vector<double> molecularWeightsInv(nSpec_);
  for( size_t n=0; n<nSpec_; ++n )
    molecularWeightsInv[n] = 1 / molecularWeights[n];

//  for( size_t i=0; i<nSpec_; ++i){
//    d <<= 0.0;
//    for( size_t j=0; j<nSpec_; ++j){
//      if( j!=i){
//        const std::vector<double>& coeffs = diffusionPolyCoeffs[indices_[i][j]]; // coefficients for pair [i][j]
//        d <<= d + *massFracs_[j] * *mmw_ / molecularWeights[j] / ( tThreeHalves * (coeffs[0] + coeffs[1] * logt + coeffs[2] * logtt + coeffs[3] * logttt + coeffs[4] * logt4) ); // polynomial in t for binary diffusion coefficients
//      }
//    }
//    const std::vector<double>& coeffs = diffusionPolyCoeffs[indices_[i][i]];
//    *results[i] <<= cond( d >= 0.0, ( 1 - *massFracs_[i] ) / (pressure_ * d) ) // use mixing rule if temporary sum is non-negative
//                        (           tThreeHalves * (coeffs[0] + coeffs[1] * logt + coeffs[2] * logtt + coeffs[3] * logttt + coeffs[4] * logt4) / pressure_); // otherwise use self-diffusion coefficient
//  }

  for( size_t i=0; i<nSpec_; ++i){
    sum1 <<= 0.0;
    sum2 <<= 0.0;
    for( size_t j=0; j<nSpec_; ++j){
      if( j != i){
        const std::vector<double>& coeffs = diffusionPolyCoeffs[indices_[i][j]]; // coefficients for pair [i][j]
        d <<= *massFracs_[j] / ( tThreeHalves * (coeffs[0] + coeffs[1] * logt + coeffs[2] * logtt + coeffs[3] * logttt + coeffs[4] * logt4) ); // polynomial in t for binary diffusion coefficients
        sum1 <<= sum1 + d * molecularWeightsInv[j];
        sum2 <<= sum2 + d;
      }
    }
    const std::vector<double>& coeffs = diffusionPolyCoeffs[indices_[i][i]];
    *results[i] <<= 1 / ( p * *mmw_ * ( sum1 + sum2 * *massFracs_[i] / ( molecularWeights[i] - molecularWeights[i] * *massFracs_[i] ) ) );
  }

# ifdef TIMINGS
  std::cout<<"D time "<<timer.elapsed()<<std::endl;
# endif

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
