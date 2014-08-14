#ifndef ThermalConductivity_Expr_h
#define ThermalConductivity_Expr_h

//#define TIMINGS

#include <expression/Expression.h>
#ifdef TIMINGS
#include <boost/timer.hpp>
#endif

#include <pokitt/CanteraObjects.h> //include cantera wrapper

#include <cantera/transport.h>

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
  Cantera::MixTransport* const trans_; // transport for mixture to be evaluated
  const int nSpec_; //number of species to iterate over

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
  void bind_operators( const SpatialOps::OperatorDatabase& opDB );
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
      mmwTag_( mmwTag ),
      trans_( dynamic_cast<Cantera::MixTransport*>( CanteraObjects::self().get_transport() )),
      nSpec_( trans_->thermo().nSpecies() )
      {
  this->set_gpu_runnable( true );

  massFracTags_.clear();
  for( size_t n=0; n<nSpec_; ++n ){
    std::ostringstream name;
    name << massFracTag.name() << "_" << n;
    massFracTags_.push_back( Expr::Tag(name.str(),massFracTag.context()) );
  }
      }

//--------------------------------------------------------------------

template< typename FieldT >
ThermalConductivity<FieldT>::
~ThermalConductivity()
{
  CanteraObjects::self().restore_transport( trans_ );
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
  for (size_t n=0; n<nSpec_; ++n) massFracs_.push_back(&fm.field_ref( massFracTags_[n] ));
}

//--------------------------------------------------------------------

template< typename FieldT >
void
ThermalConductivity<FieldT>::
bind_operators( const SpatialOps::OperatorDatabase& opDB )
{

}

//--------------------------------------------------------------------

template< typename FieldT >
void
ThermalConductivity<FieldT>::
evaluate()
{
#ifdef TIMINGS
boost::timer timer;
#endif
  using namespace SpatialOps;

  FieldT& result = this->value();

  int mtype;
# ifdef ENABLE_CUDA
  mtype = GPU_INDEX;
# else
  mtype = CPU_INDEX;
# endif

  // pre-compute powers of temperature used in polynomial evaluations
  SpatFldPtr<FieldT> logtPtr  = SpatialFieldStore::get<FieldT>(*temperature_,mtype); // log(t)
  SpatFldPtr<FieldT> logttPtr  = SpatialFieldStore::get<FieldT>(*temperature_,mtype); // log(t)*log(t)
  SpatFldPtr<FieldT> logtttPtr  = SpatialFieldStore::get<FieldT>(*temperature_,mtype); // log(t)*log(t)*log(t)
  SpatFldPtr<FieldT> sqrtTPtr; // sqrt(t)
  SpatFldPtr<FieldT> logt4Ptr; // pow( log(t),4 )
  if( trans_->model() == Cantera::cMixtureAveraged ) { // as opposed to CK mode
    logt4Ptr  = SpatialFieldStore::get<FieldT>(*temperature_,mtype);
    sqrtTPtr  = SpatialFieldStore::get<FieldT>(*temperature_,mtype);
  }

  SpatFldPtr<FieldT> speciesTCondPtr  = SpatialFieldStore::get<FieldT>(*temperature_,mtype); // temporary to store the thermal conductivity for an individual species

  SpatFldPtr<FieldT> sumPtr  = SpatialFieldStore::get<FieldT>(*temperature_,mtype); // for mixing rule
  SpatFldPtr<FieldT> inverseSumPtr  = SpatialFieldStore::get<FieldT>(*temperature_,mtype);

  FieldT& speciesTCond = *speciesTCondPtr;

  FieldT& logt = *logtPtr;
  FieldT& logtt = *logttPtr;
  FieldT& logttt = *logtttPtr;
  FieldT& logt4 = *logt4Ptr;
  FieldT& sqrtT = *sqrtTPtr;

  FieldT& sum = *sumPtr;
  FieldT& inverseSum = *inverseSumPtr;

  logt <<= log( *temperature_ );
  logtt <<= logt * logt;
  logttt <<= logtt * logt;
  if( trans_->model() == Cantera::cMixtureAveraged ) {
    logt4 <<= logttt * logt;
    sqrtT <<= sqrt( *temperature_ );
  }

  const std::vector< std::vector<double> >& conductivityCoefs = trans_->getConductivityCoefficients();

  const std::vector<double>& molecularWeights = trans_->thermo().molecularWeights();

  std::vector<double> molecularWeightsInv(nSpec_);
  for( size_t n=0; n<nSpec_; ++n)
    molecularWeightsInv[n] = 1 / molecularWeights[n];

  sum <<= 0.0; // set sum to 0 before loop
  inverseSum <<= 0.0; // set inverse sum to 0 before loop
  for( size_t n = 0; n < nSpec_; ++n){
    if( trans_->model() == Cantera::cMixtureAveraged ) // Cantera uses a 5 coefficient polynomial in temperature
      speciesTCond <<= sqrtT * ( conductivityCoefs[n][0] + conductivityCoefs[n][1] * logt + conductivityCoefs[n][2] * logtt + conductivityCoefs[n][3] * logttt + conductivityCoefs[n][4] * logt4 );
    else // CK mode
      speciesTCond <<= exp ( conductivityCoefs[n][0] + conductivityCoefs[n][1] * logt + conductivityCoefs[n][2] * logtt + conductivityCoefs[n][3] * logttt );
    sum <<= sum + *massFracs_[n] * speciesTCond * molecularWeightsInv[n];
    inverseSum <<= inverseSum + *massFracs_[n] / speciesTCond * molecularWeightsInv[n];
  }

  result <<= 0.5 * ( sum * *mmw_ + 1 / (inverseSum * *mmw_) ); // mixing rule
#ifdef TIMINGS
    std::cout<<"tc time "<<timer.elapsed()<<std::endl;
#endif
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


#endif // ThermalConductivity_Expr_h
