#ifndef MixtureMolWeight_Expr_h
#define MixtureMolWeight_Expr_h

#include <expression/Expression.h>

#include <pokitt/CanteraObjects.h> //include cantera wrapper

namespace pokitt{

/**
 *  \class MixtureMolWeight
 */
template< typename FieldT >
class MixtureMolWeight
 : public Expr::Expression<FieldT>
{
  const std::vector<double>& specMW_;
  std::vector<double> molecularWeightsInv_;
  const int nSpec_;

  DECLARE_VECTOR_OF_FIELDS( FieldT, massFracs_ )

  MixtureMolWeight( const Expr::TagList& massFracTags,
                    const std::vector<double>& specMw  );
public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a MixtureMolWeight expression
     *  @param resultTag the tag for mixture molecular weight
     *  @param massFracTags tag for mass fractions of each species, ordering is consistent with specMW
     *  @param specMW vector of species molecular weights
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::TagList& massFracTags );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::TagList massFracTags_;
    std::vector<double> specMW_;
  };

  ~MixtureMolWeight(){}
  void evaluate();
};



// ###################################################################
//
//                          Implementation
//
// ###################################################################



template< typename FieldT >
MixtureMolWeight<FieldT>::
MixtureMolWeight( const Expr::TagList& massFracTags,
                  const std::vector<double>& specMW )
  : Expr::Expression<FieldT>(),
    specMW_(specMW),
    nSpec_( specMW.size() )
{
  this->set_gpu_runnable( true );

  assert( massFracTags.size() == nSpec_ );

  molecularWeightsInv_.resize(nSpec_);
  for( size_t n=0; n<nSpec_; ++n)
    molecularWeightsInv_[n] = 1 / specMW_[n];

  this->template create_field_vector_request<FieldT>( massFracTags, massFracs_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
MixtureMolWeight<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  FieldT& mixMW = this->value();

  const FieldT& y0 = massFracs_[0]->field_ref();

  mixMW <<= y0 * molecularWeightsInv_[0];
  for( size_t n=1; n<nSpec_; ++n ){
    const FieldT& yi = massFracs_[n]->field_ref();
    mixMW <<= mixMW + yi * molecularWeightsInv_[n];
  }
  mixMW <<= 1.0 / mixMW;
}

//--------------------------------------------------------------------

template< typename FieldT >
MixtureMolWeight<FieldT>::
Builder::Builder( const Expr::Tag& resultTag,
                  const Expr::TagList& massFracTags )
  : ExpressionBuilder( resultTag ),
    massFracTags_( massFracTags )
{
  Cantera_CXX::IdealGasMix* const gasMix = CanteraObjects::get_gasmix();
  specMW_ = gasMix->molecularWeights();
  CanteraObjects::restore_gasmix(gasMix);
}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
MixtureMolWeight<FieldT>::
Builder::build() const
{
  return new MixtureMolWeight<FieldT>( massFracTags_, specMW_ );
}

} // namespace pokitt

#endif // MixtureMolWeight_Expr_h
