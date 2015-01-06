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
  const Expr::TagList massFracTags_;
  const std::vector<double>& specMW_;
  std::vector<double> molecularWeightsInv_;
  const int nSpec_;
  std::vector<const FieldT*> massFracs_;

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
MixtureMolWeight<FieldT>::
MixtureMolWeight( const Expr::TagList& massFracTags,
                  const std::vector<double>& specMW )
  : Expr::Expression<FieldT>(),
    massFracTags_( massFracTags ),
    specMW_(specMW),
    nSpec_( specMW.size() )
{
  assert( massFracTags.size() == nSpec_ );
  molecularWeightsInv_.resize(nSpec_);
  for( size_t n=0; n<nSpec_; ++n)
    molecularWeightsInv_[n] = 1 / specMW_[n];
}

//--------------------------------------------------------------------

template< typename FieldT >
void
MixtureMolWeight<FieldT>::
advertise_dependents( Expr::ExprDeps& exprDeps )
{
  exprDeps.requires_expression( massFracTags_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
MixtureMolWeight<FieldT>::
bind_fields( const Expr::FieldManagerList& fml )
{
  const typename Expr::FieldMgrSelector<FieldT>::type& fm = fml.field_manager<FieldT>();
  massFracs_.clear();
  BOOST_FOREACH( const Expr::Tag& tag, massFracTags_ ){
    massFracs_.push_back( &fm.field_ref(tag) );
  }
}

//--------------------------------------------------------------------

template< typename FieldT >
void
MixtureMolWeight<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  FieldT& mixMW = this->value();

  mixMW <<= *massFracs_[0] * molecularWeightsInv_[0];
  for( size_t n=1; n<nSpec_; ++n ){
    mixMW <<= mixMW + *massFracs_[n] * molecularWeightsInv_[n];
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
