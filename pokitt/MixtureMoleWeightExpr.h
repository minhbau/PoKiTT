#ifndef MixtureMolWeight_Expr_h
#define MixtureMolWeight_Expr_h

#include <expression/Expression.h>

/**
 *  \class MixtureMolWeight
 */
template< typename FieldT >
class MixtureMolWeight
 : public Expr::Expression<FieldT>
{
  const Expr::Tag yiTag_;
  const FieldT* yi_;
  const std::vector<double> specMW_;

    MixtureMolWeight( const Expr::Tag& yiTag,
                      const std::vector<double>& specMw  );
public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a MixtureMolWeight expression
     *  @param resultTag the tag for mixture molecular weight
     *  @param yiTag tag for mass fractions
     *  @param specMW vector of species molecular weights
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::Tag& yiTag,
             const std::vector<double>& specMW );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag yiTag_;
    const std::vector<double> specMW_;
  };

  ~MixtureMolWeight();
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
MixtureMolWeight<FieldT>::
MixtureMolWeight( const Expr::Tag& yiTag,
                  const std::vector<double>& specMW )
  : Expr::Expression<FieldT>(),
    yiTag_( yiTag ),
    specMW_(specMW)
{
  specTags_.clear();
  for( int i=0; i<nspec_; ++i ){
    std::ostringstream name;
    name << yiTag.name() << "_" << i;
    specTags_.push_back( Expr::Tag(name.str(),yiTag.context()) );
  }
}

//--------------------------------------------------------------------

template< typename FieldT >
MixtureMolWeight<FieldT>::
~MixtureMolWeight()
{}

//--------------------------------------------------------------------

template< typename FieldT >
void
MixtureMolWeight<FieldT>::
advertise_dependents( Expr::ExprDeps& exprDeps )
{
  for( Expr::TagList::const_iterator itag=specTags_.begin(); itag!=specTags_.end(); ++itag ){
    exprDeps.requires_expression( *itag );
  }
}

//--------------------------------------------------------------------

template< typename FieldT >
void
MixtureMolWeight<FieldT>::
bind_fields( const Expr::FieldManagerList& fml )
{
  const Expr::FieldMgrSelector<VolField>::type& fm = fml.field_manager<VolField>();

  yi_.clear();
  for( Expr::TagList::const_iterator itag=specTags_.begin(); itag!=specTags_.end(); ++itag ){
    yi_.push_back( &fm.field_ref( *itag ) );
  }
}

//--------------------------------------------------------------------

template< typename FieldT >
void
MixtureMolWeight<FieldT>::
bind_operators( const SpatialOps::OperatorDatabase& opDB )
{

}

//--------------------------------------------------------------------

template< typename FieldT >
void
MixtureMolWeight<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  FieldT& mixMW = this->value();

  result <<= 0.0;
    for( size_t n=0; n<nspec_; ++n ){
      result <<= result + *yi_[n] / specMW_[n];
    }
    result <<= 1.0 / result;
}

//--------------------------------------------------------------------

template< typename FieldT >
MixtureMolWeight<FieldT>::
Builder::Builder( const Expr::Tag& resultTag,
                  const Expr::Tag& yiTag,
                  const std::vector<double>& specMW )
  : ExpressionBuilder( resultTag ),
    yiTag_( yiTag ),
    specMW_( specMW )
{}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
MixtureMolWeight<FieldT>::
Builder::build() const
{
  return new MixtureMolWeight<FieldT>( yiTag_, specMW_ );
}


#endif // MixtureMolWeight_Expr_h
