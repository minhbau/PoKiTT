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
  const std::vector<double> specMW_;
  const int nSpec_;

  Expr::TagList massFracTags_;
  std::vector<const FieldT*> massFracs_;

  MixtureMolWeight( const Expr::Tag& massFracTag,
                    const std::vector<double>& specMw  );
public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a MixtureMolWeight expression
     *  @param resultTag the tag for mixture molecular weight
     *  @param massFracTag tag for mass fractions of each species, ordering is consistent with specMW
     *  @param specMW vector of species molecular weights
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::Tag& massFracTag,
             const std::vector<double>& specMW );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag massFracTag_;
    const std::vector<double> specMW_;
  };

  ~MixtureMolWeight();
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
MixtureMolWeight( const Expr::Tag& massFracTag,
                  const std::vector<double>& specMW )
  : Expr::Expression<FieldT>(),
    specMW_(specMW),
    nSpec_( specMW.size() )
{
  massFracTags_.clear();
  for( size_t n=0; n<nSpec_; ++n ){
    std::ostringstream name;
    name << massFracTag.name() << "_" << n;
    massFracTags_.push_back( Expr::Tag(name.str(),massFracTag.context()) );
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

  mixMW <<= *massFracs_[0] / specMW_[0];
  for( size_t n=1; n<nSpec_; ++n ){
    mixMW <<= mixMW + *massFracs_[n] / specMW_[n];
  }
  mixMW <<= 1.0 / mixMW;
}

//--------------------------------------------------------------------

template< typename FieldT >
MixtureMolWeight<FieldT>::
Builder::Builder( const Expr::Tag& resultTag,
                  const Expr::Tag& massFracTag,
                  const std::vector<double>& specMW )
  : ExpressionBuilder( resultTag ),
    massFracTag_( massFracTag ),
    specMW_( specMW )
{}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
MixtureMolWeight<FieldT>::
Builder::build() const
{
  return new MixtureMolWeight<FieldT>( massFracTag_, specMW_ );
}


#endif // MixtureMolWeight_Expr_h
