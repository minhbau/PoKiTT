#ifndef TemperaturePowers_Expr_h
#define TemperaturePowers_Expr_h

#include <expression/Expression.h>

/**
 *  \class TemperaturePowers
 */
template< typename FieldT >
class TemperaturePowers
 : public Expr::Expression<FieldT>
{
  typedef std::vector<FieldT*> SpecT;
  const Expr::Tag tTag_;
  const FieldT* t_;

  /* declare operators associated with this expression here */

    TemperaturePowers( const Expr::Tag& tTag );
public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a TemperaturePowers expression
     *  @param resultTag the tag for the value that this expression computes
     */
    Builder( const Expr::TagList& resultTags,
             const Expr::Tag& tTag );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag tTag_;
  };

  ~TemperaturePowers();
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
TemperaturePowers<FieldT>::
TemperaturePowers( const Expr::Tag& tTag )
  : Expr::Expression<FieldT>(),
    tTag_( tTag )
{}

//--------------------------------------------------------------------

template< typename FieldT >
TemperaturePowers<FieldT>::
~TemperaturePowers()
{}

//--------------------------------------------------------------------

template< typename FieldT >
void
TemperaturePowers<FieldT>::
advertise_dependents( Expr::ExprDeps& exprDeps )
{
  exprDeps.requires_expression( tTag_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
TemperaturePowers<FieldT>::
bind_fields( const Expr::FieldManagerList& fml )
{

  t_ = &fml.template field_ref< FieldT >( tTag_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
TemperaturePowers<FieldT>::
evaluate()
{
  SpecT& powers = this->get_value_vec();

  *powers[0] <<= *t_ * *t_; // t^2
  *powers[1] <<= *t_ * *powers[0]; // t^3
  *powers[2] <<= *t_ * *powers[1]; // t^4
  *powers[3] <<= *t_ * *powers[2]; // t^5

  *powers[4] <<= 1 / *t_; // t^-1
  *powers[5] <<= 1 / *powers[0]; // t^-2

}

//--------------------------------------------------------------------

template< typename FieldT >
TemperaturePowers<FieldT>::
Builder::Builder( const Expr::TagList& resultTags,
                  const Expr::Tag& tTag )
  : ExpressionBuilder( resultTags ),
    tTag_( tTag )
{}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
TemperaturePowers<FieldT>::
Builder::build() const
{
  return new TemperaturePowers<FieldT>( tTag_ );
}


#endif // TemperaturePowers_Expr_h
