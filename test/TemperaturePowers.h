#ifndef TemperaturePowers_Expr_h
#define TemperaturePowers_Expr_h

#include <expression/Expression.h>

#include <pokitt/CanteraObjects.h> //include cantera wrapper

#include <cantera/kernel/speciesThermoTypes.h> // contains definitions for which polynomial is being used
#include <cantera/IdealGasMix.h>

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

  TemperaturePowers( const Expr::Tag& tTag );

public:

  static const Expr::TagList& temperature_powers_tags();

  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a TemperaturePowers expression
     *  @param tTag the temperature
     */
    Builder( const Expr::Tag& tTag );

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
const Expr::TagList&
TemperaturePowers<FieldT>::temperature_powers_tags()
{
  using namespace Expr;
  static TagList tags = tag_list( Tag("T^2"  ,STATE_NONE),
                                  Tag("T^3"  ,STATE_NONE),
                                  Tag("T^4"  ,STATE_NONE),
                                  Tag("T^5"  ,STATE_NONE),
                                  Tag("1/T"  ,STATE_NONE),
                                  Tag("1/T^2",STATE_NONE) );
  return tags;
}


template< typename FieldT >
TemperaturePowers<FieldT>::
TemperaturePowers( const Expr::Tag& tTag )
  : Expr::Expression<FieldT>(),
    tTag_( tTag )
{
  this->set_gpu_runnable( true );
}

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
Builder::Builder( const Expr::Tag& tTag )
  : ExpressionBuilder( TemperaturePowers<FieldT>::temperature_powers_tags() ),
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
