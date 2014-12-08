#ifndef Pressure_Expr_h
#define Pressure_Expr_h

#include <expression/Expression.h>

#include <pokitt/CanteraObjects.h> //include cantera wrapper

#include <cantera/kernel/ct_defs.h> // contains value of gas constant

namespace pokitt{

/**
 *  \class  Pressure
 *  \author Nathan Yonkee
 *  \date October, 2014
 *
 *  \brief Calculates the Pressure of a mixture using the
 *  ideal gas equation of state [Pa]
 *
 * \f[
 * P = \frac{\rho RT}{MW_{mix}}
 * \f]
 *
 * Where \f$ MW_{mix}\f$ is the mixture molecular weight
 *
 */

template< typename FieldT >
class Pressure
    : public Expr::Expression<FieldT>
{
  const Expr::Tag tTag_;
  const Expr::Tag rhoTag_;
  const Expr::Tag mmwTag_;
  const FieldT* t_;
  const FieldT* rho_;
  const FieldT* mmw_;

  Pressure( const Expr::Tag& tTag,
            const Expr::Tag& rhoTag,
            const Expr::Tag& mmwTag );
public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a Pressure expression
     *  @param resultTag tag for the density
     *  @param tTag tag for temperature
     *  @param rhoTag tag for density
     *  @param mmwTag tag for mixture molecular weight
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::Tag& tTag,
             const Expr::Tag& rhoTag,
             const Expr::Tag& mmwTag );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag tTag_;
    const Expr::Tag rhoTag_;
    const Expr::Tag mmwTag_;
  };

  ~Pressure();
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
Pressure<FieldT>::
Pressure( const Expr::Tag& tTag,
          const Expr::Tag& rhoTag,
          const Expr::Tag& mmwTag )
  : Expr::Expression<FieldT>(),
    tTag_( tTag ),
    rhoTag_( rhoTag ),
    mmwTag_( mmwTag )
{
  this->set_gpu_runnable( true );
}

//--------------------------------------------------------------------

template< typename FieldT >
Pressure<FieldT>::
~Pressure()
{

}

//--------------------------------------------------------------------

template< typename FieldT >
void
Pressure<FieldT>::
advertise_dependents( Expr::ExprDeps& exprDeps )
{
  exprDeps.requires_expression( tTag_ );
  exprDeps.requires_expression( rhoTag_ );
  exprDeps.requires_expression( mmwTag_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
Pressure<FieldT>::
bind_fields( const Expr::FieldManagerList& fml )
{
  const typename Expr::FieldMgrSelector<FieldT>::type& fm = fml.field_manager<FieldT>();

  t_ = &fm.field_ref( tTag_ );
  rho_ = &fm.field_ref( rhoTag_ );
  mmw_ = &fm.field_ref( mmwTag_ );

}

//--------------------------------------------------------------------

template< typename FieldT >
void
Pressure<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  FieldT& p = this->value();

  p <<= *rho_ * Cantera::GasConstant * *t_ / *mmw_;
}
//--------------------------------------------------------------------

template< typename FieldT >
Pressure<FieldT>::
Builder::Builder( const Expr::Tag& resultTag,
         const Expr::Tag& tTag,
         const Expr::Tag& rhoTag,
         const Expr::Tag& mmwTag )
: ExpressionBuilder( resultTag ),
  tTag_( tTag ),
  rhoTag_( rhoTag ),
  mmwTag_( mmwTag )
{}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
Pressure<FieldT>::
Builder::build() const
{
  return new Pressure<FieldT>( tTag_, rhoTag_, mmwTag_);
}

} // namespace pokitt

#endif // Pressure_Expr_h
