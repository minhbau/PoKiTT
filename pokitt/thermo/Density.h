#ifndef Density_Expr_h
#define Density_Expr_h

#include <expression/Expression.h>

#include <pokitt/CanteraObjects.h> //include cantera wrapper

#include <cantera/kernel/ct_defs.h> // contains value of gas constant

namespace pokitt{

/**
 *  \class  Density
 *  \author Nathan Yonkee
 *  \date October, 2014
 *
 *  \brief Calculates the density of a mixture using the
 *  ideal gas equation of state [kg/m^3]
 *
 * \f[
 * \rho = \frac{P MW_{mix}}{RT}
 * \f]
 *
 * Where \f$ MW_{mix}\f$ is the mixture molecular weight
 *
 */

template< typename FieldT >
class Density
    : public Expr::Expression<FieldT>
{
  const Expr::Tag tTag_;
  const Expr::Tag pTag_;
  const Expr::Tag mmwTag_;
  const FieldT* t_;
  const FieldT* p_;
  const FieldT* mmw_;

  Density( const Expr::Tag& tTag,
           const Expr::Tag& pTag,
           const Expr::Tag& mmwTag );
public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a Density expression
     *  @param resultTag tag for the density
     *  @param tTag tag for temperature
     *  @param pTag tag for pressure
     *  @param mmwTag tag for mixture molecular weight
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::Tag& tTag,
             const Expr::Tag& pTag,
             const Expr::Tag& mmwTag );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag tTag_;
    const Expr::Tag pTag_;
    const Expr::Tag mmwTag_;
  };

  ~Density();
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
Density<FieldT>::
Density( const Expr::Tag& tTag,
         const Expr::Tag& pTag,
         const Expr::Tag& mmwTag )
  : Expr::Expression<FieldT>(),
    tTag_( tTag ),
    pTag_( pTag ),
    mmwTag_( mmwTag )
{
  this->set_gpu_runnable( true );
}

//--------------------------------------------------------------------

template< typename FieldT >
Density<FieldT>::
~Density()
{

}

//--------------------------------------------------------------------

template< typename FieldT >
void
Density<FieldT>::
advertise_dependents( Expr::ExprDeps& exprDeps )
{
  exprDeps.requires_expression( tTag_ );
  exprDeps.requires_expression( pTag_ );
  exprDeps.requires_expression( mmwTag_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
Density<FieldT>::
bind_fields( const Expr::FieldManagerList& fml )
{
  const typename Expr::FieldMgrSelector<FieldT>::type& fm = fml.field_manager<FieldT>();

  t_ = &fm.field_ref( tTag_ );
  p_ = &fm.field_ref( pTag_ );
  mmw_ = &fm.field_ref( mmwTag_ );

}

//--------------------------------------------------------------------

template< typename FieldT >
void
Density<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  FieldT& rho = this->value();

  rho <<= *p_ * *mmw_ / *t_ / Cantera::GasConstant;
}
//--------------------------------------------------------------------

template< typename FieldT >
Density<FieldT>::
Builder::Builder( const Expr::Tag& resultTag,
         const Expr::Tag& tTag,
         const Expr::Tag& pTag,
         const Expr::Tag& mmwTag )
: ExpressionBuilder( resultTag ),
  tTag_( tTag ),
  pTag_( pTag ),
  mmwTag_( mmwTag )
{}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
Density<FieldT>::
Builder::build() const
{
  return new Density<FieldT>( tTag_, pTag_, mmwTag_);
}

} // namespace pokitt

#endif // Density_Expr_h
