#ifndef SpecificVol_Expr_h
#define SpecificVol_Expr_h

#include <expression/Expression.h>

#include <pokitt/CanteraObjects.h> //include cantera wrapper

#include <cantera/kernel/ct_defs.h> // contains value of gas constant

namespace pokitt{

/**
 *  \class  SpecificVol
 *  \author Nathan Yonkee
 *  \date October, 2014
 *
 *  \brief Calculates the specific volume of a mixture using the
 *  ideal gas equation of state [m^3/kg]
 *
 * \f[
 * \nu = \frac{RT}{P MW_{mix}}
 * \f]
 *
 * Where \f$ MW_{mix}\f$ is the mixture molecular weight
 *
 */

template< typename FieldT >
class SpecificVol
    : public Expr::Expression<FieldT>
{
  DECLARE_FIELDS( FieldT, t_, p_, mmw_ )

  SpecificVol( const Expr::Tag& tTag,
               const Expr::Tag& pTag,
               const Expr::Tag& mmwTag );
public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a SpecificVol expression
     *  @param resultTag tag for the specific volume
     *  @param tTag tag for temperature
     *  @param pTag tag for pressure
     *  @param mmwTag tag for mixture molecular weight
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::Tag& tTag,
             const Expr::Tag& pTag,
             const Expr::Tag& mmwTag,
             const int nghost = DEFAULT_NUMBER_OF_GHOSTS );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag tTag_;
    const Expr::Tag pTag_;
    const Expr::Tag mmwTag_;
  };

  ~SpecificVol();
  void evaluate();

};

// ###################################################################
//
//                          Implementation
//
// ###################################################################

template< typename FieldT >
SpecificVol<FieldT>::
SpecificVol( const Expr::Tag& tTag,
             const Expr::Tag& pTag,
             const Expr::Tag& mmwTag )
  : Expr::Expression<FieldT>()
{
  this->set_gpu_runnable( true );

  t_   = this->template create_field_request<FieldT>( tTag   );
  p_   = this->template create_field_request<FieldT>( pTag   );
  mmw_ = this->template create_field_request<FieldT>( mmwTag );
}

//--------------------------------------------------------------------

template< typename FieldT >
SpecificVol<FieldT>::
~SpecificVol()
{}

//--------------------------------------------------------------------

template< typename FieldT >
void
SpecificVol<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  FieldT& nu = this->value();
  const FieldT& t   =   t_->field_ref();
  const FieldT& p   =   p_->field_ref();
  const FieldT& mmw = mmw_->field_ref();
  nu <<= t * Cantera::GasConstant / ( p * mmw  );
}
//--------------------------------------------------------------------

template< typename FieldT >
SpecificVol<FieldT>::
Builder::Builder( const Expr::Tag& resultTag,
                  const Expr::Tag& tTag,
                  const Expr::Tag& pTag,
                  const Expr::Tag& mmwTag,
                  const int nghost )
: ExpressionBuilder( resultTag, nghost ),
  tTag_( tTag ),
  pTag_( pTag ),
  mmwTag_( mmwTag )
{}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
SpecificVol<FieldT>::
Builder::build() const
{
  return new SpecificVol<FieldT>( tTag_, pTag_, mmwTag_);
}

} // namespace pokitt

#endif // SpecificVol_Expr_h