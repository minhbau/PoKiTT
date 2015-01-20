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
  DECLARE_FIELDS( FieldT, t_, rho_, mmw_ )

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
             const Expr::Tag& mmwTag,
             const int nghost = DEFAULT_NUMBER_OF_GHOSTS );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag tTag_;
    const Expr::Tag rhoTag_;
    const Expr::Tag mmwTag_;
  };

  ~Pressure();
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
  : Expr::Expression<FieldT>()
{
  this->set_gpu_runnable( true );

  t_   = this->template create_field_request<FieldT>(   tTag );
  rho_ = this->template create_field_request<FieldT>( rhoTag );
  mmw_ = this->template create_field_request<FieldT>( mmwTag );
}

//--------------------------------------------------------------------

template< typename FieldT >
Pressure<FieldT>::
~Pressure()
{}

//--------------------------------------------------------------------

template< typename FieldT >
void
Pressure<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  const FieldT& rho = rho_->field_ref();
  const FieldT& t   = t_  ->field_ref();
  const FieldT& mmw = mmw_->field_ref();
  this->value() <<= rho * Cantera::GasConstant * t / mmw;
}
//--------------------------------------------------------------------

template< typename FieldT >
Pressure<FieldT>::
Builder::Builder( const Expr::Tag& resultTag,
                  const Expr::Tag& tTag,
                  const Expr::Tag& rhoTag,
                  const Expr::Tag& mmwTag,
                  const int nghost )
: ExpressionBuilder( resultTag, nghost ),
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
