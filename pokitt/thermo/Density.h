/*
 * The MIT License
 *
 * Copyright (c) 2015 The University of Utah
 *
  * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

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
  DECLARE_FIELDS( FieldT, t_, p_, mmw_ )

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
     *  @param nghost number of ghost cells
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

  ~Density();
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
  : Expr::Expression<FieldT>()
{
  this->set_gpu_runnable( true );

  t_   = this->template create_field_request<FieldT>( tTag   );
  p_   = this->template create_field_request<FieldT>( pTag   );
  mmw_ = this->template create_field_request<FieldT>( mmwTag );
}

//--------------------------------------------------------------------

template< typename FieldT >
Density<FieldT>::
~Density()
{}

//--------------------------------------------------------------------

template< typename FieldT >
void
Density<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  FieldT& rho = this->value();
  const FieldT& t   =   t_->field_ref();
  const FieldT& p   =   p_->field_ref();
  const FieldT& mmw = mmw_->field_ref();
  rho <<= p * mmw / ( t * Cantera::GasConstant );
}
//--------------------------------------------------------------------

template< typename FieldT >
Density<FieldT>::
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
Density<FieldT>::
Builder::build() const
{
  return new Density<FieldT>( tTag_, pTag_, mmwTag_);
}

} // namespace pokitt

#endif // Density_Expr_h
