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

#ifndef SpecificVol_Expr_h
#define SpecificVol_Expr_h

#include <expression/Expression.h>

#include <pokitt/CanteraObjects.h> //include cantera wrapper

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

  const double gasConstant_;

  SpecificVol( const Expr::Tag& tTag,
               const Expr::Tag& pTag,
               const Expr::Tag& mmwTag );
public:
  class Builder : public Expr::ExpressionBuilder
  {
    const Expr::Tag tTag_;
    const Expr::Tag pTag_;
    const Expr::Tag mmwTag_;
  public:
    /**
     *  @brief Build a SpecificVol expression
     *  @param resultTag tag for the specific volume
     *  @param tTag tag for temperature
     *  @param pTag tag for pressure
     *  @param mmwTag tag for mixture molecular weight
     *  @param nghost number of ghost cells
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::Tag& tTag,
             const Expr::Tag& pTag,
             const Expr::Tag& mmwTag,
             const SpatialOps::GhostData nghost = DEFAULT_NUMBER_OF_GHOSTS )
    : ExpressionBuilder( resultTag, nghost ),
      tTag_( tTag ),
      pTag_( pTag ),
      mmwTag_( mmwTag )
    {}

    Expr::ExpressionBase* build() const{
      return new SpecificVol<FieldT>( tTag_, pTag_, mmwTag_);
    }
  };

  ~SpecificVol(){};
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
  : Expr::Expression<FieldT>(),
    gasConstant_( CanteraObjects::gas_constant() )
{
  this->set_gpu_runnable( true );

  t_   = this->template create_field_request<FieldT>( tTag   );
  p_   = this->template create_field_request<FieldT>( pTag   );
  mmw_ = this->template create_field_request<FieldT>( mmwTag );
}

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
  nu <<= t * gasConstant_ / ( p * mmw  );
}

//--------------------------------------------------------------------

} // namespace pokitt

#endif // SpecificVol_Expr_h
