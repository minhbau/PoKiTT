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

  const double gasConstant_;

  Density( const Expr::Tag& tTag,
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
     *  @brief Build a Density expression
     *  @param resultTag tag for the density  (kg/m^3)
     *  @param tTag tag for temperature (K)
     *  @param pTag tag for pressure (Pa)
     *  @param mmwTag tag for mixture molecular weight (kg/kmol)
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
      return new Density<FieldT>( tTag_, pTag_, mmwTag_);
    }
  };

  ~Density();
  void evaluate();
  void sensitivity( const Expr::Tag& var );
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
    gasConstant_( CanteraObjects::gas_constant() )
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
  rho <<= p * mmw / ( t * gasConstant_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
Density<FieldT>::
sensitivity( const Expr::Tag& var )
{
  using namespace SpatialOps;
  if( var == this->get_tag() ){
    this->sensitivity_result( var ) <<= 0.0;
  }
  else{
    this->sensitivity_result( var ) <<= this->value() * ( p_  ->sens_field_ref( var ) / p_  ->field_ref() +
                                                          mmw_->sens_field_ref( var ) / mmw_->field_ref() - /* yes, this is a minus */
                                                          t_  ->sens_field_ref( var ) / t_  ->field_ref() );
  }
}

//--------------------------------------------------------------------


} // namespace pokitt

#endif // Density_Expr_h
