/*
 * The MIT License
 *
 * Copyright (c) 2016 The University of Utah
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

/**
 *  \file   SpecificToVolumetric.h
 *  \date   Oct 20, 2016
 *  \author mike
 */

#ifndef SPECIFICTOVOLUMETRIC_H_
#define SPECIFICTOVOLUMETRIC_H_

#include <expression/Expression.h>

namespace pokitt{

  /*
   * class: SpecificToVolumetric
   *
   * Convert a specific quantity \f$\phi\f$ to a volumetric quantity \f$\rho\phi\f$ using the density \f$\rho\f$
   */
  template< typename FieldT >
  class SpecificToVolumetric : public Expr::Expression<FieldT>
  {
  private:
    DECLARE_FIELD( FieldT, rho_ )
    DECLARE_FIELD( FieldT, phi_ )

    SpecificToVolumetric( const Expr::Tag& rhoTag, const Expr::Tag& phiTag )
    : Expr::Expression<FieldT>()
      {
      this->set_gpu_runnable(true);
      rho_ = this->template create_field_request<FieldT>( rhoTag );
      phi_ = this->template create_field_request<FieldT>( phiTag );
      }

  public:

    struct Builder : public Expr::ExpressionBuilder
    {

      Builder( const Expr::Tag& rhoPhiTag,
               const Expr::Tag& rhoTag,
               const Expr::Tag& phiTag )
            : Expr::ExpressionBuilder( rhoPhiTag ),
              rhoTag_( rhoTag ),
              phiTag_( phiTag )
      {}

      ~Builder(){}
      Expr::ExpressionBase* build() const{
        return new SpecificToVolumetric<FieldT>( rhoTag_, phiTag_ );
      }

    private:
      const Expr::Tag rhoTag_, phiTag_;
    };


    void evaluate()
    {
      using namespace SpatialOps;
      this->value() <<= rho_->field_ref() * phi_->field_ref();
    }

    void sensitivity( const Expr::Tag& var )
    {
      FieldT& dfdv = this->sensitivity_result( var );
      if( var == this->get_tag() ){
        dfdv <<= 1.0;
      }
      else{
        const FieldT& drhodv = this->rho_->sens_field_ref( var );
        const FieldT& dphidv = this->phi_->sens_field_ref( var );
        dfdv <<= rho_->field_ref() * dphidv + phi_->field_ref() * drhodv;
      }
    }

  };

} // namespace pokitt


#endif /* SPECIFICTOVOLUMETRIC_H_ */
