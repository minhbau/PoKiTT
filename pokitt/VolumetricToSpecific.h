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
 *  \file   VolumetricToSpecific.h
 *  \date   Oct 20, 2016
 *  \author mike
 */

#ifndef VOLUMETRICTOSPECIFIC_H_
#define VOLUMETRICTOSPECIFIC_H_

#include <expression/Expression.h>

namespace pokitt{

  /*
   * class: VolumetricToSpecific
   *
   * Convert a volumetric quantity \f$\rho\phi\f$ to a specific quantity \f$\phi\f$ using the density \f$\rho\f$
   */
  template< typename FieldT >
  class VolumetricToSpecific : public Expr::Expression<FieldT>
  {
  private:
    DECLARE_FIELD( FieldT, rhoPhi_ )
    DECLARE_FIELD( FieldT, rho_ )

    VolumetricToSpecific( const Expr::Tag& rhoPhiTag, const Expr::Tag& rhoTag )
    : Expr::Expression<FieldT>()
      {
      this->set_gpu_runnable(true);
      rhoPhi_ = this->template create_field_request<FieldT>( rhoPhiTag );
      rho_ = this->template create_field_request<FieldT>( rhoTag );
      }

  public:

    class Builder : public Expr::ExpressionBuilder
    {
      const Expr::Tag rhoTag_, rhoPhiTag_;
    public:
      Builder( const Expr::Tag& phiTag,
               const Expr::Tag& rhoTag,
               const Expr::Tag& rhoPhiTag,
               const SpatialOps::GhostData nghost = DEFAULT_NUMBER_OF_GHOSTS )
      : Expr::ExpressionBuilder( phiTag ),
        rhoTag_( rhoTag ),
        rhoPhiTag_( rhoPhiTag )
      {}

      ~Builder(){}
      Expr::ExpressionBase* build() const{
        return new VolumetricToSpecific<FieldT>( rhoPhiTag_, rhoTag_ );
      }
    };


    void evaluate()
    {
      using namespace SpatialOps;
      this->value() <<= rhoPhi_->field_ref() / rho_->field_ref();
    }

    void sensitivity( const Expr::Tag& var )
    {
      FieldT& dfdv = this->sensitivity_result( var );
      if( var == this->get_tag() ){
        dfdv <<= 1.0;
      }
      else{
        const FieldT& rho = rho_->field_ref();
        const FieldT& drhodv = this->rho_->sens_field_ref( var );
        const FieldT& drhophidv = this->rhoPhi_->sens_field_ref( var );
        dfdv <<= drhophidv / rho - rhoPhi_->field_ref() * drhodv / ( rho * rho );
      }
    }

  };

} // namespace pokitt



#endif /* VOLUMETRICTOSPECIFIC_H_ */
