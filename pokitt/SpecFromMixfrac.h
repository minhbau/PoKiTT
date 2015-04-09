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

#ifndef SpecFromMixfrac_Expr_h
#define SpecFromMixfrac_Expr_h

#include <expression/Expression.h>

#include <pokitt/MixtureFraction.h>

namespace PoKiTT{

/**
 *  \class SpecFromMixfrac
 */
template< typename FieldT >
class SpecFromMixfrac
 : public Expr::Expression<FieldT>
{
  DECLARE_FIELD( FieldT, mixfr_ )

  const MixtureFraction mixfrac_;
  
  SpecFromMixfrac( const Expr::Tag& mixfrTag )
    : Expr::Expression<FieldT>()
  {
    mixfr_ = this->template create_field_request<FieldT>( mixfrTag );
  }


public:

  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a SpecFromMixfrac expression
     *  @param resultTag the tag for the value that this expression computes
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::Tag& mixfrTag,
             const int nghost = DEFAULT_NUMBER_OF_GHOSTS )
    : ExpressionBuilder( resultTag, nghost ),
      mixfrTag_( mixfrTag )
  {}


    Expr::ExpressionBase* build() const{
      return new SpecFromMixfrac<FieldT>( mixfrTag_ );
    }

  private:
    const Expr::Tag mixfrTag_;
  };

  ~SpecFromMixfrac(){}

  void evaluate(){
    mixfrac_.mixfrac_to_species( mixfr_->field_ref(), this->value_vec() );
  }
};

} // namespace PoKiTT

#endif // SpecFromMixfrac_Expr_h
