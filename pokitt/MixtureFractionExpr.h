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

#ifndef MixtureFractionExpr_h
#define MixtureFractionExpr_h

#include <expression/ExprLib.h>

#include <pokitt/MixtureFraction.h>  // mixture fraction support


namespace PoKiTT{


/**
 * \class MixtureFractionExpr
 * \author James C. Sutherland
 * \brief Calculate the mixture fraction given the species compositions.
 */
template<typename FieldT>
class MixtureFractionExpr
  : public Expr::Expression<FieldT>
{
  DECLARE_VECTOR_OF_FIELDS( FieldT, species_ )

  const MixtureFraction mixfrac_;

  MixtureFractionExpr( const Expr::TagList& specTags,
                       const std::vector<double>& fuelMassFracs,
                       const std::vector<double>& oxidMassFracs )
  : Expr::Expression<FieldT>(),
    mixfrac_( oxidMassFracs, fuelMassFracs, true )
  {
    this->template create_field_vector_request<FieldT>( specTags, species_ );
  }

  ~MixtureFractionExpr(){}

public:

  class Builder : public Expr::ExpressionBuilder
  {
    const Expr::TagList specTags_;
    const std::vector<double> fuel_, oxid_;
  public:
    Builder( const Expr::Tag& result,
             const Expr::TagList& specTags,
             const std::vector<double>& fuelMassFracs,
             const std::vector<double>& oxidMassFracs )
    : ExpressionBuilder(result),
      specTags_( specTags ),
      fuel_( fuelMassFracs ),
      oxid_( oxidMassFracs )
    {}

    Expr::ExpressionBase* build() const{
      return new MixtureFractionExpr<FieldT>( specTags_, fuel_, oxid_ );
    }

  };

  void evaluate(){
    mixfrac_.species_to_mixfrac( species_, this->value() );
  }

};

//--------------------------------------------------------------------

} // namespace PoKiTT

#endif // MixtureFractionExpr_h
