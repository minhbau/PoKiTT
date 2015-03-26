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

#ifndef MixtureMolWeight_Expr_h
#define MixtureMolWeight_Expr_h

#include <expression/Expression.h>

#include <pokitt/CanteraObjects.h> //include cantera wrapper

namespace pokitt{

/**
 *  \class MixtureMolWeight
 */
template< typename FieldT >
class MixtureMolWeight
 : public Expr::Expression<FieldT>
{

  std::vector<double> molecularWeightsInv_;
  const int nSpec_;

  DECLARE_VECTOR_OF_FIELDS( FieldT, massFracs_ )

  MixtureMolWeight( const Expr::TagList& massFracTags );
public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a MixtureMolWeight expression
     *  @param resultTag the tag for mixture molecular weight
     *  @param massFracTags species mass fraction, indexing is consistent with Cantera input
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::TagList& massFracTags );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::TagList massFracTags_;
  };

  ~MixtureMolWeight(){}
  void evaluate();
};



// ###################################################################
//
//                          Implementation
//
// ###################################################################



template< typename FieldT >
MixtureMolWeight<FieldT>::
MixtureMolWeight( const Expr::TagList& massFracTags )
  : Expr::Expression<FieldT>(),
    nSpec_( CanteraObjects::number_species() )
{
  this->set_gpu_runnable( true );

  assert( massFracTags.size() == nSpec_ );

  const std::vector<double>& specMW = CanteraObjects::molecular_weights();
  molecularWeightsInv_.resize(nSpec_);
  for( size_t n=0; n<nSpec_; ++n)
    molecularWeightsInv_[n] = 1 / specMW[n];

  this->template create_field_vector_request<FieldT>( massFracTags, massFracs_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
MixtureMolWeight<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  FieldT& mixMW = this->value();

  const FieldT& y0 = massFracs_[0]->field_ref();

  mixMW <<= y0 * molecularWeightsInv_[0];
  for( size_t n=1; n<nSpec_; ++n ){
    const FieldT& yi = massFracs_[n]->field_ref();
    mixMW <<= mixMW + yi * molecularWeightsInv_[n];
  }
  mixMW <<= 1.0 / mixMW;
}

//--------------------------------------------------------------------

template< typename FieldT >
MixtureMolWeight<FieldT>::
Builder::Builder( const Expr::Tag& resultTag,
                  const Expr::TagList& massFracTags )
  : ExpressionBuilder( resultTag ),
    massFracTags_( massFracTags )
{

}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
MixtureMolWeight<FieldT>::
Builder::build() const
{
  return new MixtureMolWeight<FieldT>( massFracTags_ );
}

} // namespace pokitt

#endif // MixtureMolWeight_Expr_h
