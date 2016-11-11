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

#ifndef MASSTOMOLEFRACS_H_
#define MASSTOMOLEFRACS_H_

#include <expression/Expression.h>
#include <pokitt/CanteraObjects.h>

namespace pokitt{

/**
 *  \class MassToMoleFracs
 */
template< typename FieldT >
class MassToMoleFracs : public Expr::Expression<FieldT>
{
  std::vector<double> molecularWeightsInv_;
  const int nSpec_;

  DECLARE_VECTOR_OF_FIELDS( FieldT, massFracs_ )
  DECLARE_FIELD( FieldT, mixtureMolWeight_ )

  MassToMoleFracs( const Expr::TagList& massFracTags, const Expr::Tag& mmwTag );
public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a MassToMoleFracs expression
     *  @param resultTags the tags for the mole fractions
     *  @param massFracTags species mass fractions
     *  @param mmwTag the mixture molecular weight tag
     */
    Builder( const Expr::TagList& resultTags,
             const Expr::TagList& massFracTags,
             const Expr::Tag&     mmwTag );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::TagList massFracTags_;
    const Expr::Tag     mmwTag_;
  };

  ~MassToMoleFracs(){}
  void evaluate();
};



// ###################################################################
//
//                          Implementation
//
// ###################################################################



template< typename FieldT >
MassToMoleFracs<FieldT>::
MassToMoleFracs( const Expr::TagList& massFracTags,
                 const Expr::Tag&     mmwTag )
  : Expr::Expression<FieldT>(),
    nSpec_( CanteraObjects::number_species() )
{
  this->set_gpu_runnable( true );

  assert( massFracTags.size() == nSpec_ );

  const std::vector<double>& specMW = CanteraObjects::molecular_weights();
  molecularWeightsInv_.resize(nSpec_);
  for( size_t n=0; n<nSpec_; ++n){
    molecularWeightsInv_[n] = 1 / specMW[n];
  }

  mixtureMolWeight_ = this->template create_field_request<FieldT>( mmwTag );
  this->template create_field_vector_request<FieldT>( massFracTags, massFracs_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
MassToMoleFracs<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  typename Expr::Expression<FieldT>::ValVec& Xi = this->get_value_vec();
  const FieldT& mmw = mixtureMolWeight_->field_ref();
  for( size_t n=0; n<nSpec_; ++n )
    *Xi[n] <<= massFracs_[n]->field_ref() * mmw * molecularWeightsInv_[n];
}

//--------------------------------------------------------------------

template< typename FieldT >
MassToMoleFracs<FieldT>::
Builder::Builder( const Expr::TagList& resultTags,
                  const Expr::TagList& massFracTags,
                  const Expr::Tag&     mmwTag )
  : ExpressionBuilder( resultTags ),
    massFracTags_( massFracTags ),
    mmwTag_( mmwTag )
{}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
MassToMoleFracs<FieldT>::
Builder::build() const
{
  return new MassToMoleFracs<FieldT>( massFracTags_, mmwTag_ );
}

} // namespace pokitt



#endif /* MASSTOMOLEFRACS_H_ */
