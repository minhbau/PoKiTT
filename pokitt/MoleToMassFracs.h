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

#ifndef MOLETOMASSFRACS_H_
#define MOLETOMASSFRACS_H_

#include <expression/Expression.h>
#include <pokitt/CanteraObjects.h>

namespace pokitt{

/**
 *  \class MoleToMassFracs
 */
template< typename FieldT >
class MoleToMassFracs : public Expr::Expression<FieldT>
{
  std::vector<double> molecularWeights_;
  const int nSpec_;

  DECLARE_VECTOR_OF_FIELDS( FieldT, moleFracs_ )
  DECLARE_FIELD( FieldT, mixtureMolWeight_ )

  MoleToMassFracs( const Expr::TagList& massFracTags, const Expr::Tag& mmwTag );
public:
  class Builder : public Expr::ExpressionBuilder
  {
    const Expr::TagList moleFracTags_;
    const Expr::Tag     mmwTag_;
  public:
    /**
     *  @brief Build a MoleToMassFracs expression
     *  @param resultTags the tags for the mass fractions
     *  @param moleFracTags species mole fractions
     *  @param mmwTag the mixture molecular weight tag
     */
    Builder( const Expr::TagList& resultTags,
             const Expr::TagList& moleFracTags,
             const Expr::Tag&     mmwTag,
             const SpatialOps::GhostData nghost = DEFAULT_NUMBER_OF_GHOSTS )
    : ExpressionBuilder( resultTags. nghost ),
      moleFracTags_( moleFracTags ),
      mmwTag_( mmwTag )
    {}
    Expr::ExpressionBase* build() const{
      return new MoleToMassFracs<FieldT>( moleFracTags_, mmwTag_ );
    }
  };

  ~MoleToMassFracs(){}
  void evaluate();
};



// ###################################################################
//
//                          Implementation
//
// ###################################################################



template< typename FieldT >
MoleToMassFracs<FieldT>::
MoleToMassFracs( const Expr::TagList& moleFracTags,
                 const Expr::Tag&     mmwTag )
  : Expr::Expression<FieldT>(),
    nSpec_( CanteraObjects::number_species() )
{
  this->set_gpu_runnable( true );

  assert( moleFracTags.size() == nSpec_ );

  const std::vector<double>& specMW = CanteraObjects::molecular_weights();
  molecularWeights_.resize(nSpec_);
  for( size_t n=0; n<nSpec_; ++n){
    molecularWeights_[n] = specMW[n];
  }

  mixtureMolWeight_ = this->template create_field_request<FieldT>( mmwTag );
  this->template create_field_vector_request<FieldT>( moleFracTags, moleFracs_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
MoleToMassFracs<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  typename Expr::Expression<FieldT>::ValVec& Yi = this->get_value_vec();
  const FieldT& mmw = mixtureMolWeight_->field_ref();
  for( size_t n=0; n<nSpec_; ++n )
    *Yi[n] <<= moleFracs_[n]->field_ref() / mmw * molecularWeights_[n];

}

//--------------------------------------------------------------------

} // namespace pokitt



#endif /* MOLETOMASSFRACS_H_ */
