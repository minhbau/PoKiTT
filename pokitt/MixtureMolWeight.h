/*
 * The MIT License
 *
 * Copyright (c) 2015-2017 The University of Utah
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
#include <pokitt/CanteraObjects.h>

namespace pokitt{

/**
 * \class FractionType
 */
enum FractionType
{
  MASS,
  MOLE
};


/**
 *  \class MixtureMolWeight
 */
template< typename FieldT >
class MixtureMolWeight
 : public Expr::Expression<FieldT>
{

  std::vector<double> molecularWeights_, molecularWeightsInv_;
  const int nSpec_;
  const FractionType fracType_;

  DECLARE_VECTOR_OF_FIELDS( FieldT, fracs_ )

  MixtureMolWeight( const Expr::TagList& fracTags,
                    const FractionType fracType );
public:
  class Builder : public Expr::ExpressionBuilder
  {
    const Expr::TagList fracTags_;
    const FractionType fracType_;
  public:
    /**
     *  @brief Build a MixtureMolWeight expression
     *  @param resultTag the tag for mixture molecular weight
     *  @param fracTags species mass fraction, indexing is consistent with Cantera input
     *  @param fracType if mole (MOLE) or mass (MASS) fractions are provided
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::TagList& fracTags,
             const FractionType fracType,
             const SpatialOps::GhostData nghost = DEFAULT_NUMBER_OF_GHOSTS )
    : ExpressionBuilder( resultTag, nghost ),
      fracTags_( fracTags ),
      fracType_( fracType )
    {}

    Expr::ExpressionBase* build() const{
      return new MixtureMolWeight<FieldT>( fracTags_, fracType_ );
    }
  };

  ~MixtureMolWeight(){}
  void evaluate();
  void sensitivity( const Expr::Tag& var );
};



// ###################################################################
//
//                          Implementation
//
// ###################################################################



template< typename FieldT >
MixtureMolWeight<FieldT>::
MixtureMolWeight( const Expr::TagList& fracTags,
                  const FractionType fracType )
  : Expr::Expression<FieldT>(),
    nSpec_( CanteraObjects::number_species() ),
    fracType_( fracType )
{
  this->set_gpu_runnable( true );

  assert( fracTags.size() == nSpec_ );

  const std::vector<double>& specMW = CanteraObjects::molecular_weights();
  molecularWeights_.resize(nSpec_);
  molecularWeightsInv_.resize(nSpec_);
  for( size_t n=0; n<nSpec_; ++n){
    molecularWeights_[n]    = specMW[n];
    molecularWeightsInv_[n] = 1 / specMW[n];
  }

  this->template create_field_vector_request<FieldT>( fracTags, fracs_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
MixtureMolWeight<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  FieldT& mixMW = this->value();

  if( fracType_ == MOLE ){
    const FieldT& x0 = fracs_[0]->field_ref();

    mixMW <<= x0 * molecularWeights_[0];
    for( size_t n=1; n<nSpec_; ++n ){
      const FieldT& xi = fracs_[n]->field_ref();
      mixMW <<= mixMW + xi * molecularWeights_[n];
    }
  }
  else{
    const FieldT& y0 = fracs_[0]->field_ref();

    mixMW <<= y0 * molecularWeightsInv_[0];
    for( size_t n=1; n<nSpec_; ++n ){
      const FieldT& yi = fracs_[n]->field_ref();
      mixMW <<= mixMW + yi * molecularWeightsInv_[n];
    }
    mixMW <<= 1.0 / mixMW;
  }
}

//--------------------------------------------------------------------

template< typename FieldT >
void
MixtureMolWeight<FieldT>::
sensitivity( const Expr::Tag& var )
{
  FieldT& mixMW = this->value();
  FieldT& dfdv = this->sensitivity_result( var );
  if( var == this->get_tag() ){
    dfdv <<= 1.0;
  }
  else{
    if( fracType_ == MOLE ){
      dfdv <<= fracs_[0]->sens_field_ref( var ) * (molecularWeights_[0] - molecularWeights_[nSpec_-1]);
      for( size_t n=1; n<nSpec_-1; ++n ){
        dfdv <<= dfdv + fracs_[n]->sens_field_ref( var ) * (molecularWeights_[n] - molecularWeights_[nSpec_-1]);
      }
    }
    else{
      dfdv <<= fracs_[0]->sens_field_ref( var ) * (molecularWeightsInv_[0]- molecularWeightsInv_[nSpec_-1]);
      for( size_t n=1; n<nSpec_-1; ++n ){
        dfdv <<= dfdv + fracs_[n]->sens_field_ref( var ) * (molecularWeightsInv_[n] - molecularWeightsInv_[nSpec_-1]);
      }
      dfdv <<= - mixMW * mixMW * dfdv;
    }
  }
}

//--------------------------------------------------------------------

} // namespace pokitt

#endif // MixtureMolWeight_Expr_h
