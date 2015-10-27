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

#ifndef LinearMassFracs_Expr_h
#define LinearMassFracs_Expr_h

#include <expression/Expression.h>

/**
 *  \class LinearMassFracs
 */
template< typename FieldT >
class LinearMassFracs
 : public Expr::Expression<FieldT>
{
  DECLARE_FIELD( FieldT, xCoord_ )
  const int nSpec_;
  LinearMassFracs( const Expr::Tag& xCoordTag,
                   const int nSpec );

public:

  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a LinearMassFracs expression
     *  @param resultTag the tag for the value that this expression computes
     *  @param xCoordTag independent variable for linear function
     */
    Builder( const Expr::TagList& resultTags,
             const Expr::Tag& xCoordTag,
             const int nghost = DEFAULT_NUMBER_OF_GHOSTS )
  : ExpressionBuilder( resultTags, nghost ),
    xCoordTag_( xCoordTag ),
    nSpec_( resultTags.size() )
  {}

    Expr::ExpressionBase* build() const{
      return new LinearMassFracs<FieldT>( xCoordTag_, nSpec_ );
    }

  private:
    const Expr::Tag xCoordTag_;
    const int nSpec_;
  };

  ~LinearMassFracs(){}
  void evaluate();
};



// ###################################################################
//
//                          Implementation
//
// ###################################################################



template< typename FieldT >
LinearMassFracs<FieldT>::
LinearMassFracs( const Expr::Tag& xCoordTag,
                 const int nSpec )
  : Expr::Expression<FieldT>(),
    nSpec_( nSpec )
{
  this->set_gpu_runnable( true );
  xCoord_ = this->template create_field_request<FieldT>( xCoordTag );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
LinearMassFracs<FieldT>::
evaluate()
{
  std::vector< FieldT* >& massFracs = this->get_value_vec();

  const FieldT& xCoord = xCoord_->field_ref();
  SpatialOps::SpatFldPtr<FieldT> sum = SpatialOps::SpatialFieldStore::get<FieldT>(xCoord);

  *sum <<= 0.0;
  for( size_t n=0; n<nSpec_; ++n ){
    FieldT& yi = *massFracs[n];
    yi <<= n + 1 + xCoord;
    *sum <<= *sum + yi;
  }
  for( size_t n=0; n<nSpec_; ++n ){
    FieldT& yi = *massFracs[n];
    yi <<= yi / *sum;
  }
}

//--------------------------------------------------------------------

#endif // LinearMassFracs_Expr_h
