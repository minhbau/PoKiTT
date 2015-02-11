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

#ifndef SpeciesRHS_Expr_h
#define SpeciesRHS_Expr_h

#include <expression/Expression.h>

/**
 *  \class SpeciesRHS
 */
template< typename DivOp >
class SpeciesRHS
 : public Expr::Expression< typename DivOp::DestFieldType >
{
  typedef typename DivOp::SrcFieldType  FluxT; ///< the scalar field type implied by the gradient operator
  typedef typename DivOp::DestFieldType ScalarT;   ///< the flux field type implied by the gradient operator

  const DivOp* divOp_;

  DECLARE_FIELD( FluxT, flux_ )
  DECLARE_FIELD( ScalarT, source_ )
  
  SpeciesRHS( const Expr::Tag& fluxTag,
              const Expr::Tag& sourceTag );

public:

  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a SpeciesRHS expression
     *  @param resultTag the tag for the value that this expression computes
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::Tag& fluxTag,
             const Expr::Tag& sourceTag,
             const int nghost = DEFAULT_NUMBER_OF_GHOSTS );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag fluxTag_, sourceTag_;
  };

  ~SpeciesRHS();
  void bind_operators( const SpatialOps::OperatorDatabase& opDB );
  void evaluate();
};



// ###################################################################
//
//                          Implementation
//
// ###################################################################



template< typename DivOp >
SpeciesRHS<DivOp>::
SpeciesRHS( const Expr::Tag& fluxTag,
            const Expr::Tag& sourceTag )
  : Expr::Expression< typename DivOp::DestFieldType >()
{
  flux_ = this->template create_field_request<FluxT>( fluxTag );
  source_ = this->template create_field_request<ScalarT>( sourceTag );
}

//--------------------------------------------------------------------

template< typename DivOp >
SpeciesRHS<DivOp>::
~SpeciesRHS()
{}

//--------------------------------------------------------------------

template< typename DivOp >
void
SpeciesRHS<DivOp>::
bind_operators( const SpatialOps::OperatorDatabase& opDB )
{
  divOp_ = opDB.retrieve_operator<DivOp>();
}

//--------------------------------------------------------------------

template< typename DivOp >
void
SpeciesRHS<DivOp>::
evaluate()
{
  ScalarT& RHS = this->value();

  const FluxT& flux = flux_->field_ref();
  const ScalarT& source = source_->field_ref();

  RHS <<= - (*divOp_)( flux ) + source;
}

//--------------------------------------------------------------------

template< typename DivOp >
SpeciesRHS<DivOp>::
Builder::Builder( const Expr::Tag& resultTag,
                  const Expr::Tag& fluxTag,
                  const Expr::Tag& sourceTag,
                  const int nghost )
  : ExpressionBuilder( resultTag, nghost ),
    fluxTag_( fluxTag ),
    sourceTag_( sourceTag )
{}

//--------------------------------------------------------------------

template< typename DivOp >
Expr::ExpressionBase*
SpeciesRHS<DivOp>::
Builder::build() const
{
  return new SpeciesRHS<DivOp>( fluxTag_,sourceTag_ );
}


#endif // SpeciesRHS_Expr_h
