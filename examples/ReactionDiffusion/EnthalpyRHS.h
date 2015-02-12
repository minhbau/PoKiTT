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

#ifndef EnthalpyRHS_Expr_h
#define EnthalpyRHS_Expr_h

#include <expression/Expression.h>

/**
 *  \class EnthalpyRHS
 */
template< typename FieldT >
class EnthalpyRHS
 : public Expr::Expression< FieldT >
{
  typedef typename SpatialOps::FaceTypes<FieldT> FaceTypes;
  typedef typename FaceTypes::XFace XFluxT; ///< The type of field for the x-face variables.
  typedef typename FaceTypes::YFace YFluxT; ///< The type of field for the y-face variables.
  typedef typename FaceTypes::ZFace ZFluxT; ///< The type of field for the z-face variables.

  typedef typename SpatialOps::BasicOpTypes<FieldT> OpTypes;
  typedef typename OpTypes::DivX    DivX; ///< Divergence operator (surface integral) in the x-direction
  typedef typename OpTypes::DivY    DivY; ///< Divergence operator (surface integral) in the y-direction
  typedef typename OpTypes::DivZ    DivZ; ///< Divergence operator (surface integral) in the z-direction

  const DivX* divXOp_;
  const DivY* divYOp_;
  const DivZ* divZOp_;

  DECLARE_FIELD( XFluxT, xFlux_ )
  DECLARE_FIELD( YFluxT, yFlux_ )
  DECLARE_FIELD( ZFluxT, zFlux_ )

  DECLARE_FIELD( FieldT, rho_ )
  
  EnthalpyRHS( const Expr::Tag& rhoTag,
               const Expr::TagList& fluxTags );

public:

  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a EnthalpyRHS expression
     *  @param resultTag the tag for the value that this expression computes
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::Tag& rhoTag,
             const Expr::TagList& fluxTags,
             const int nghost = DEFAULT_NUMBER_OF_GHOSTS );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag rhoTag_;
    const Expr::TagList fluxTags_;
  };

  ~EnthalpyRHS();
  void bind_operators( const SpatialOps::OperatorDatabase& opDB );
  void evaluate();
};



// ###################################################################
//
//                          Implementation
//
// ###################################################################



template< typename FieldT >
EnthalpyRHS<FieldT>::
EnthalpyRHS( const Expr::Tag& rhoTag,
             const Expr::TagList& fluxTags )
  : Expr::Expression<FieldT>()
{
  xFlux_ = this->template create_field_request<XFluxT>( fluxTags[0] );
  yFlux_ = this->template create_field_request<YFluxT>( fluxTags[1] );
  zFlux_ = this->template create_field_request<ZFluxT>( fluxTags[2] );

  rho_ = this->template create_field_request<FieldT>( rhoTag );
}

//--------------------------------------------------------------------

template< typename FieldT >
EnthalpyRHS<FieldT>::
~EnthalpyRHS()
{}

//--------------------------------------------------------------------

template< typename FieldT >
void
EnthalpyRHS<FieldT>::
bind_operators( const SpatialOps::OperatorDatabase& opDB )
{
  divXOp_ = opDB.retrieve_operator<DivX>();
  divYOp_ = opDB.retrieve_operator<DivY>();
  divZOp_ = opDB.retrieve_operator<DivZ>();
}

//--------------------------------------------------------------------

template< typename FieldT >
void
EnthalpyRHS<FieldT>::
evaluate()
{
  FieldT& RHS = this->value();

  const XFluxT& xFlux = xFlux_->field_ref();
  const YFluxT& yFlux = yFlux_->field_ref();
  const ZFluxT& zFlux = zFlux_->field_ref();

  const FieldT& rho = rho_->field_ref();

  RHS <<= ( - (*divXOp_)( xFlux ) - (*divYOp_)( yFlux ) - (*divZOp_)( zFlux ) ) / rho;
//  RHS <<= 0.0;
}

//--------------------------------------------------------------------

template< typename FieldT >
EnthalpyRHS<FieldT>::
Builder::Builder( const Expr::Tag& resultTag,
                  const Expr::Tag& rhoTag,
                  const Expr::TagList& fluxTags,
                  const int nghost )
  : ExpressionBuilder( resultTag, nghost ),
    rhoTag_( rhoTag ),
    fluxTags_( fluxTags )
{}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
EnthalpyRHS<FieldT>::
Builder::build() const
{
  return new EnthalpyRHS<FieldT>( rhoTag_, fluxTags_ );
}


#endif // EnthalpyRHS_Expr_h
