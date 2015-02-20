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

#ifndef HeatFlux_h
#define HeatFlux_h


#include <expression/ExprLib.h>

template< typename FluxT >
class HeatFlux : public Expr::Expression< FluxT >
{
  typedef typename SpatialOps::VolType<FluxT>::VolField ScalarT;
  typedef typename SpatialOps::OperatorTypeBuilder<SpatialOps::Gradient,   ScalarT,FluxT>::type GradT;
  typedef typename SpatialOps::OperatorTypeBuilder<SpatialOps::Interpolant,ScalarT,FluxT>::type InterpT;

public:

  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a HeatFlux expression
     *  @param heatFlux the tag for the heat flux
     *  @param tTag temperature
     *  @param tCondTag thermal conductivity
     *  @param hTags species enthalpies
     *  @param massFluxTags mass flux of each species
     */
    Builder( const Expr::Tag& heatFlux,
             const Expr::Tag& tTag,
             const Expr::Tag& tCondTag,
             const Expr::TagList& hTags,
             const Expr::TagList& massFluxTags,
             const int nghost = DEFAULT_NUMBER_OF_GHOSTS );
    ~Builder(){}
    Expr::ExpressionBase* build() const;
  private:
    const Expr::Tag tTag_, tCondTag_;
    const Expr::TagList hTags_, massFluxTags_;
  };

  void evaluate();
  void bind_operators( const SpatialOps::OperatorDatabase& opDB );

  protected:

  HeatFlux( const Expr::Tag& tTag,
            const Expr::Tag& tCondTag,
            const Expr::TagList& hTags,
            const Expr::TagList& massFluxTags );

  ~HeatFlux(){}

  const GradT*   gradOp_;
  const InterpT* interpOp_;

  DECLARE_FIELDS( ScalarT, t_, tCond_ )

  DECLARE_VECTOR_OF_FIELDS( ScalarT, h_ )
  DECLARE_VECTOR_OF_FIELDS( FluxT, j_ )
};



// ###################################################################
//
//                         Implementation
//
// ###################################################################




//--------------------------------------------------------------------
template<typename FluxT>
HeatFlux<FluxT>::
HeatFlux( const Expr::Tag& tTag,
          const Expr::Tag& tCondTag,
          const Expr::TagList& hTags,
          const Expr::TagList& massFluxTags )
  : Expr::Expression<FluxT>()
{
  t_ = this->template create_field_request<ScalarT>( tTag );
  tCond_ = this->template create_field_request<ScalarT>( tCondTag  );

  this->template create_field_vector_request<ScalarT>( hTags,        h_ );
  this->template create_field_vector_request<FluxT>(   massFluxTags, j_ );
}
//--------------------------------------------------------------------
template<typename FluxT>
void
HeatFlux<FluxT>::bind_operators( const SpatialOps::OperatorDatabase& opDB )
{
  gradOp_   = opDB.retrieve_operator<GradT  >();
  interpOp_ = opDB.retrieve_operator<InterpT>();
}
//--------------------------------------------------------------------
template<typename FluxT>
void
HeatFlux<FluxT>::evaluate()
{
  using namespace SpatialOps;
  FluxT& flux = this->value();

  const ScalarT& t = t_->field_ref();
  const ScalarT& tCond = tCond_->field_ref();

  flux <<= 0.0;
  flux <<= - (*interpOp_) ( tCond ) * (*gradOp_) ( t );

  for( size_t i=0; i<h_.size(); ++i ){
    const ScalarT& h = h_[i]->field_ref();
    const FluxT&   j = j_[i]->field_ref();

    flux <<= flux + (*interpOp_)( h ) * j;
  }
}


//====================================================================


//--------------------------------------------------------------------
template<typename FluxT>
HeatFlux<FluxT>::Builder::
Builder( const Expr::Tag& heatFlux,
         const Expr::Tag& tTag,
         const Expr::Tag& tCondTag,
         const Expr::TagList& hTags,
         const Expr::TagList& massFluxTags,
         const int nghost )
: ExpressionBuilder(heatFlux, nghost),
  tTag_ ( tTag ),
  tCondTag_( tCondTag ),
  hTags_( hTags ),
  massFluxTags_(massFluxTags)
{}
//--------------------------------------------------------------------
template<typename FluxT>
Expr::ExpressionBase*
HeatFlux<FluxT>::Builder::build() const
{
  return new HeatFlux<FluxT>( tTag_, tCondTag_, hTags_, massFluxTags_ );
}

#endif
