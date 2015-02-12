#ifndef HeatFlux_h
#define HeatFlux_h


#include <expression/ExprLib.h>

template< typename ScalarT, typename FluxT >
class HeatFlux : public Expr::Expression< FluxT >
{
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
template<typename ScalarT,typename FluxT>
HeatFlux<ScalarT,FluxT>::
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
template<typename ScalarT,typename FluxT>
void
HeatFlux<ScalarT,FluxT>::bind_operators( const SpatialOps::OperatorDatabase& opDB )
{
  gradOp_   = opDB.retrieve_operator<GradT  >();
  interpOp_ = opDB.retrieve_operator<InterpT>();
}
//--------------------------------------------------------------------
template<typename ScalarT,typename FluxT>
void
HeatFlux<ScalarT,FluxT>::evaluate()
{
  using namespace SpatialOps;
  FluxT& flux = this->value();

  const ScalarT& t = t_->field_ref();
  const ScalarT& tCond = tCond_->field_ref();

  flux <<= - (*interpOp_) ( tCond ) * (*gradOp_) ( t );

  for( size_t i=0; i<h_.size(); ++i ){
    const ScalarT& h = h_[i]->field_ref();
    const FluxT&   j = j_[i]->field_ref();

    flux <<= flux + (*interpOp_)( h ) * j;
  }
}


//====================================================================


//--------------------------------------------------------------------
template<typename ScalarT,typename FluxT>
HeatFlux<ScalarT,FluxT>::Builder::
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
template<typename ScalarT,typename FluxT>
Expr::ExpressionBase*
HeatFlux<ScalarT,FluxT>::Builder::build() const
{
  return new HeatFlux<ScalarT,FluxT>( tTag_, tCondTag_, hTags_, massFluxTags_ );
}

#endif
