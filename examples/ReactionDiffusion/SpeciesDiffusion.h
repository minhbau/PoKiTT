#ifndef SpeciesDiffusion_h
#define SpeciesDiffusion_h


#include <expression/ExprLib.h>


/**
 *  \class  SpeciesDiffFlux
 *  \author Nathan Yonkee
 *  \date December, 2014
 *
 *  \brief Calculates the mass diffusive flux for each species relative to the mass averaged velocity [kg/m^2].
 *
 *  This class calculates the mass flux for each species relative to the mass averaged velocity.
 *  Mixture averaged diffusion coefficients are calculated such that
 *
 *  \f[
 *  j_i = -n D_i \nabla x_i
 *  \f]
 *
 *  where \f$ j_i \f$ is the mass flux, \f$ n \f$ is the molar density, \f$ D_i \f$ is the diffusion coefficient
 *  and \f$ x_i \f$ is the gradient of the mol fraction.
 *
 *  Using mass terms on the RHS, this is equivalent to
 *
 *  \f[
 *  j_i = -\frac{\rho}{MMW} D_i \nabla (y_i MMW)
 *  \f]
 *
 *  where \f$ MMW \f$ is the mixture molecular weight.
 *
 *  \tparam ScalarT the type of the scalar field (diffusivity, mole fraction)
 *  \tparam FluxT the type of the flux field.
 */
template< typename ScalarT, typename FluxT >
class SpeciesDiffFlux : public Expr::Expression< FluxT >
{
  typedef typename SpatialOps::OperatorTypeBuilder<SpatialOps::Gradient,   ScalarT,FluxT>::type GradT;
  typedef typename SpatialOps::OperatorTypeBuilder<SpatialOps::Interpolant,ScalarT,FluxT>::type InterpT;

public:

  typedef std::vector<FluxT*> SpecFluxT;

  typedef SpatialOps::SpatFldPtr<  FluxT> FluxPtrT;
  typedef SpatialOps::SpatFldPtr<ScalarT> ScalarPtrT;

  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a SpeciesDiffFlux expression
     *  @param fluxes the tags for the mass diffusive flux
     *  @param massFracs mass fractions
     *  @param densTag mass density
     *  @param mmwTag mixture molecular weight
     *  @param diffCoeffTags diffusion coefficients
     */
    Builder( const Expr::TagList& fluxes,
             const Expr::TagList& massFracs,
             const Expr::Tag& densTag,
             const Expr::Tag& mmwTag,
             const Expr::TagList& diffCoeffTags,
             const int nghost = DEFAULT_NUMBER_OF_GHOSTS );
    ~Builder(){}
    Expr::ExpressionBase* build() const;
  private:
    const Expr::Tag densTag_, mmwTag_;
    const Expr::TagList yiTags_, diffCoeffTags_;
  };

  void evaluate();
  void bind_operators( const SpatialOps::OperatorDatabase& opDB );

  protected:

  SpeciesDiffFlux( const Expr::TagList& massFracs,
                   const Expr::Tag& densTag,
                   const Expr::Tag& mmwTag,
                   const Expr::TagList& diffCoeffTags);

  ~SpeciesDiffFlux(){}

  const GradT*   gradOp_;
  const InterpT* interpOp_;

  DECLARE_FIELDS( ScalarT, density_, mmw_ )

  DECLARE_VECTOR_OF_FIELDS( ScalarT, species_    )
  DECLARE_VECTOR_OF_FIELDS( ScalarT, diffCoeffs_ )
};



// ###################################################################
//
//                         Implementation
//
// ###################################################################




//--------------------------------------------------------------------
template<typename ScalarT,typename FluxT>
SpeciesDiffFlux<ScalarT,FluxT>::
SpeciesDiffFlux( const Expr::TagList& yiTags,
                 const Expr::Tag& densTag,
                 const Expr::Tag& mmwTag,
                 const Expr::TagList& diffCoeffTags )
  : Expr::Expression<FluxT>()
{
  density_ = this->template create_field_request<ScalarT>( densTag );
  mmw_     = this->template create_field_request<ScalarT>( mmwTag  );

  this->template create_field_vector_request<ScalarT>( yiTags,        species_    );
  this->template create_field_vector_request<ScalarT>( diffCoeffTags, diffCoeffs_ );
}
//--------------------------------------------------------------------
template<typename ScalarT,typename FluxT>
void
SpeciesDiffFlux<ScalarT,FluxT>::bind_operators( const SpatialOps::OperatorDatabase& opDB )
{
  gradOp_   = opDB.retrieve_operator<GradT  >();
  interpOp_ = opDB.retrieve_operator<InterpT>();
}
//--------------------------------------------------------------------
template<typename ScalarT,typename FluxT>
void
SpeciesDiffFlux<ScalarT,FluxT>::evaluate()
{
  using namespace SpatialOps;
  SpecFluxT& fluxes = this->get_value_vec();

  const ScalarT& density = density_->field_ref();
  const ScalarT& mmw     = mmw_    ->field_ref();


  ScalarPtrT molConc = SpatialFieldStore::get<ScalarT>( density );
  *molConc <<= density / mmw;

  FluxPtrT fluxSum = SpatialFieldStore::get<FluxT>( *(fluxes[0]) );
  *fluxSum <<= 0.0;

  for( size_t i=0; i<species_.size(); ++i ){
    const ScalarT& yi       = species_   [i]->field_ref();
    const ScalarT& diffCoef = diffCoeffs_[i]->field_ref();
    FluxT& flux = *fluxes[i];

    // convert mass fraction gradients to mole fraction gradients
    // $ M_i $ cancel out when converting from molar flux to mass flux
    flux <<= - (*interpOp_) ( *molConc * diffCoef ) * (*gradOp_) ( yi * mmw );
    *fluxSum <<= *fluxSum + flux;
  }

  // enforce the sum of fluxes to equal 0
  for( size_t i=0; i<species_.size(); ++i ){
    const ScalarT& yi = species_[i]->field_ref();
    FluxT& flux = *fluxes[i];
    flux <<= flux - *fluxSum * (*interpOp_) ( yi );
  }
}


//====================================================================


//--------------------------------------------------------------------
template<typename ScalarT,typename FluxT>
SpeciesDiffFlux<ScalarT,FluxT>::Builder::
Builder( const Expr::TagList& result,
         const Expr::TagList& massFracs,
         const Expr::Tag& densTag,
         const Expr::Tag& mmwTag,
         const Expr::TagList& diffCoeffTags,
         const int nghost )
: ExpressionBuilder(result, nghost),
  densTag_ ( densTag ),
  yiTags_( massFracs    ),
  mmwTag_(mmwTag),
  diffCoeffTags_(diffCoeffTags)
{}
//--------------------------------------------------------------------
template<typename ScalarT,typename FluxT>
Expr::ExpressionBase*
SpeciesDiffFlux<ScalarT,FluxT>::Builder::build() const
{
  return new SpeciesDiffFlux<ScalarT,FluxT>( yiTags_, densTag_, mmwTag_, diffCoeffTags_ );
}

#endif
