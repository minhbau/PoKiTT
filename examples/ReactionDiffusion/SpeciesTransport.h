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

#ifndef SpeciesTransport_h
#define SpeciesTransport_h

#include "TagManager.h"
#include "MassFractions.h"
#include "SpeciesRHS.h"
#include "SpeciesDiffusion.h"
#include <pokitt/thermo/Pressure.h>
#include <pokitt/thermo/Density.h>
#include <pokitt/MixtureMolWeight.h>
#include <pokitt/kinetics/ReactionRates.h>
#include <pokitt/transport/DiffusionCoeffMix.h>

#include <expression/ExprLib.h>
#include <expression/TransportEquation.h>
#include <expression/BoundaryConditionExpression.h>


namespace pokitt{

template< typename FieldT >
class RhoYiIC : public Expr::Expression<FieldT>
{
public:

  struct Builder : public Expr::ExpressionBuilder
  {

    Builder( const Expr::Tag& rhoYi,
             const Expr::Tag& rho,
             const Expr::Tag& yi )
    : Expr::ExpressionBuilder( rhoYi ),
      rho_( rho ),
      yi_( yi )
    {}

    ~Builder(){}
    Expr::ExpressionBase* build() const{
      return new RhoYiIC<FieldT>( rho_, yi_ );
    }

  private:
    const Expr::Tag rho_, yi_;
  };

  void evaluate()
  {
    using namespace SpatialOps;
    this->value() <<= yi_->field_ref() * rho_->field_ref();;
  }

private:
  RhoYiIC( const Expr::Tag& rhoTag,
                  const Expr::Tag& yiTag )
  : Expr::Expression<FieldT>(  )
  {
    this->set_gpu_runnable(true);
    rho_ = this->template create_field_request<FieldT>( rhoTag );
    yi_ = this->template create_field_request<FieldT>( yiTag );
  }

  DECLARE_FIELD( FieldT, rho_ )
  DECLARE_FIELD( FieldT, yi_ )
};

template< typename FieldT >
class SpeciesTransport : public Expr::TransportEquation
{

  const TagManager& tagM_;
  const int n_;

public:

  SpeciesTransport( Expr::ExpressionFactory& execFactory,
                     const TagManager& tagM,
                     const int n);

  ~SpeciesTransport() {};

  void setup_boundary_conditions(const SpatialOps::Grid& grid, Expr::ExpressionFactory& execFactory );

  void register_one_time_expressions( Expr::ExpressionFactory& initFactory, Expr::ExpressionFactory& execFactory );

  Expr::ExpressionID initial_condition( Expr::ExpressionFactory& initFactory,
                                        const double yi,
                                        const double slope );

private:

  // overwrite base class virtual functions
  Expr::ExpressionID initial_condition( Expr::ExpressionFactory& execFactory ){ return Expr::ExpressionID::null_id(); }
  void setup_boundary_conditions( Expr::ExpressionFactory& ){}

};

//====================================================================


//--------------------------------------------------------------------
template< typename FieldT >
SpeciesTransport<FieldT>::
SpeciesTransport( Expr::ExpressionFactory& execFactory,
                   const TagManager& tagM,
                   const int n)
                  : Expr::TransportEquation( tagM[RHOYI_N][n].name(),
                    Expr::Tag( tagM[RHOYI_N][n], "_RHS" ) ),
                    tagM_( tagM ),
                    n_( n )
{
  typedef typename SpeciesRHS <FieldT>::Builder SpeciesRHS;

  execFactory.register_expression( new SpeciesRHS ( this->get_rhs_tag(), tagM_[R_N][n_], tag_list( tagM_[JX_N][n_], tagM_[JY_N][n_] ) ) );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
SpeciesTransport<FieldT>::
setup_boundary_conditions( const SpatialOps::Grid& grid, Expr::ExpressionFactory& execFactory )
{
  typedef typename SpatialOps::FaceTypes<FieldT> FaceTypes;
  typedef typename FaceTypes::XFace XFluxT;
  typedef typename FaceTypes::YFace YFluxT;

  SpatialOps::GhostData ghosts( SpatialOps::IntVec( 1, 1, 0), SpatialOps::IntVec( 1, 1, 0 ) ); // 1 on +-x and +- y and 0 on z
  std::vector<SpatialOps::IntVec> xminusPts, xplusPts, yminusPts, yplusPts;

  for( int i=0; i < grid.extent(1); ++i ){
    xminusPts.push_back( SpatialOps::IntVec(0,              i, 0) );
    xplusPts. push_back( SpatialOps::IntVec(grid.extent(0), i, 0) );
  }
  for(int i = 0; i < grid.extent(0); i++){
    yminusPts.push_back( SpatialOps::IntVec(i,              0, 0) );
    yplusPts. push_back( SpatialOps::IntVec(i, grid.extent(1), 0) );
  }

  const typename SpatialOps::BoundaryCellInfo xInfo = SpatialOps::BoundaryCellInfo::build<XFluxT>(true,true,true);
  const typename SpatialOps::BoundaryCellInfo yInfo = SpatialOps::BoundaryCellInfo::build<YFluxT>(true,true,true);
  const SpatialOps::MemoryWindow xwindow( SpatialOps::get_window_with_ghost( grid.extent(), ghosts, xInfo ) );
  const SpatialOps::MemoryWindow ywindow( SpatialOps::get_window_with_ghost( grid.extent(), ghosts, yInfo ) );

  typedef typename Expr::ConstantBCOpExpression<FieldT,SpatialOps::XDIR>::Builder XBC;
  typedef typename Expr::ConstantBCOpExpression<FieldT,SpatialOps::YDIR>::Builder YBC;

  typename XBC::MaskPtr xminus( new typename XBC::MaskType( xwindow, xInfo, ghosts, xminusPts ) );
  typename XBC::MaskPtr xplus ( new typename XBC::MaskType( xwindow, xInfo, ghosts, xplusPts  ) );
  typename YBC::MaskPtr yminus( new typename YBC::MaskType( ywindow, yInfo, ghosts, yminusPts ) );
  typename YBC::MaskPtr yplus ( new typename YBC::MaskType( ywindow, yInfo, ghosts, yplusPts  ) );

# ifdef ENABLE_CUDA
  // Masks are created on CPU so we need to explicitly transfer them to GPU
  xminus->add_consumer( GPU_INDEX );
  xplus ->add_consumer( GPU_INDEX );
  yminus->add_consumer( GPU_INDEX );
  yplus ->add_consumer( GPU_INDEX );
# endif

  Expr::Tag xmbcTag( this->solution_variable_name() + "_xmbc", Expr::STATE_NONE );
  Expr::Tag xpbcTag( this->solution_variable_name() + "_xpbc", Expr::STATE_NONE );
  Expr::Tag ymbcTag( this->solution_variable_name() + "_ymbc", Expr::STATE_NONE );
  Expr::Tag ypbcTag( this->solution_variable_name() + "_ypbc", Expr::STATE_NONE );

  execFactory.register_expression( new XBC( xmbcTag, xminus, SpatialOps::NEUMANN, SpatialOps::MINUS_SIDE, 0.0 ) );
  execFactory.register_expression( new XBC( xpbcTag, xplus , SpatialOps::NEUMANN, SpatialOps::PLUS_SIDE,  0.0 ) );
  execFactory.register_expression( new YBC( ymbcTag, yminus, SpatialOps::NEUMANN, SpatialOps::MINUS_SIDE, 0.0 ) );
  execFactory.register_expression( new YBC( ypbcTag, yplus , SpatialOps::NEUMANN, SpatialOps::PLUS_SIDE,  0.0 ) );

  execFactory.attach_modifier_expression( xmbcTag, Expr::Tag( this->solution_variable_name(), Expr::STATE_N ) );
  execFactory.attach_modifier_expression( xpbcTag, Expr::Tag( this->solution_variable_name(), Expr::STATE_N ) );
  execFactory.attach_modifier_expression( ymbcTag, Expr::Tag( this->solution_variable_name(), Expr::STATE_N ) );
  execFactory.attach_modifier_expression( ypbcTag, Expr::Tag( this->solution_variable_name(), Expr::STATE_N ) );
}

template< typename FieldT >
Expr::ExpressionID
SpeciesTransport<FieldT>::
initial_condition( Expr::ExpressionFactory& initFactory,
                   const double yi,
                   const double slope )
{
  typedef typename RhoYiIC <FieldT>::Builder RhoYi;
  typedef typename Expr::LinearFunction <FieldT>::Builder Yi;

  initFactory.register_expression(   new Yi(  tagM_[YI_N][n_], tagM_[XCOORD], slope, yi   ) );
  return initFactory.register_expression( new RhoYi ( Expr::Tag( this->solution_variable_name(), Expr::STATE_N), tagM_[YI_N][n_], tagM_[RHO] ) );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
SpeciesTransport<FieldT>::
register_one_time_expressions( Expr::ExpressionFactory& initFactory, Expr::ExpressionFactory& execFactory )
{
  typedef typename Expr::ConstantExpr <FieldT>::Builder P0;
  typedef typename Density< FieldT >::Builder Rho;
  typedef typename MixtureMolWeight <FieldT>::Builder MixMolWeight;
  initFactory.register_expression( new P0         ( tagM_[P], 101000   ) );
  initFactory.register_expression( new MixMolWeight( tagM_[MMW], tagM_[YI_N]          ) );
  initFactory.register_expression( new Rho  ( tagM_[RHO], tagM_[T], tagM_[P], tagM_[MMW]) );

  if( execFactory.have_entry(tagM_[RHO]) )
    return;

  typedef typename SpatialOps::FaceTypes<FieldT> FaceTypes;
  typedef typename FaceTypes::XFace XFluxT;
  typedef typename FaceTypes::YFace YFluxT;

  typedef typename Expr::PlaceHolder<FieldT>::Builder Density;
  typedef typename MassFractions    <FieldT>::Builder MassFracs;
  typedef typename Pressure         <FieldT>::Builder Pressure;
  typedef typename DiffusionCoeff   <FieldT>::Builder DiffusionCoeffs;
  typedef typename ReactionRates    <FieldT>::Builder ReactionRates;
  typedef typename SpeciesDiffFlux  <XFluxT>::Builder MassFluxX;
  typedef typename SpeciesDiffFlux  <YFluxT>::Builder MassFluxY;

  execFactory.register_expression( new Density         ( tagM_[RHO]   ) );
  execFactory.register_expression( new MixMolWeight    ( tagM_[MMW],  tagM_[YI_N]  ) );
  execFactory.register_expression( new MassFracs       ( tagM_[YI_N], tagM_[RHO],  tagM_[RHOYI_N]         ) );
  execFactory.register_expression( new Pressure        ( tagM_[P],    tagM_[T],    tagM_[RHO], tagM_[MMW] ) );
  execFactory.register_expression( new DiffusionCoeffs ( tagM_[D_N],  tagM_[T],    tagM_[P],   tagM_[YI_N], tagM_[MMW] ) );
  execFactory.register_expression( new ReactionRates   ( tagM_[R_N],  tagM_[T],    tagM_[RHO],   tagM_[YI_N], tagM_[MMW] ) );
  execFactory.register_expression( new MassFluxX       ( tagM_[JX_N], tagM_[YI_N], tagM_[RHO], tagM_[MMW],  tagM_[D_N] ) );
  execFactory.register_expression( new MassFluxY       ( tagM_[JY_N], tagM_[YI_N], tagM_[RHO], tagM_[MMW],  tagM_[D_N] ) );
}

//--------------------------------------------------------------------

} // namespace pokitt

#endif // SpeciesTransport_h
