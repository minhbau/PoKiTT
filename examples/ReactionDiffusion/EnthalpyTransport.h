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

#ifndef EnthalpyTransport_h
#define EnthalpyTransport_h

#include "TagManager.h"
#include "EnthalpyRHS.h"
#include <pokitt/transport/HeatFlux.h>
#include <pokitt/transport/ThermalCondMix.h>
#include <pokitt/thermo/Enthalpy.h>
#include <pokitt/thermo/Temperature.h>
#include <pokitt/SpeciesN.h>
#include <pokitt/MixtureMolWeight.h>

#include <expression/ExprLib.h>
#include <expression/TransportEquation.h>
#include <expression/BoundaryConditionExpression.h>

namespace pokitt{

template< typename FieldT >
class EnthalpyTransport : public Expr::TransportEquation
{

  const TagManager& tagMgr_;

public:

  EnthalpyTransport( Expr::ExpressionFactory& execFactory,
                     const TagManager& tagM );

  ~EnthalpyTransport() {}

  void setup_boundary_conditions( const SpatialOps::Grid& grid, Expr::ExpressionFactory& execFactory );

  Expr::ExpressionID initial_condition( Expr::ExpressionFactory& initFactory,
                                        const double tMax,
                                        const double tDev,
                                        const SpatialOps::DoubleVec tMean,
                                        const SpatialOps::DoubleVec tBase );

private:

  // overwrite base class virtual functions
  Expr::ExpressionID initial_condition( Expr::ExpressionFactory& execFactory ){ return Expr::ExpressionID::null_id(); }
  void setup_boundary_conditions( Expr::ExpressionFactory& execFactory ){}

};

//====================================================================


//--------------------------------------------------------------------
template< typename FieldT >
EnthalpyTransport<FieldT>::
EnthalpyTransport( Expr::ExpressionFactory& execFactory,
                   const TagManager& tagM )
                  : Expr::TransportEquation( tagM[H].name(),
                    Expr::Tag( tagM[H], "_RHS" ) ),
                    tagMgr_( tagM )
{
  typedef typename SpatialOps::FaceTypes<FieldT> FaceTypes;
  typedef typename FaceTypes::XFace XFluxT;
  typedef typename FaceTypes::YFace YFluxT;

  typedef typename Temperature        <FieldT>::Builder Temperature;
  typedef typename ThermalConductivity<FieldT>::Builder ThermCond;
  typedef typename SpeciesEnthalpy    <FieldT>::Builder SpecEnthalpy;
  typedef typename HeatFlux           <XFluxT>::Builder HeatFluxX;
  typedef typename HeatFlux           <YFluxT>::Builder HeatFluxY;
  typedef typename EnthalpyRHS        <FieldT>::Builder EnthalpyRHS;

  const Expr::Tag emptyTag;
  execFactory.register_expression( new Temperature( tagMgr_[T],   tagMgr_[YI_N], tagMgr_[H],    emptyTag                                 ) );
  execFactory.register_expression( new ThermCond  ( tagMgr_[LAM], tagMgr_[T],    tagMgr_[YI_N], tagMgr_[MMW]                ) );
  execFactory.register_expression( new HeatFluxX  ( tagMgr_[QX],  tagMgr_[T],    tagMgr_[LAM],  tagMgr_[H_N], tagMgr_[JX_N] ) );
  execFactory.register_expression( new HeatFluxY  ( tagMgr_[QY],  tagMgr_[T],    tagMgr_[LAM],  tagMgr_[H_N], tagMgr_[JY_N] ) );
  for( size_t n = 0; n < tagMgr_[H_N].size(); ++n ){
    execFactory.register_expression( new SpecEnthalpy ( tagMgr_[H_N][n], tagMgr_[T], n ) );
  }
  execFactory.register_expression( new EnthalpyRHS ( this->get_rhs_tag(), tagMgr_[RHO], tag_list( tagMgr_[QX], tagMgr_[QY] ) ) );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
EnthalpyTransport<FieldT>::
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

  const SpatialOps::IntVec hasBC(true,true,true);
  const typename SpatialOps::BoundaryCellInfo xInfo = SpatialOps::BoundaryCellInfo::build<XFluxT>(hasBC,hasBC);
  const typename SpatialOps::BoundaryCellInfo yInfo = SpatialOps::BoundaryCellInfo::build<YFluxT>(hasBC,hasBC);
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

  Expr::Tag xmbcTag( Expr::Tag( tagMgr_[T], "_xmbc") );
  Expr::Tag xpbcTag( Expr::Tag( tagMgr_[T], "_xpbc") );
  Expr::Tag ymbcTag( Expr::Tag( tagMgr_[T], "_ymbc") );
  Expr::Tag ypbcTag( Expr::Tag( tagMgr_[T], "_ypbc") );

  execFactory.register_expression( new XBC( xmbcTag, xminus, SpatialOps::NEUMANN, SpatialOps::MINUS_SIDE,  0.0 ) );
  execFactory.register_expression( new XBC( xpbcTag, xplus , SpatialOps::NEUMANN, SpatialOps::PLUS_SIDE,   0.0 ) );
  execFactory.register_expression( new YBC( ymbcTag, yminus, SpatialOps::NEUMANN, SpatialOps::MINUS_SIDE,  0.0 ) );
  execFactory.register_expression( new YBC( ypbcTag, yplus , SpatialOps::NEUMANN, SpatialOps::PLUS_SIDE,   0.0 ) );

  execFactory.attach_modifier_expression( xmbcTag, tagMgr_[T] );
  execFactory.attach_modifier_expression( xpbcTag, tagMgr_[T] );
  execFactory.attach_modifier_expression( ymbcTag, tagMgr_[T] );
  execFactory.attach_modifier_expression( ypbcTag, tagMgr_[T] );

  Expr::Tag hxmbcTag( Expr::Tag( tagMgr_[H], "_xmbc") );
  Expr::Tag hxpbcTag( Expr::Tag( tagMgr_[H], "_xpbc") );
  Expr::Tag hymbcTag( Expr::Tag( tagMgr_[H], "_ymbc") );
  Expr::Tag hypbcTag( Expr::Tag( tagMgr_[H], "_ypbc") );

  execFactory.register_expression( new XBC( hxmbcTag, xminus, SpatialOps::NEUMANN, SpatialOps::MINUS_SIDE,  0.0 ) );
  execFactory.register_expression( new XBC( hxpbcTag, xplus , SpatialOps::NEUMANN, SpatialOps::PLUS_SIDE,   0.0 ) );
  execFactory.register_expression( new YBC( hymbcTag, yminus, SpatialOps::NEUMANN, SpatialOps::MINUS_SIDE,  0.0 ) );
  execFactory.register_expression( new YBC( hypbcTag, yplus , SpatialOps::NEUMANN, SpatialOps::PLUS_SIDE,   0.0 ) );

  execFactory.attach_modifier_expression( hxmbcTag, tagMgr_[H] );
  execFactory.attach_modifier_expression( hxpbcTag, tagMgr_[H] );
  execFactory.attach_modifier_expression( hymbcTag, tagMgr_[H] );
  execFactory.attach_modifier_expression( hypbcTag, tagMgr_[H] );
}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionID
EnthalpyTransport<FieldT>::
initial_condition( Expr::ExpressionFactory& initFactory,
                   const double tMax,
                   const double tBase,
                   const SpatialOps::DoubleVec tMean,
                   const SpatialOps::DoubleVec tDev )
{
  typedef typename SpeciesN                <FieldT>::Builder SpecN;
  typedef typename Expr::GaussianFunction2D<FieldT>::Builder Temperature0;
  typedef typename Expr::PlaceHolder       <FieldT>::Builder XCoord;
  typedef typename Expr::PlaceHolder       <FieldT>::Builder YCoord;
  typedef typename Enthalpy                <FieldT>::Builder Enthalpy;

  initFactory.register_expression(        new XCoord      ( tagMgr_[XCOORD] ) );
  initFactory.register_expression(        new YCoord      ( tagMgr_[YCOORD] ) );
  initFactory.register_expression(        new SpecN       ( tagMgr_[YI_N].back(), tagMgr_[YI_N] ) );
  initFactory.register_expression(        new Temperature0( tagMgr_[T], tagMgr_[XCOORD],  tagMgr_[YCOORD], tMax-tBase, tDev[0], tDev[1], tMean[0], tMean[1], tBase ) );
  return initFactory.register_expression( new Enthalpy    ( tagMgr_[H], tagMgr_[T], tagMgr_[YI_N] ) );
}

//--------------------------------------------------------------------

} // namespace pokitt

#endif // EnthalpyTransport_h
