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

#include <expression/ExprLib.h>
#include <expression/TransportEquation.h>
#include <expression/BoundaryConditionExpression.h>

#include <pokitt/transport/ThermalCondMix.h>
#include "EnthalpyRHS.h"
#include "HeatFlux.h"
#include <pokitt/thermo/Enthalpy.h>
#include <pokitt/thermo/Temperature.h>
#include <pokitt/SpeciesN.h>

namespace SO = SpatialOps;
using namespace pokitt;
using Expr::Tag;
using Expr::TagList;
using SO::IntVec;

template< typename FieldT >
class EnthalpyTransport : public Expr::TransportEquation
{

  TagList yiTags_;
  const Tag tTag_;
  const Tag hTag_;

public:

  EnthalpyTransport( Expr::ExpressionFactory& execFactory,
                     const int nSpec,
                     const Tag& hTag,
                     const Tag& tTag,
                     const Tag& rhoTag,
                     const Tag& yiTag,
                     const Tag& jTag,
                     const Tag& mmwTag );

  ~EnthalpyTransport();

  void setup_boundary_conditions( const SO::Grid& grid, Expr::ExpressionFactory& execFactory );

  Expr::ExpressionID initial_condition( Expr::ExpressionFactory& initFactory,
                                        const Tag xTag,
                                        const Tag yTag,
                                        const double tMax,
                                        const double tDev,
                                        const double tMean,
                                        const double tBase );

private:

  void setup_boundary_conditions( Expr::ExpressionFactory& execFactory ){}

};

//====================================================================


//--------------------------------------------------------------------
template< typename FieldT >
EnthalpyTransport<FieldT>::
EnthalpyTransport( Expr::ExpressionFactory& execFactory,
                   const int nSpec,
                   const Tag& hTag,
                   const Tag& tTag,
                   const Tag& rhoTag,
                   const Tag& yiTag,
                   const Tag& jTag,
                   const Tag& mmwTag)
                  : Expr::TransportEquation( hTag.name(),
                    Expr::Tag( hTag, "_RHS" ) ),
                    tTag_( tTag ),
                    hTag_( hTag )
{
  typedef typename SO::FaceTypes<FieldT> FaceTypes;
  typedef typename FaceTypes::XFace XFluxT;
  typedef typename FaceTypes::YFace YFluxT;

  const Tag lamTag( "thermal conductivity", Expr::STATE_NONE);
  const Tag qXTag ( "Heat Flux X",          Expr::STATE_NONE);
  const Tag qYTag ( "Heat Flux Y",          Expr::STATE_NONE);

  yiTags_.clear();
  TagList jxTags, jyTags, hiTags;
  for( size_t n=0; n<nSpec; ++n ){
    std::string spec = boost::lexical_cast<std::string>(n);
    yiTags_.push_back( Tag( yiTag, spec ) );
    hiTags.push_back(  Tag( "h" + spec, Expr::STATE_NONE ) );
    jxTags.push_back(  Tag( jTag, "x" + spec ) );
    jyTags.push_back(  Tag( jTag, "y" + spec ) );
  }

  typedef typename Temperature        <FieldT>::Builder Temperature;
  typedef typename ThermalConductivity<FieldT>::Builder ThermCond;
  typedef typename SpeciesEnthalpy    <FieldT>::Builder SpecEnthalpy;
  typedef typename HeatFlux           <XFluxT>::Builder HeatFluxX;
  typedef typename HeatFlux           <YFluxT>::Builder HeatFluxY;
  typedef typename EnthalpyRHS        <FieldT>::Builder EnthalpyRHS;

  execFactory.register_expression( new Temperature    ( tTag_,     yiTags_, hTag     ) );
  execFactory.register_expression( new ThermCond      ( lamTag,    tTag,    yiTags_, mmwTag        ) );
  execFactory.register_expression( new HeatFluxX      ( qXTag,     tTag_,   lamTag,  hiTags, jxTags ) );
  execFactory.register_expression( new HeatFluxY      ( qYTag,     tTag_,   lamTag,  hiTags, jyTags ) );
  for( size_t n = 0; n < hiTags.size(); ++n ){
    execFactory.register_expression( new SpecEnthalpy ( hiTags[n], tTag_,   n        ) );
  }
  execFactory.register_expression( new EnthalpyRHS ( this->get_rhs_tag(), rhoTag, tag_list( qXTag, qYTag ) ) );
}

//--------------------------------------------------------------------

template< typename FieldT >
EnthalpyTransport<FieldT>::~EnthalpyTransport()
{}

//--------------------------------------------------------------------

template< typename FieldT >
void
EnthalpyTransport<FieldT>::
setup_boundary_conditions( const SO::Grid& grid, Expr::ExpressionFactory& execFactory )
{
  typedef typename SO::FaceTypes<FieldT> FaceTypes;
  typedef typename FaceTypes::XFace XFluxT;
  typedef typename FaceTypes::YFace YFluxT;

  SO::GhostData ghosts( IntVec( 1, 1, 0), IntVec( 1, 1, 0 ) ); // 1 on +-x and +- y and 0 on z
  std::vector<IntVec> xminusPts, xplusPts, yminusPts, yplusPts;

  for( int i=0; i < grid.extent(1); ++i ){
    xminusPts.push_back( IntVec(0,              i, 0) );
    xplusPts. push_back( IntVec(grid.extent(0), i, 0) );
  }
  for(int i = 0; i < grid.extent(0); i++){
    yminusPts.push_back( IntVec(i,              0, 0) );
    yplusPts. push_back( IntVec(i, grid.extent(1), 0) );
  }

  const typename SO::BoundaryCellInfo xInfo = SO::BoundaryCellInfo::build<XFluxT>(true,true,true);
  const typename SO::BoundaryCellInfo yInfo = SO::BoundaryCellInfo::build<YFluxT>(true,true,true);
  const SO::MemoryWindow xwindow( SO::get_window_with_ghost( grid.extent(), ghosts, xInfo ) );
  const SO::MemoryWindow ywindow( SO::get_window_with_ghost( grid.extent(), ghosts, yInfo ) );

  typedef typename Expr::ConstantBCOpExpression<FieldT,SO::XDIR>::Builder XBC;
  typedef typename Expr::ConstantBCOpExpression<FieldT,SO::YDIR>::Builder YBC;

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

  Tag xmbcTag( Tag( tTag_, "_xmbc") );
  Tag xpbcTag( Tag( tTag_, "_xpbc") );
  Tag ymbcTag( Tag( tTag_, "_ymbc") );
  Tag ypbcTag( Tag( tTag_, "_ypbc") );

  execFactory.register_expression( new XBC( xmbcTag, xminus, SO::NEUMANN, SO::MINUS_SIDE,  0.0 ) );
  execFactory.register_expression( new XBC( xpbcTag, xplus , SO::NEUMANN, SO::PLUS_SIDE,  -2e6  ) );
  execFactory.register_expression( new YBC( ymbcTag, yminus, SO::NEUMANN, SO::MINUS_SIDE,  0.0 ) );
  execFactory.register_expression( new YBC( ypbcTag, yplus , SO::NEUMANN, SO::PLUS_SIDE,  -2e6  ) );

  execFactory.attach_modifier_expression( xmbcTag, tTag_ );
  execFactory.attach_modifier_expression( xpbcTag, tTag_ );
  execFactory.attach_modifier_expression( ymbcTag, tTag_ );
  execFactory.attach_modifier_expression( ypbcTag, tTag_ );

  Tag hTag( this->solution_variable_name(), Expr::STATE_N );
  Tag hxmbcTag( Tag( hTag, "_xmbc") );
  Tag hxpbcTag( Tag( hTag, "_xpbc") );
  Tag hymbcTag( Tag( hTag, "_ymbc") );
  Tag hypbcTag( Tag( hTag, "_ypbc") );

  execFactory.register_expression( new XBC( hxmbcTag, xminus, SO::NEUMANN, SO::MINUS_SIDE,  0.0 ) );
  execFactory.register_expression( new XBC( hxpbcTag, xplus , SO::NEUMANN, SO::PLUS_SIDE,   0.0 ) );
  execFactory.register_expression( new YBC( hymbcTag, yminus, SO::NEUMANN, SO::MINUS_SIDE,  0.0 ) );
  execFactory.register_expression( new YBC( hypbcTag, yplus , SO::NEUMANN, SO::PLUS_SIDE,   0.0 ) );

  execFactory.attach_modifier_expression( hxmbcTag, hTag );
  execFactory.attach_modifier_expression( hxpbcTag, hTag );
  execFactory.attach_modifier_expression( hymbcTag, hTag );
  execFactory.attach_modifier_expression( hypbcTag, hTag );
}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionID
EnthalpyTransport<FieldT>::
initial_condition( Expr::ExpressionFactory& initFactory,
                   const Tag xTag,
                   const Tag yTag,
                   const double tMax,
                   const double tDev,
                   const double tMean,
                   const double tBase )
{
  typedef typename SpeciesN                <FieldT>::Builder SpecN;
  typedef typename Expr::GaussianFunction2D<FieldT>::Builder Temperature0;
  typedef typename Expr::PlaceHolder       <FieldT>::Builder XCoord;
  typedef typename Expr::PlaceHolder       <FieldT>::Builder YCoord;
  typedef typename Enthalpy                <FieldT>::Builder Enthalpy;

  initFactory.register_expression(        new XCoord      ( xTag     ) );
  initFactory.register_expression(        new YCoord      ( yTag     ) );
  initFactory.register_expression(        new SpecN       ( yiTags_.back(), yiTags_ ) );
  initFactory.register_expression(        new Temperature0( tTag_,   xTag,  yTag,   tMax-tBase, tDev, tDev, tMean/5, tMean, tBase ) );
  return initFactory.register_expression( new Enthalpy    ( hTag_,   tTag_, yiTags_ ) );
}

//--------------------------------------------------------------------


#endif // EnthalpyTransport_h
