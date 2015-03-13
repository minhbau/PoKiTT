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

#include <expression/ExprLib.h>
#include <expression/TransportEquation.h>

#include "SpeciesRHS.h"
#include "MassFractions.h"
#include "SpeciesDiffusion.h"
#include "TagManager.h"
#include <pokitt/thermo/Pressure.h>
#include <pokitt/MixtureMolWeight.h>
#include <pokitt/kinetics/ReactionRates.h>
#include <pokitt/transport/DiffusionCoeffMix.h>

using namespace pokitt;
using Expr::Tag;
using Expr::TagList;

template< typename FieldT >
class SpeciesTransport : public Expr::TransportEquation
{

  const TagManager& tagM_;
  const int n_;

public:

  SpeciesTransport( Expr::ExpressionFactory& execFactory,
                     const TagManager& tagM,
                     const int n);

  ~SpeciesTransport();

  void setup_boundary_conditions(const SO::Grid& grid, Expr::ExpressionFactory& execFactory );

  void register_one_time_expressions(Expr::ExpressionFactory& execFactory );

  Expr::ExpressionID initial_condition( Expr::ExpressionFactory& initFactory,
                                        const double yi,
                                        const double rho0);

private:

  void setup_boundary_conditions( Expr::ExpressionFactory& ){}

};

//====================================================================


//--------------------------------------------------------------------
template< typename FieldT >
SpeciesTransport<FieldT>::
SpeciesTransport( Expr::ExpressionFactory& execFactory,
                   const TagManager& tagM,
                   const int n)
                  : Expr::TransportEquation( tagM.rhoYiN(n).name(),
                    Tag( tagM.rhoYiN(n), "_RHS" ) ),
                    tagM_( tagM ),
                    n_( n )
{
  typedef typename SpeciesRHS <FieldT>::Builder SpeciesRHS;

  execFactory.register_expression( new SpeciesRHS ( this->get_rhs_tag(), tagM_.rN(n_), tag_list( tagM_.jxN(n_), tagM_.jyN(n_) ) ) );
}

//--------------------------------------------------------------------

template< typename FieldT >
SpeciesTransport<FieldT>::~SpeciesTransport()
{}

//--------------------------------------------------------------------

template< typename FieldT >
void
SpeciesTransport<FieldT>::
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

  Tag xmbcTag( this->solution_variable_name() + "_xmbc", Expr::STATE_NONE );
  Tag xpbcTag( this->solution_variable_name() + "_xpbc", Expr::STATE_NONE );
  Tag ymbcTag( this->solution_variable_name() + "_ymbc", Expr::STATE_NONE );
  Tag ypbcTag( this->solution_variable_name() + "_ypbc", Expr::STATE_NONE );

  execFactory.register_expression( new XBC( xmbcTag, xminus, SO::NEUMANN, SO::MINUS_SIDE, 0.0 ) );
  execFactory.register_expression( new XBC( xpbcTag, xplus , SO::NEUMANN, SO::PLUS_SIDE,  0.0 ) );
  execFactory.register_expression( new YBC( ymbcTag, yminus, SO::NEUMANN, SO::MINUS_SIDE, 0.0 ) );
  execFactory.register_expression( new YBC( ypbcTag, yplus , SO::NEUMANN, SO::PLUS_SIDE,  0.0 ) );

  execFactory.attach_modifier_expression( xmbcTag, Tag( this->solution_variable_name(), Expr::STATE_N ) );
  execFactory.attach_modifier_expression( xpbcTag, Tag( this->solution_variable_name(), Expr::STATE_N ) );
  execFactory.attach_modifier_expression( ymbcTag, Tag( this->solution_variable_name(), Expr::STATE_N ) );
  execFactory.attach_modifier_expression( ypbcTag, Tag( this->solution_variable_name(), Expr::STATE_N ) );
}

template< typename FieldT >
Expr::ExpressionID
SpeciesTransport<FieldT>::
initial_condition( Expr::ExpressionFactory& initFactory,
                   const double yi,
                   const double rho0 )
{
  typedef typename Expr::LinearFunction <FieldT>::Builder RhoYi;
  typedef typename Expr::ConstantExpr   <FieldT>::Builder Rho;
  typedef typename Expr::ConstantExpr   <FieldT>::Builder Yi;

  if( !initFactory.have_entry(tagM_[RHO]) )
    initFactory.register_expression( new Rho( tagM_[RHO],    rho0 ) );
  initFactory.register_expression(   new Yi(  tagM_.yiN(n_), yi   ) );
  return initFactory.register_expression( new RhoYi ( Tag( this->solution_variable_name(), Expr::STATE_N), tagM_[RHO], yi, 0.0  ) );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
SpeciesTransport<FieldT>::
register_one_time_expressions( Expr::ExpressionFactory& execFactory )
{
  if( execFactory.have_entry(tagM_[RHO]) )
    return;

  typedef typename SpatialOps::FaceTypes<FieldT> FaceTypes;
  typedef typename FaceTypes::XFace XFluxT;
  typedef typename FaceTypes::YFace YFluxT;

  typedef typename Expr::PlaceHolder<FieldT>::Builder Density;
  typedef typename MassFractions    <FieldT>::Builder MassFracs;
  typedef typename MixtureMolWeight <FieldT>::Builder MixMolWeight;
  typedef typename Pressure         <FieldT>::Builder Pressure;
  typedef typename DiffusionCoeff   <FieldT>::Builder DiffusionCoeffs;
  typedef typename ReactionRates    <FieldT>::Builder ReactionRates;
  typedef typename SpeciesDiffFlux  <XFluxT>::Builder MassFluxX;
  typedef typename SpeciesDiffFlux  <YFluxT>::Builder MassFluxY;

  execFactory.register_expression( new Density         ( tagM_[RHO]   ) );
  execFactory.register_expression( new MixMolWeight    ( tagM_[MMW],  tagM_.yiN()  ) );
  execFactory.register_expression( new MassFracs       ( tagM_.yiN(), tagM_[RHO],  tagM_.rhoYiN()         ) );
  execFactory.register_expression( new Pressure        ( tagM_[P],    tagM_[T],    tagM_[RHO], tagM_[MMW] ) );
  execFactory.register_expression( new DiffusionCoeffs ( tagM_.dN(),  tagM_[T],    tagM_[P],   tagM_.yiN(), tagM_[MMW] ) );
  execFactory.register_expression( new ReactionRates   ( tagM_.rN(),  tagM_[T],    tagM_[RHO],   tagM_.yiN(), tagM_[MMW] ) );
  execFactory.register_expression( new MassFluxX       ( tagM_.jxN(), tagM_.yiN(), tagM_[RHO], tagM_[MMW],  tagM_.dN() ) );
  execFactory.register_expression( new MassFluxY       ( tagM_.jyN(), tagM_.yiN(), tagM_[RHO], tagM_[MMW],  tagM_.dN() ) );
}

//--------------------------------------------------------------------


#endif // SpeciesTransport_h
