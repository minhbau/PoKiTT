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

  const Tag rTag_;
  const Tag rhoYi_;
  const Tag jTag_;

public:

  SpeciesTransport( Expr::ExpressionFactory& execFactory,
                     const Tag& rhoYi,
                     const Tag& jTag,
                     const std::string n);

  ~SpeciesTransport();

  void setup_boundary_conditions(const SO::Grid& grid, Expr::ExpressionFactory& execFactory );

  void register_one_time_expressions(Expr::ExpressionFactory& execFactory,
                                     const int nSpec,
                                     const Tag& rhoTag,
                                     const Tag& yiTag,
                                     const Tag& mmwTag,
                                     const Tag& tTag );

  Expr::ExpressionID initial_condition( Expr::ExpressionFactory& initFactory,
                                        const Tag& rhoTag,
                                        const Tag& yiTag,
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
                   const Tag& rhoYi,
                   const Tag& jTag,
                   const std::string n)
                  : Expr::TransportEquation( rhoYi.name() + n,
                    Tag( rhoYi, n +"RHS" ) ),
                    rTag_( Tag( "r", Expr::STATE_NONE)),
                    rhoYi_( rhoYi ),
                    jTag_( jTag )
{
  typedef typename SpeciesRHS <FieldT>::Builder SpeciesRHS;

  execFactory.register_expression( new SpeciesRHS ( this->get_rhs_tag(), Tag( rTag_, n), tag_list( Tag( jTag,"x" +n), Tag( jTag,"y" +n) ) ) );
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
                   const Tag& rhoTag,
                   const Tag& yiTag,
                   const double yi,
                   const double rho0 )
{
  typedef typename Expr::LinearFunction <FieldT>::Builder RhoYi;
  typedef typename Expr::ConstantExpr   <FieldT>::Builder Rho;
  typedef typename Expr::ConstantExpr   <FieldT>::Builder Yi;

  if( !initFactory.have_entry(rhoTag) )
    initFactory.register_expression( new Rho( rhoTag, rho0 ) );
  initFactory.register_expression( new Yi( yiTag, yi ) );
  return initFactory.register_expression( new RhoYi ( Tag( this->solution_variable_name(), Expr::STATE_N), rhoTag, yi, 0.0  ) );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
SpeciesTransport<FieldT>::
register_one_time_expressions( Expr::ExpressionFactory& execFactory,
                               const int nSpec,
                               const Tag& rhoTag,
                               const Tag& yiTag,
                               const Tag& mmwTag,
                               const Tag& tTag )
{
  if( execFactory.have_entry(rhoTag) )
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

  Tag pTag( "pressure", Expr::STATE_NONE );
  TagList jxTags, jyTags, rhoYiTags, rTags, dTags, yiTags;
  for( size_t n=0; n<nSpec; ++n ){
    std::string spec = boost::lexical_cast<std::string>(n);
    jxTags.push_back( Tag(jTag_,"x" + spec ) );
    jyTags.push_back( Tag(jTag_,"y" + spec ) );
    rTags.push_back( Tag(rTag_, spec ) );
    dTags.push_back( Tag("D"+spec, Expr::STATE_NONE ) );
    yiTags.push_back( Tag( yiTag, spec ));
    if( n != (nSpec - 1) ){
      rhoYiTags.push_back( Tag( rhoYi_, spec ) );
    }
  }

  execFactory.register_expression( new Density         ( rhoTag                                ) );
  execFactory.register_expression( new MassFracs       ( yiTags, rhoTag, rhoYiTags             ) );
  execFactory.register_expression( new MixMolWeight    ( mmwTag, yiTags                        ) );
  execFactory.register_expression( new Pressure        ( pTag,   tTag,   rhoTag, mmwTag        ) );
  execFactory.register_expression( new DiffusionCoeffs ( dTags,  tTag,   pTag,   yiTags, mmwTag) );
  execFactory.register_expression( new ReactionRates   ( rTags,  tTag,   pTag,   yiTags, mmwTag) );
  execFactory.register_expression( new MassFluxX       ( jxTags, yiTags, rhoTag, mmwTag, dTags ) );
  execFactory.register_expression( new MassFluxY       ( jyTags, yiTags, rhoTag, mmwTag, dTags ) );
}

//--------------------------------------------------------------------


#endif // SpeciesTransport_h
