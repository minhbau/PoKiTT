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

  const TagList yiTags_;
  const Tag      tTag_;

public:

  EnthalpyTransport( Expr::ExpressionFactory& execFactory,
                     const Tag& hTag,
                     const Tag& tTag,
                     const Tag& rhoTag,
                     const TagList& yiTags,
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

  void setup_boundary_conditions( Expr::ExpressionFactory& execFactory );

};

//====================================================================


//--------------------------------------------------------------------
template< typename FieldT >
EnthalpyTransport<FieldT>::
EnthalpyTransport( Expr::ExpressionFactory& execFactory,
                   const Tag& hTag,
                   const Tag& tTag,
                   const Tag& rhoTag,
                   const TagList& yiTags,
                   const Tag& jTag,
                   const Tag& mmwTag )
                  : Expr::TransportEquation( hTag.name(),
                    Expr::Tag( hTag.name()+"RHS", Expr::STATE_N ) ),
                    yiTags_( yiTags ),
                    tTag_( tTag )
{
  typedef typename SO::FaceTypes<FieldT> FaceTypes;
  typedef typename FaceTypes::XFace XFluxT;
  typedef typename FaceTypes::YFace YFluxT;

  const Tag lamTag( "thermal conductivity", Expr::STATE_NONE);
  const Tag qXTag ( "Heat Flux X",          Expr::STATE_NONE);
  const Tag qYTag ( "Heat Flux Y",          Expr::STATE_NONE);

  TagList hiTags;
  TagList jxTags, jyTags;
  for( size_t n=0; n<yiTags_.size(); ++n ){
    std::string spec = boost::lexical_cast<std::string>(n);
    hiTags.push_back( Tag( "hi" + spec, Expr::STATE_NONE ) );
    jxTags.push_back( Tag( jTag, "x" + spec ) );
    jyTags.push_back( Tag( jTag, "y" + spec ) );
  }

  typedef typename Temperature        <FieldT>::Builder Temperature;
  typedef typename ThermalConductivity<FieldT>::Builder ThermCond;
  typedef typename SpeciesEnthalpy    <FieldT>::Builder SpecEnthalpy;
  typedef typename HeatFlux           <XFluxT>::Builder HeatFluxX;
  typedef typename HeatFlux           <YFluxT>::Builder HeatFluxY;
  typedef typename EnthalpyRHS        <FieldT>::Builder EnthalpyRHS;

  execFactory.register_expression( new Temperature     ( tTag_,   yiTags_, hTag                  ) );
  execFactory.register_expression( new ThermCond ( lamTag, tTag, yiTags_, mmwTag ) );
  for( size_t n = 0; n < hiTags.size(); ++n ){
    execFactory.register_expression( new SpecEnthalpy ( hiTags[n], tTag, n ) );
  }
  execFactory.register_expression( new HeatFluxX( qXTag, tTag_, lamTag, hiTags, jxTags ) );
  execFactory.register_expression( new HeatFluxY( qYTag, tTag_, lamTag, hiTags, jyTags ) );
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
setup_boundary_conditions( Expr::ExpressionFactory& execFactory )
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

  Tag xmbcTag( tTag_.name() + "xmbc", Expr::STATE_NONE );
  Tag xpbcTag( tTag_.name() + "xpbc", Expr::STATE_NONE );
  Tag ymbcTag( tTag_.name() + "ymbc", Expr::STATE_NONE );
  Tag ypbcTag( tTag_.name() + "ypbc", Expr::STATE_NONE );

  execFactory.register_expression( new XBC( xmbcTag, xminus, SO::NEUMANN, SO::MINUS_SIDE,  0.0 ) );
  execFactory.register_expression( new XBC( xpbcTag, xplus , SO::NEUMANN, SO::PLUS_SIDE,  0.0 ) );
  execFactory.register_expression( new YBC( ymbcTag, yminus, SO::NEUMANN, SO::MINUS_SIDE,  0.0 ) );
  execFactory.register_expression( new YBC( ypbcTag, yplus , SO::NEUMANN, SO::PLUS_SIDE,  0.0 ) );

  execFactory.attach_modifier_expression( xmbcTag, tTag_ );
  execFactory.attach_modifier_expression( xpbcTag, tTag_ );
  execFactory.attach_modifier_expression( ymbcTag, tTag_ );
  execFactory.attach_modifier_expression( ypbcTag, tTag_ );
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
  typedef typename Expr::GaussianFunction2D<FieldT>::Builder Temperature;
  typedef typename Expr::PlaceHolder       <FieldT>::Builder XCoord;
  typedef typename Expr::PlaceHolder       <FieldT>::Builder YCoord;
  typedef typename Enthalpy                <FieldT>::Builder Enthalpy;
  std::cout<<"here042\n";
  Tag t0Tag( "T0", Expr::STATE_NONE );
  initFactory.register_expression( new SpecN( yiTags_, yiTags_.size()-1 ) );
  std::cout<<"here043\n";
  initFactory.register_expression( new XCoord( xTag )    );
  std::cout<<"here044\n";
  initFactory.register_expression( new YCoord( yTag )    );
  std::cout<<"here045\n";
  initFactory.register_expression( new Temperature( t0Tag, xTag, yTag, tMax-tBase, tDev, tDev, tMean, tMean, tBase )    );
  std::cout<<"here046\n";
  return initFactory.register_expression( new Enthalpy ( Tag( this->solution_variable_name(), Expr::STATE_N) , t0Tag, yiTags_ ) );
}

//--------------------------------------------------------------------


#endif // EnthalpyTransport_h
