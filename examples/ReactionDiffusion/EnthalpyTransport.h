#ifndef EnthalpyTransport_h
#define EnthalpyTransport_h

#include <expression/ExprLib.h>
#include <expression/TransportEquation.h>
#include <expression/BoundaryConditionExpression.h>

#include <pokitt/transport/ThermalCondMix.h>
#include "EnthalpyRHS.h"
#include "HeatFlux.h"
#include <pokitt/thermo/Enthalpy.h>
#include "MassFractions.h"

namespace SO = SpatialOps;
using namespace pokitt;
using Expr::Tag;
using Expr::TagList;

template< typename FieldT >
class EnthalpyTransport : public Expr::TransportEquation
{

public:

  EnthalpyTransport( Expr::ExpressionFactory& execFactory,
                     const Tag& varTag,
                     const Tag& tTag,
                     const Tag& rhoTag,
                     const TagList& yiTags,
                     const Tag& mmwTag,
                     const TagList& jXTags,
                     const TagList& jYTags );

  ~EnthalpyTransport();

  void set_grid( const SO::Grid& grid );

  void setup_boundary_conditions( Expr::ExpressionFactory& execFactory );

  Expr::ExpressionID initial_condition( Expr::ExpressionFactory& initFactory,
                                        const TagList yiTags,
                                        const TagList rhoYiTags,
                                        const Tag rhoTag,
                                        const Tag hTag,
                                        const Tag xTag,
                                        const Tag yTag,
                                        const Tag tTag,
                                        const double rho0,
                                        const double tMax,
                                        const double tDev,
                                        const double tMean,
                                        const double tBase );

};

//====================================================================


//--------------------------------------------------------------------
template< typename FieldT >
EnthalpyTransport<FieldT>::
EnthalpyTransport( Expr::ExpressionFactory& execFactory,
                   const Tag& varTag,
                   const Tag& tTag,
                   const Tag& rhoTag,
                   const TagList& yiTags,
                   const Tag& mmwTag,
                   const TagList& jXTags,
                   const TagList& jYTags )
                  : Expr::TransportEquation( varTag.name(),
                    Expr::Tag( varTag.name()+"RHS", Expr::STATE_N ) )
{
  namespace SO = SpatialOps;
  typedef typename SO::FaceTypes<FieldT> FaceTypes;
  typedef typename FaceTypes::XFace XFluxT;
  typedef typename FaceTypes::YFace YFluxT;

  const Tag lamTag( "thermal conductivity", Expr::STATE_NONE);
  const Tag qXTag ( "Heat Flux X",          Expr::STATE_NONE);
  const Tag qYTag ( "Heat Flux Y",          Expr::STATE_NONE);
  const Tag hRHSTag( varTag.name()+ "RHS",        Expr::STATE_N);
  TagList hiTags;
  for( size_t n=0; n<yiTags.size(); ++n ){
    std::string spec = boost::lexical_cast<std::string>(n);
    hiTags.push_back( Tag( "hi_" + spec, Expr::STATE_NONE ) );
  }

  typedef typename ThermalConductivity<FieldT>::Builder ThermCond;
  typedef typename SpeciesEnthalpy     <FieldT>::Builder SpecEnthalpy;
  typedef typename HeatFlux        <XFluxT>::Builder HeatFluxX;
  typedef typename HeatFlux        <YFluxT>::Builder HeatFluxY;
  typedef typename EnthalpyRHS <FieldT>::Builder EnthalpyRHS;

  execFactory.register_expression( new ThermCond ( lamTag, tTag, yiTags, mmwTag ) );
  for( size_t n = 0; n < hiTags.size(); ++n ){
    execFactory.register_expression( new SpecEnthalpy ( hiTags[n], tTag, n ) );
  }
  execFactory.register_expression( new HeatFluxX( qXTag, tTag, lamTag, hiTags, jXTags ) );
  execFactory.register_expression( new HeatFluxY( qYTag, tTag, lamTag, hiTags, jYTags ) );
  execFactory.register_expression( new EnthalpyRHS ( hRHSTag, rhoTag, tag_list( qXTag, qYTag ) ) );
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
Expr::ExpressionID
EnthalpyTransport<FieldT>::
initial_condition( Expr::ExpressionFactory& initFactory,
                   const TagList yiTags,
                   const TagList rhoYiTags,
                   const Tag rhoTag,
                   const Tag hTag,
                   const Tag xTag,
                   const Tag yTag,
                   const Tag tTag,
                   const double rho0,
                   const double tMax,
                   const double tDev,
                   const double tMean,
                   const double tBase )
{
  typedef typename Expr::PlaceHolder       <FieldT>::Builder XCoord;
  typedef typename Expr::PlaceHolder       <FieldT>::Builder YCoord;
  typedef typename Expr::GaussianFunction2D<FieldT>::Builder Temperature;
  typedef typename MassFractions           <FieldT>::Builder MassFracs;
  typedef typename Enthalpy                <FieldT>::Builder Enthalpy;
  typedef typename Expr::ConstantExpr      <FieldT>::Builder Density;
  typedef typename MassFractions           <FieldT>::Builder MassFracs;

  initFactory.register_expression( new Density( rhoTag, rho0 )    );
  initFactory.register_expression( new MassFracs( yiTags, rhoTag, rhoYiTags) );

  initFactory.register_expression( new XCoord( xTag )    );
  initFactory.register_expression( new YCoord( yTag )    );
  initFactory.register_expression( new Temperature( tTag, xTag, yTag, tMax-tBase, tDev, tDev, tMean, tMean, tBase )    );
  return initFactory.register_expression( new Enthalpy ( hTag , tTag, yiTags ) );
}

//--------------------------------------------------------------------


#endif // EnthalpyTransport_h
