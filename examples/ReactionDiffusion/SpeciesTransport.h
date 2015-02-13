#ifndef SpeciesTransport_h
#define SpeciesTransport_h

#include <expression/ExprLib.h>
#include <expression/TransportEquation.h>

#include "SpeciesRHS.h"

using namespace pokitt;
using Expr::Tag;
using Expr::TagList;

template< typename FieldT >
class SpeciesTransport : public Expr::TransportEquation
{
public:

  SpeciesTransport( Expr::ExpressionFactory& execFactory,
                     const Tag& varTag,
                     const Tag& rTag,
                     const Tag& jXTag,
                     const Tag& jYTag );

  ~SpeciesTransport();

  void setup_boundary_conditions( Expr::ExpressionFactory& );

  Expr::ExpressionID initial_condition( Expr::ExpressionFactory& initFactory,
                                        const Tag& rhoTag,
                                        const Tag& rhoYiTag,
                                        const double yi );

};

//====================================================================


//--------------------------------------------------------------------
template< typename FieldT >
SpeciesTransport<FieldT>::
SpeciesTransport( Expr::ExpressionFactory& execFactory,
                   const Tag& varTag,
                   const Tag& rTag,
                   const Tag& jXTag,
                   const Tag& jYTag )
                  : Expr::TransportEquation( varTag.name(),
                    Tag( varTag.name()+"RHS", Expr::STATE_N ) )
{
  const Tag specRHSTag( varTag.name()+"RHS", Expr::STATE_N );

  typedef typename SpeciesRHS <FieldT>::Builder SpeciesRHS;

  execFactory.register_expression( new SpeciesRHS ( specRHSTag, rTag, tag_list( jXTag, jYTag ) ) );
}

//--------------------------------------------------------------------

template< typename FieldT >
SpeciesTransport<FieldT>::~SpeciesTransport()
{}

//--------------------------------------------------------------------

template< typename FieldT >
void
SpeciesTransport<FieldT>::
setup_boundary_conditions( Expr::ExpressionFactory& factory )
{}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionID
SpeciesTransport<FieldT>::
initial_condition( Expr::ExpressionFactory& initFactory,
                   const Tag& rhoTag,
                   const Tag& rhoYiTag,
                   const double yi)
{
  typedef typename Expr::LinearFunction    <FieldT>::Builder rhoYi;

  return initFactory.register_expression( new rhoYi ( rhoYiTag, rhoTag, yi, 0.0  ) );
}

//--------------------------------------------------------------------


#endif // SpeciesTransport_h
