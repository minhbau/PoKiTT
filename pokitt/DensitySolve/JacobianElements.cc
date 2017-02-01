#include "JacobianElements.h"
#include <pokitt/CanteraObjects.h>
#include <pokitt/MixtureMolWeight.h>
#include <spatialops/structured/SpatialFieldStore.h>
#include <spatialops/structured/MatVecOps.h>
/*
* JacobianElements.cc
 *
 *  Created on: Jan 30, 2017
 *      Author: josh
 */




// ###################################################################
//
//                          Implementation
//
// ###################################################################

template<typename FieldT>
PartialJacobian_Species<FieldT>::
PartialJacobian_Species( const Expr::Tag& rhoGuessTag,
                         const Expr::Tag& phiTag,
                         const Expr::Tag& mmwTag,
                         const double     mw_i,
                         const double     mw_n )
  : Expr::Expression<FieldT>(),
    mwTerm_( 1.0/mw_i - 1.0/mw_n )
{
  this->set_gpu_runnable(true);
  rhoGuess_ = this->template create_field_request<FieldT>( rhoGuessTag );
  phi_      = this->template create_field_request<FieldT>( phiTag );
  mmw_      = this->template create_field_request<FieldT>( mmwTag );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
PartialJacobian_Species<FieldT>::
evaluate()
{
  FieldT& result = this->value();

  const FieldT& rhoGuess = rhoGuess_->field_ref();
  const FieldT& phi      = phi_     ->field_ref();
  const FieldT& mmw      = mmw_     ->field_ref();

  result <<= phi*rhoGuess*mmw*mwTerm_;
}

//--------------------------------------------------------------------

template<typename FieldT>
PartialJacobian_Species<FieldT>::
Builder::Builder( const Expr::Tag& resultTag,
                  const Expr::Tag& rhoGuessTag,
                  const Expr::Tag& phiTag,
                  const Expr::Tag& mmwTag,
                  const double     mw_i,
                  const double     mw_n )
  : ExpressionBuilder( resultTag ),
    rhoOldTag_( rhoGuessTag ),
    phiTag_     ( phiTag      ),
    mmwTag_     ( mmwTag      ),
    mw_i_       ( mw_i        ),
    mw_n_       ( mw_n        )
{}

//====================================================================

// ###################################################################
//
//                          Implementation
//
// ###################################################################

template<typename FieldT>
PartialJacobian_Enthalpy<FieldT>::
PartialJacobian_Enthalpy( const Expr::Tag& rhoGuessTag,
                          const Expr::Tag& phiTag,
                          const Expr::Tag& cpTag,
                          const Expr::Tag& temperatureTag )
: Expr::Expression<FieldT>()
{
  this->set_gpu_runnable(true);
  rhoGuess_ = this->template create_field_request<FieldT>( rhoGuessTag    );
  phi_      = this->template create_field_request<FieldT>( phiTag         );
  cp_       = this->template create_field_request<FieldT>( cpTag          );
  t_        = this->template create_field_request<FieldT>( temperatureTag );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
PartialJacobian_Enthalpy<FieldT>::
evaluate()
{
  FieldT& result = this->value();

  const FieldT& rhoGuess = rhoGuess_->field_ref();
  const FieldT& phi      = phi_     ->field_ref();
  const FieldT& cp       = cp_      ->field_ref();
  const FieldT& t        = t_       ->field_ref();

  result <<= phi*rhoGuess / ( t*cp );
}

//--------------------------------------------------------------------

template<typename FieldT>
PartialJacobian_Enthalpy<FieldT>::
Builder::Builder( const Expr::Tag& resultTag,
                  const Expr::Tag& rhoGuessTag,
                  const Expr::Tag& phiTag,
                  const Expr::Tag& cpTag,
                  const Expr::Tag& temperatureTag )
  : ExpressionBuilder( resultTag ),
    rhoOldTag_   ( rhoGuessTag    ),
    phiTag_        ( phiTag         ),
    cpTag_         ( cpTag          ),
    temperatureTag_( temperatureTag )
{}

// explicit template instantiation
#include <spatialops/structured/FVStaggeredFieldTypes.h>
template class PartialJacobian_Species <SpatialOps::SVolField>;
template class PartialJacobian_Enthalpy<SpatialOps::SVolField>;
