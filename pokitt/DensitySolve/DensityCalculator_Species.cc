#include "DensityCalculator_Species.h"

#include <pokitt/CanteraObjects.h>
#include <pokitt/MixtureMolWeight.h>
#include <spatialops/structured/SpatialFieldStore.h>
#include <spatialops/structured/MatVecOps.h>

namespace WasatchCore{

// ###################################################################
//
//                          Implementation
//
// ###################################################################

template<typename FieldT>
DensityFromSpecies<FieldT>::
DensityFromSpecies( const Expr::Tag&     rhoOldTag,
                    const Expr::Tag&     hOldTag,
                    const Expr::Tag&     rhoHTag,
                    const Expr::TagList& yiOldTags,
                    const Expr::TagList& rhoYiTags )
  : Expr::Expression<FieldT>(),
    nSpec_( CanteraObjects::number_species()),
    nEq_     ( nSpec_ ), // number of independent species (nSpec - 1) + enthalpy
    mw_      ( CanteraObjects::molecular_weights() )
{
  std::cout<<"Constructor... ";
  this->set_gpu_runnable(true);

  // calculate inverse of species molecular weights
  mwInv_.clear();
  for(int i = 0; i<nSpec_; ++i){  mwInv_.push_back(1./mw_[i]); };

//  Expr::Tag mmwTag("mmw_Jac", Expr::STATE_NONE);
//  builder_ = new typename pokitt::MixtureMolWeight<FieldT>::Builder(mmwTag,yiOldTags,pokitt::MASS);
//  mmwExpr_ = builder_->build();

  rhoOld_ = this->template create_field_request<FieldT>(rhoOldTag);
  hOld_   = this->template create_field_request<FieldT>(hOldTag  );
  rhoH_   = this->template create_field_request<FieldT>(rhoHTag  );

  this->template create_field_vector_request<FieldT>(yiOldTags, yiOld_);
  this->template create_field_vector_request<FieldT>(rhoYiTags, rhoYi_);
  std::cout<<"done\n";
}

//--------------------------------------------------------------------

template<typename FieldT>
DensityFromSpecies<FieldT>::
Builder::Builder( const Expr::Tag&     resultTag,
                  const Expr::Tag&     rhoOldTag,
                  const Expr::Tag&     hOldTag,
                  const Expr::Tag&     rhoHTag,
                  const Expr::TagList& yiOldTags,
                  const Expr::TagList& rhoYiTags,
                  const double         rTol,
                  const unsigned       maxIter )
  : ExpressionBuilder( resultTag ),
    rhoOldTag_( rhoOldTag ),
    hOldTag_  ( hOldTag   ),
    rhoHTag_  ( rhoHTag   ),
    yiOldTags_( yiOldTags ),
    rhoYiTags_( rhoYiTags ),
    rTol_     ( rTol      ),
    maxIter_  ( maxIter   )
{
  std::cout<<"Builder\n";
}

//--------------------------------------------------------------------

template< typename FieldT >
void
DensityFromSpecies<FieldT>::
evaluate()
{
  std::cout<<"evaluate... ";
  using namespace SpatialOps;
  FieldT& rho = this->value();

  const MemoryWindow&     window   = rho.window_with_ghost();
  const BoundaryCellInfo& bcInfo   = rho.boundary_info();
  const GhostData&        ghosts   = rho.get_ghost_data();
  const short int         devIndex = rho.active_device_index();

  // how are we going to update these?
  SpatFldPtr<FieldT> mixMW       = SpatialFieldStore::get<FieldT>( rho );
  SpatFldPtr<FieldT> cp          = SpatialFieldStore::get<FieldT>( rho );
  SpatFldPtr<FieldT> temperature = SpatialFieldStore::get<FieldT>( rho );
  *mixMW       <<= 29;
  *cp          <<= 0.5;
  *temperature <<= 700;

  FieldVector<FieldT> deltaPhi     (nEq_,window,bcInfo,ghosts, devIndex );
  FieldVector<FieldT> phi          (nEq_,window,bcInfo,ghosts, devIndex );
  FieldVector<FieldT> dRhodPhi     (nEq_,window,bcInfo,ghosts, devIndex );
  FieldVector<FieldT> negOfResidual(nEq_,window,bcInfo,ghosts, devIndex );
  FieldMatrix<FieldT> jacobian     (nEq_,window,bcInfo,ghosts, devIndex );

//  ================== Initial guesses ================== //
  // density
  rho <<= rhoOld_->field_ref();

  // species mass fractions
    for( size_t i=0; i<nSpec_-1; ++i ){
      phi(i) <<= yiOld_[i]->field_ref();
    }

    //enthalpy
    phi(nEq_-1) = hOld_->field_ref();


//  ================= nonlinear solve =================== //

    //  obtain d(rho)/d(phi_i)
    drho_dphi( dRhodPhi, mixMW, cp, temperature );

    // obtain the jacobian
    calc_jac_and_res( jacobian, negOfResidual, dRhodPhi, phi );

    // solve for deltaPhi
    deltaPhi = jacobian.solve(negOfResidual);

    // update guess for phi
    phi = phi + deltaPhi;

    // update guess for mmw, density, temperature, etc.

//  ===================================================== //

  rho <<= 42;
  std::cout<<"done\n";
}

//--------------------------------------------------------------------

// it would be just awesome if I could make the arguments here members of DensityFromSpecies...
template<typename FieldT>
void
DensityFromSpecies<FieldT>::
calc_jac_and_res( SpatialOps::FieldMatrix<FieldT>& jacobian,
                  SpatialOps::FieldVector<FieldT>& negOfResidual,
                  const SpatialOps::FieldVector<FieldT>& dRhodPhi,
                  const SpatialOps::FieldVector<FieldT>& phi )
{
  std::cout<<"jacobian... ";
  FieldT& rho = this->value();

  // residuals multiplied by (-1)
  // species mass fractions
  for( int i = 0; i<nSpec_-1; ++i){
    negOfResidual(i) <<= -rhoYi_[i]->field_ref() + rho * phi(i);
  }
  // enthalpy
  negOfResidual(nEq_-1) <<= -rhoH_->field_ref() + rho * phi(nEq_-1);

  // [Jacobian] - rho*[identity]
  for( int i = 0; i<nEq_-1; ++i){
    for( int j = 0; j<nEq_-1; ++j){

      // non-diagonal elements
      jacobian(i,j) <<= phi(i) * dRhodPhi(j);

      // diagonal elements
      if(i == j){ jacobian(i,j) <<= jacobian(i,j) - rho; }
    }
  }
  std::cout<<"done\n";
}

//--------------------------------------------------------------------

template<typename FieldT>
void
DensityFromSpecies<FieldT>::
drho_dphi( SpatialOps::FieldVector<FieldT>& dRhodPhi,
           const SpatialOps::SpatFldPtr<FieldT> mixMW,
           const SpatialOps::SpatFldPtr<FieldT> cp,
           const SpatialOps::SpatFldPtr<FieldT> temperature )
{
  FieldT& rho = this->value();

  //  d(rho)/d(Y_i)
  for( int i = 0; i<nSpec_-1; ++i){
    dRhodPhi(i) <<= -rho * *mixMW * (mwInv_[i] - mwInv_[nSpec_]);
  }

  //  d(rho)/d(h)
  dRhodPhi(nEq_-1) <<= -rho / ( *temperature * *cp );

}

//====================================================================

// explicit template instantiation
#include <spatialops/structured/FVStaggeredFieldTypes.h>
template class DensityFromSpecies      <SpatialOps::SVolField>;
}// namespace WasatchCore
