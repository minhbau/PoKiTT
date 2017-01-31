#include "DensityCalculator_Species.h"

#include <pokitt/CanteraObjects.h>
#include <pokitt/MixtureMolWeight.h>
#include <pokitt/thermo/Enthalpy.h>
#include <pokitt/thermo/Density.h>

#include <pokitt/SpeciesN.h>
#include <spatialops/structured/SpatialFieldStore.h>
#include <spatialops/structured/MatVecOps.h>

#include <expression/Functions.h>

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
                    const Expr::Tag&     tOldTag,
                    const Expr::Tag&     pTag,
                    const Expr::TagList& yiOldTags,
                    const Expr::TagList& rhoYiTags )
  : Expr::Expression<FieldT>(),
    localFactory_(),
    patchPtr_     ( nullptr ),
    integratorPtr_( nullptr ),
    setupHasRun_  ( false   ),
    hGuessTag_    ( hOldTag.name()  +"_local", Expr::STATE_NONE ),
    rhoGuessTag_  ( rhoOldTag.name()+"_local", Expr::STATE_NONE ),
    mmwTag_       ( "mixMW_local"            , Expr::STATE_NONE ),
    tGuessTag_    ( tOldTag.name()+"_local"  , Expr::STATE_NONE ),
    pTag_         ( pTag.name()+"_local"     , Expr::STATE_NONE ),
    nSpec_        ( CanteraObjects::number_species()),
    nEq_          ( nSpec_ ), // number of independent species (nSpec - 1) + enthalpy
    mw_           ( CanteraObjects::molecular_weights() )
{
  std::cout<<"Constructor... ";
  this->set_gpu_runnable(true);

  /*
   * set tags for guesses for species mass fractions
   * and calculate inverse of species molecular weights
   */
  mwInv_.clear(); yiGuessTags_.clear();
  for(int i = 0; i<nSpec_; ++i)
  {
     const Expr::Tag yTag(yiOldTags[i].name()+"_local", Expr::STATE_NONE);
     yiGuessTags_.push_back(yTag);
     mwInv_.push_back(1.0/mw_[i]);
  };

  rhoOld_   = this->template create_field_request<FieldT>( rhoOldTag );
  hOld_     = this->template create_field_request<FieldT>( hOldTag   );
  rhoH_     = this->template create_field_request<FieldT>( rhoHTag   );
  temp_     = this->template create_field_request<FieldT>( tOldTag   );
  pressure_ = this->template create_field_request<FieldT>( rhoHTag   );

  this->template create_field_vector_request<FieldT>( yiOldTags, yiOld_ );
  this->template create_field_vector_request<FieldT>( rhoYiTags, rhoYi_ );

  std::cout<<"done\n";
}

//--------------------------------------------------------------------

template<typename FieldT>
DensityFromSpecies<FieldT>::
Builder::Builder( const Expr::Tag&     resultTag,
                  const Expr::Tag&     rhoOldTag,
                  const Expr::Tag&     hOldTag,
                  const Expr::Tag&     rhoHTag,
                  const Expr::Tag&     tOldTag,
                  const Expr::Tag&     pTag,
                  const Expr::TagList& yiOldTags,
                  const Expr::TagList& rhoYiTags,
                  const double         rTol,
                  const unsigned       maxIter )
  : ExpressionBuilder( resultTag ),
    rhoOldTag_  ( rhoOldTag  ),
    hOldTag_    ( hOldTag    ),
    rhoHTag_    ( rhoHTag    ),
    tOldTag_    ( tOldTag    ),
    pTag_       ( pTag       ),
    yiOldTags_  ( yiOldTags  ),
    rhoYiTags_  ( rhoYiTags  ),
    rTol_       ( rTol       ),
    maxIter_    ( maxIter    )
{
  std::cout<<"Builder\n";
}

//--------------------------------------------------------------------

template< typename FieldT >
void
DensityFromSpecies<FieldT>::
setup()
{
  using namespace SpatialOps;
  const IntVec extent = rhoOld_->field_ref().window_with_ghost().extent();

// define a pointer to a patch, integrator, and a local field manager list
  patchPtr_ = new Expr::ExprPatch(extent[0], extent[1], extent[2]);
  integratorPtr_ = new Expr::TimeStepper(localFactory_, Expr::TSMethod::FORWARD_EULER, "local_time_stepper", patchPtr_->id());
  Expr::FieldManagerList& localfml = patchPtr_->field_manager_list();

  Expr::ExpressionTree* tree = integratorPtr_->get_tree();
  typedef typename Expr::PlaceHolder<FieldT>::Builder PlcHldr;
  typedef typename Expr::LinearFunction<FieldT>::Builder Builder;

  /*
   *  register placeholders for guesses for first n-1 species mass fractions.
   *  These are roots of the local tree, as are guesses for density and
   *  temperature.
   *
   *  NOTE: It may be advantageous to include enthalpy in the Newton solve
   *        rather than enthalpy since enthalpy can be computed directly from
   *        heat capacity and temperature, whereas calculating temperature
   *        requires a Newton solve of its own.
   */
  for(int i = 0; i<nSpec_-1; ++i)
  {
    /*
     *  It appears I need to call tree->insert_tree for each root ExprID rather than call
     *  tree->insert_tree for a set of ExprIDs to avoid a segfault when I call
     *  tree->register_fields. I wonder why this is...
     */
    tree->insert_tree(localFactory_.register_expression(new PlcHldr(yiGuessTags_[i])) );
  }
  // register placeholders for temperature and pressure
  tree->insert_tree(localFactory_.register_expression(new PlcHldr(tGuessTag_)) );
  tree->insert_tree(localFactory_.register_expression(new PlcHldr(pTag_     )) );

  // register an expression for the nth species, mixture molecular weight, enthalpy, and density.
  localFactory_.register_expression(new typename pokitt::SpeciesN        <FieldT>::Builder(yiGuessTags_[nSpec_-1], yiGuessTags_     ));
  localFactory_.register_expression(new typename pokitt::MixtureMolWeight<FieldT>::Builder(mmwTag_, yiGuessTags_, pokitt::MASS      ));
  localFactory_.register_expression(new typename pokitt::Enthalpy        <FieldT>::Builder(hGuessTag_,tGuessTag_, yiGuessTags_      ));
  localFactory_.register_expression(new typename pokitt::Density         <FieldT>::Builder(rhoGuessTag_, tGuessTag_, pTag_, mmwTag_ ));

  //==============================================================================
  //==============================================================================
  /*
   * todo: register fields for Jacobian elements and residuals.
   */
  //==============================================================================
  //==============================================================================

  tree->register_fields( localfml );
  tree->lock_fields( localfml );
  localfml.allocate_fields( patchPtr_->field_info() );
  tree->bind_fields( localfml );
  integratorPtr_->finalize( localfml, patchPtr_->operator_database(), patchPtr_->field_info() );
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

  // setup() needs to be run here because we need fields to be defined before it runs
  if(!setupHasRun_){setup(); setupHasRun_=true;}


//  // how are we going to update these?
//  SpatFldPtr<FieldT> mixMW       = SpatialFieldStore::get<FieldT>( rho );
//  SpatFldPtr<FieldT> cp          = SpatialFieldStore::get<FieldT>( rho );
//  SpatFldPtr<FieldT> temperature = SpatialFieldStore::get<FieldT>( rho );
//  *mixMW       <<= 29;
//  *cp          <<= 0.5;
//  *temperature <<= 700;
//
//  FieldVector<FieldT> deltaPhi     (nEq_,window,bcInfo,ghosts, devIndex );
//  FieldVector<FieldT> phi          (nEq_,window,bcInfo,ghosts, devIndex );
//  FieldVector<FieldT> dRhodPhi     (nEq_,window,bcInfo,ghosts, devIndex );
//  FieldVector<FieldT> negOfResidual(nEq_,window,bcInfo,ghosts, devIndex );
//  FieldMatrix<FieldT> jacobian     (nEq_,window,bcInfo,ghosts, devIndex );
//
////  ================== Initial guesses ================== //
//  // density
//  rho <<= rhoOld_->field_ref();
//
//  // species mass fractions
//    for( size_t i=0; i<nSpec_-1; ++i ){
//      phi(i) <<= yiOld_[i]->field_ref();
//    }
//
//    //enthalpy
//    phi(nEq_-1) = hOld_->field_ref();
//
//
////  ================= nonlinear solve =================== //
//
//    //  obtain d(rho)/d(phi_i)
//    drho_dphi( dRhodPhi, mixMW, cp, temperature );
//
//    // obtain the jacobian
//    calc_jac_and_res( jacobian, negOfResidual, dRhodPhi, phi );
//
//    // solve for deltaPhi
//    deltaPhi = jacobian.solve(negOfResidual);
//
//    // update guess for phi
//    phi = phi + deltaPhi;
//
//    // update guess for mmw, density, temperature, etc.

//  ===================================================== //

  rho <<= 42;
  std::cout<<"done\n";
}

//--------------------------------------------------------------------

template< typename FieldT >
DensityFromSpecies<FieldT>::
~DensityFromSpecies()
{
  delete patchPtr_;
  delete integratorPtr_;
}

//====================================================================

// explicit template instantiation
#include <spatialops/structured/FVStaggeredFieldTypes.h>
template class DensityFromSpecies<SpatialOps::SVolField>;
}// namespace WasatchCore
