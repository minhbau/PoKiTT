#include "DensityCalculator_Species.h"
#include "JacobianElements.h"

#include <pokitt/CanteraObjects.h>
#include <pokitt/MixtureMolWeight.h>
#include <pokitt/thermo/Enthalpy.h>
#include <pokitt/thermo/Density.h>

#include <pokitt/SpeciesN.h>
#include <spatialops/structured/SpatialFieldStore.h>
#include <spatialops/structured/MatVecOps.h>

#include <expression/ExprPatch.h>
#include <expression/matrix-assembly/Compounds.h>
#include <expression/matrix-assembly/MatrixExpression.h>
#include <expression/matrix-assembly/MapUtilities.h>


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
    suffix_("_densSolve"),
    localFactory_(),
    aMat_( boost::make_shared<DenseMatType         >( "A" ) ),
    bMat_( boost::make_shared<ScaledIdentityMatType>( "B" ) ),
    jac_          ( nullptr ),
    patchPtr_     ( nullptr ),
    treePtr_      ( nullptr ),
    setupHasRun_  ( false   ),
    hGuessTag_    ( hOldTag.name()   + suffix_, Expr::STATE_NONE ),
    rhoGuessTag_  ( rhoOldTag.name() + suffix_, Expr::STATE_NONE ),
    mmwTag_       ( "mixMW"          + suffix_, Expr::STATE_NONE ),
    tGuessTag_    ( tOldTag.name()   + suffix_, Expr::STATE_NONE ),
    pTag_         ( pTag.name()      + suffix_, Expr::STATE_NONE ),
    nSpec_        ( CanteraObjects::number_species()),
    nEq_          ( nSpec_ ), // number of independent species (nSpec - 1) + temperature
    mw_           ( CanteraObjects::molecular_weights() )
{
  std::cout<<"Constructor... ";
  this->set_gpu_runnable(true);

  /*
   * set tags for guesses for species mass fractions, and
   * and calculate inverse of species molecular weights
   */
  mwInv_.clear(); yiGuessTags_.clear();
  for(int i = 0; i<nSpec_; ++i)
  {
     const Expr::Tag yTag(yiOldTags[i].name() + suffix_, Expr::STATE_NONE );
     yiGuessTags_.push_back(yTag);
     mwInv_.push_back(1.0/mw_[i]);
  };

  /*
   * Populate some string vectors. These will be used for tags for
   * \f[ \phi_i \times \frac{\partial \rho}{\partial \phi_j} \f]
   */
  {
    std::vector<std::string> aRowNames, aColNames, jacRowNames, jacColNames;
    const std::string aRowSuffix   = "_times_dRho";
    const std::string aColPrefix   = "d";
    const std::string jacRowPrefix = "jacobian_rho";
    const std::string enthName = "enthalpy";

    // for first n-1 species
    for(int i=0; i<nSpec_-1; ++i)
    {
      aRowNames  .push_back( yiOldTags[i].name() + aRowSuffix );
      aColNames  .push_back( aColPrefix   + yiOldTags[i].name() );
      jacRowNames.push_back( jacRowPrefix + yiOldTags[i].name() );
      jacColNames.push_back( yiOldTags[i].name() );
    }
    //for temperature and enthalpy
    aRowNames  .push_back( enthName     + aRowSuffix      );
    aColNames  .push_back( aColPrefix   + tOldTag.name()  );
    jacRowNames.push_back( jacRowPrefix + enthName );
    jacColNames.push_back( tOldTag.name()                 );

    aMatTags_     = matrix::matrix_tags( aRowNames  ,"_",aColNames  );
    jacobianTags_ = matrix::matrix_tags( jacRowNames,"_",jacColNames);
  }

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
  Expr::FieldManagerList& localfml = patchPtr_->field_manager_list();

  std::set<Expr::ExpressionID> rootIDs;
  typedef typename Expr::PlaceHolder<FieldT>::Builder PlcHldr;

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
    localFactory_.register_expression(new PlcHldr(yiGuessTags_[i]));
  }
  // register placeholders for temperature and pressure
  localFactory_.register_expression(new PlcHldr(tGuessTag_));
  localFactory_.register_expression(new PlcHldr(pTag_     ));

  // register an expression for the nth species, mixture molecular weight, enthalpy, and density.
  localFactory_.register_expression(new typename pokitt::SpeciesN        <FieldT>::Builder(yiGuessTags_[nSpec_-1], yiGuessTags_     ));
  localFactory_.register_expression(new typename pokitt::MixtureMolWeight<FieldT>::Builder(mmwTag_, yiGuessTags_, pokitt::MASS      ));
  localFactory_.register_expression(new typename pokitt::Enthalpy        <FieldT>::Builder(hGuessTag_,tGuessTag_, yiGuessTags_      ));
  localFactory_.register_expression(new typename pokitt::Density         <FieldT>::Builder(rhoGuessTag_, tGuessTag_, pTag_, mmwTag_ ));

  // register expressions for elements of matrices we will need
  rootIDs = register_matrix_element_expressions();

  treePtr_  = new Expr::ExpressionTree( rootIDs, localFactory_, patchPtr_->id(), "species_density_solve" );
  treePtr_->register_fields( localfml );
  treePtr_->lock_fields( localfml );
  localfml.allocate_fields( patchPtr_->field_info() );
  treePtr_->bind_fields( localfml );
  {
	  std::ofstream ofile("inner_tree.dot");
	  treePtr_->write_tree(ofile);
  }
}

//--------------------------------------------------------------------

template< typename FieldT >
std::set<Expr::ExpressionID>
DensityFromSpecies<FieldT>::
register_matrix_element_expressions()
{
  std::set<Expr::ExpressionID> rootIDs;
  //set up a TagList for phi={Y_j,h};
  Expr::TagList phiTags; phiTags.clear();
  phiTags.insert(phiTags.begin(),yiGuessTags_.begin(), yiGuessTags_.end()-1);
  phiTags.push_back(hGuessTag_);

  /*
   * setup aMat_ = Jacobian + rho*Identity.
   * aMat_{i,j} = (phi_i)*d(rho)/d(Y_j) for j in [0,nEq_-2],
   * aMat_{i,j} = (phi_i)*d(rho)/d(T)   for j = nEq_-1,
   * where phi = {Y_j,T}
   */
  for(int i=0; i<nEq_; ++i)
    for(int j=0; j<nEq_; ++j)
    {
      const Expr::Tag tag = aMatTags_[ flat_index(i,j) ];
      aMat_-> template element<FieldT>(i,j) = tag;
    }

  typedef typename PartialJacobian_Species    <FieldT>::Builder DrhoDYj;
  typedef typename PartialJacobian_Temperature<FieldT>::Builder DrhoDT;

  for(int j=0; j<nEq_-1; ++j)
  {
    //populate a TagList corresponding to a column of aMat_, phi_j = Y_j
    Expr::TagList tags; tags.clear();
    for(int i=0; i<nEq_; ++i) tags.push_back(aMatTags_[ flat_index(i,j) ]);

    // register an expression for the jth column vector of aMat_, j < nEq_-1
    localFactory_.register_expression(new DrhoDYj(tags, phiTags, rhoGuessTag_, mmwTag_, mw_[j], mw_[nSpec_-1]) );
  }

  //populate a TagList corresponding to a column of aMat_,  phi_j = T
  Expr::TagList tags; tags.clear();
  for(int i=0; i<nEq_; ++i) tags.push_back(aMatTags_[ flat_index(i,nEq_-1) ]);

  // register an expression for the jth column vector of aMat_, j = nEq_-1
  localFactory_.register_expression(new DrhoDT(tags, phiTags, rhoGuessTag_, tGuessTag_) );


  // setup a matrix bMat_ = rho*Identity
  bMat_-> template scale<FieldT>() = rhoGuessTag_;

  aMat_-> template finalize();
  bMat_-> template finalize();

  // calculate jacobian matrix jac_ = aMat_ - bMat_
  jac_ = aMat_ - bMat_;
  typedef typename matrix::MatrixExpression<FieldT>::Builder MatrixExprBuilder;
  rootIDs.insert(localFactory_.register_expression( new MatrixExprBuilder( jacobianTags_, jac_ ) ));

  return rootIDs;
}

//--------------------------------------------------------------------

template< typename FieldT >
void
DensityFromSpecies<FieldT>::
set_initial_guesses()
{
  Expr::FieldManagerList& localfml = patchPtr_->field_manager_list();

  FieldT& rho      = localfml.field_manager<FieldT>().field_ref( rhoGuessTag_ );
  FieldT& pressure = localfml.field_manager<FieldT>().field_ref( pTag_        );
  FieldT& temp     = localfml.field_manager<FieldT>().field_ref( tGuessTag_   );

  rho      <<= rhoOld_  ->field_ref(); // I don't think i need this...
  pressure <<= pressure_->field_ref();
  temp     <<= temp_    ->field_ref();

  for(int i=0; i<nSpec_-1; ++i)
  {
    FieldT& y = localfml.field_manager<FieldT>().field_ref( yiGuessTags_[i] );
    y <<= yiOld_[i]->field_ref();
  }

}

//--------------------------------------------------------------------

template< typename FieldT >
void
DensityFromSpecies<FieldT>::
calculate_residuals( SpatialOps::FieldVector<FieldT>& residualVec )
{
//  assert( residualVec.size() == nEq_ );

  Expr::FieldManagerList& localfml = patchPtr_->field_manager_list();

  const FieldT& rhoGuess = localfml.field_manager<FieldT>().field_ref( rhoGuessTag_ );
  const FieldT& hGuess   = localfml.field_manager<FieldT>().field_ref( hGuessTag_   );
  const FieldT& rhoH     = rhoH_->field_ref();

  //compute residuals (multiplied by -1) corresponding to nSpec_-1 species being solved for
  for(int i=0; i<nSpec_-1; ++i)
  {
    const FieldT& yiGuess = localfml.field_manager<FieldT>().field_ref( yiGuessTags_[i] );
    const FieldT& rhoYi   = rhoYi_[i]->field_ref();
    residualVec(i) <<= rhoGuess*yiGuess - rhoYi;
  }

  //compute residual (multiplied by -1) corresponding to enthalpy. This should be the last element.
  residualVec(nEq_ -1 ) <<= rhoGuess*hGuess - rhoH;
}

//--------------------------------------------------------------------

template< typename FieldT >
void
DensityFromSpecies<FieldT>::
assemble_jacobian(SpatialOps::FieldMatrix<FieldT>& jacobian )
{
  for( matrix::OrdinalType rowIdx=0; rowIdx < nEq_; ++rowIdx ){
	SpatialOps::FieldVector<FieldT> row( jacobian.row( rowIdx ) );
	jac_->assemble( row, rowIdx, matrix::Place() );
  }
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

  Expr::FieldManagerList& localfml = patchPtr_->field_manager_list();
  // set initial guesses for density, species mass fractions, and temperature;
  set_initial_guesses();

  std::cout<<"\nExecuting tree...\n";
  treePtr_->execute_tree();

  // obtain memory for the jacobian and residual
  std::vector<SpatialOps::SpatFldPtr<FieldT> > resPtrs, jacPtrs;
  for( int i=0; i<nEq_; ++i ){
    resPtrs.push_back( SpatialOps::SpatialFieldStore::get<FieldT>( rho ) );
    for( int j=0; j<nEq_; ++j )
      jacPtrs.push_back( SpatialOps::SpatialFieldStore::get<FieldT>( rho ) );
  }
  FieldVector<FieldT> residualVec( resPtrs );
  FieldMatrix<FieldT> jacobian   ( jacPtrs );

  /*
   * compute residuals (multiplied by -1) corresponding to nSpec_-1 species being
   * solved for and assemble the jacobian matrix
   */
  calculate_residuals( residualVec );
  assemble_jacobian  ( jacobian    );

  // do linear solve and dispatch updates
  residualVec = jacobian.solve( residualVec );



//  // collect residuals
//  for( size_t i=0; i<nvars_; ++i )
//    resVector(i) <<= i+2;
//
//  // build the lhs matrix
//  for( matrix::OrdinalType rowIdx=0; rowIdx < nvars_; ++rowIdx ){
//    SpatialOps::FieldVector<FieldT> row( lhsMatrix.row( rowIdx ) );
//    jac_->assemble( row, rowIdx, matrix::Place() );
//  }
//
//  // do linear solve and dispatch updates
//  resVector = lhsMatrix.solve( resVector );
//
//  std::cout<<"\nelement 0:\n";
//  print_field(resVector(0), std::cout );
//
//  std::cout<<"\nelement 1:\n";
//  print_field(resVector(1), std::cout );


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
//  std::cout<<std::endl;
//  std::cout << jac_->assembler_name() << ": " << '\n';
//  std::cout << jac_->get_structure();
//  std::cout << "--------\n";
}

//--------------------------------------------------------------------

template< typename FieldT >
const int
DensityFromSpecies<FieldT>::
flat_index( const int row, const int col)
{
  return row*nEq_ + col;
}

//--------------------------------------------------------------------

template< typename FieldT >
DensityFromSpecies<FieldT>::
~DensityFromSpecies()
{
  delete patchPtr_;
  delete treePtr_;
}

//====================================================================

// explicit template instantiation
#include <spatialops/structured/FVStaggeredFieldTypes.h>
template class DensityFromSpecies<SpatialOps::SVolField>;
}// namespace WasatchCore
