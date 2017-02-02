#include "JacobianElements.h"
#include <spatialops/structured/SpatialFieldStore.h>

// ###################################################################
//
//                          Implementation
//
// ###################################################################

template<typename FieldT>
PartialJacobian_Species<FieldT>::
PartialJacobian_Species( const Expr::TagList& phiTags,
                         const Expr::Tag&     rhoTag,
                         const Expr::Tag&     mmwTag,
                         const double         mw_i,
                         const double         mw_n )
  : Expr::Expression<FieldT>(),
    mwTerm_( 1.0/mw_i - 1.0/mw_n )
{
  this->set_gpu_runnable(true);

  this->template create_field_vector_request<FieldT>( phiTags, phi_ );

  rho_ = this->template create_field_request<FieldT>( rhoTag );
  mmw_ = this->template create_field_request<FieldT>( mmwTag );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
PartialJacobian_Species<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  typename Expr::Expression<FieldT>::ValVec&  result = this->get_value_vec();

  const FieldT& rho = rho_->field_ref();
  const FieldT& mmw = mmw_->field_ref();

  // calculate \f$x = -\frac{\partial \rho}{\partial \Y_j}\f$
  SpatFldPtr<FieldT> x = SpatialFieldStore::get<FieldT>( rho );
  *x <<= rho*mmw*mwTerm_;

  /*
   * calculate \f$A_{i,j} = - \phi_i\frac{\partial \rho}{\partial \Y_j}\f$
   * for fixed \f i \f$.
   */
  for(int i=0; i<result.size(); ++i){
  *result[i] <<= ( phi_[i]->field_ref() )*(*x);
  }
}

//--------------------------------------------------------------------

template<typename FieldT>
PartialJacobian_Species<FieldT>::
Builder::Builder( const Expr::TagList& resultTags,
                  const Expr::TagList& phiTags,
                  const Expr::Tag& rhoTag,
                  const Expr::Tag& mmwTag,
                  const double     mw_i,
                  const double     mw_n )
  : ExpressionBuilder( resultTags ),
    phiTags_( phiTags     ),
    rhoTag_ ( rhoTag      ),
    mmwTag_ ( mmwTag      ),
    mw_i_   ( mw_i        ),
    mw_n_   ( mw_n        )
{
  if( !(resultTags.size()==phiTags.size()) )
  {
    std::ostringstream msg;
    msg << __FILE__ << " : " << __LINE__
        << std::endl
        << "resultTags and phiTags must have the same number of elements." << std::endl
        <<"number of resultTags: "<< resultTags.size() << std::endl
        <<"number of phiTags   : "<< phiTags.size()    << std::endl;
    throw std::runtime_error( msg.str() );
  }
}

//====================================================================

// ###################################################################
//
//                          Implementation
//
// ###################################################################

template<typename FieldT>
PartialJacobian_Enthalpy<FieldT>::
PartialJacobian_Enthalpy( const Expr::TagList& phiTags,
                          const Expr::Tag&     rhoTag,
                          const Expr::Tag&     cpTag,
                          const Expr::Tag&     temperatureTag )
: Expr::Expression<FieldT>()
{
  this->set_gpu_runnable(true);

  this->template create_field_vector_request<FieldT>( phiTags, phi_ );

  rho_ = this->template create_field_request<FieldT>( rhoTag         );
  cp_  = this->template create_field_request<FieldT>( cpTag          );
  t_   = this->template create_field_request<FieldT>( temperatureTag );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
PartialJacobian_Enthalpy<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  typename Expr::Expression<FieldT>::ValVec&  result = this->get_value_vec();

  const FieldT& rho = rho_->field_ref();
  const FieldT& cp  = cp_ ->field_ref();
  const FieldT& t   = t_  ->field_ref();

  // calculate \f$x = -\frac{\partial \rho}{\partial h}\f$
  SpatFldPtr<FieldT> x = SpatialFieldStore::get<FieldT>( rho );
  *x <<= rho / ( t*cp );;

  /*
   * calculate \f$A_{i,j} = - \phi_i\frac{\partial \rho}{\partial h}\f$
   * for fixed \f i \f$.
   */

  for(int i=0; i<result.size(); ++i){
    *result[i] <<= ( phi_[i]->field_ref() )*(*x);
  }
}

//--------------------------------------------------------------------

template<typename FieldT>
PartialJacobian_Enthalpy<FieldT>::
Builder::Builder( const Expr::TagList& resultTags,
                  const Expr::TagList& phiTags,
                  const Expr::Tag&     rhoTag,
                  const Expr::Tag&     cpTag,
                  const Expr::Tag&     temperatureTag )
  : ExpressionBuilder( resultTags ),
    phiTags_       ( phiTags        ),
    rhoTag_        ( rhoTag         ),
    cpTag_         ( cpTag          ),
    temperatureTag_( temperatureTag )
{
  if( !(resultTags.size()==phiTags.size()) )
  {
    std::ostringstream msg;
    msg << __FILE__ << " : " << __LINE__
        << std::endl
        << "resultTags and phiTags must have the same number of elements."
        << std::endl;
    throw std::runtime_error( msg.str() );
  }
}

//====================================================================

// ###################################################################
//
//                          Implementation
//
// ###################################################################

template<typename FieldT>
PartialJacobian_Temperature<FieldT>::
PartialJacobian_Temperature( const Expr::TagList& phiTags,
                          const Expr::Tag&     rhoTag,
                          const Expr::Tag&     temperatureTag )
: Expr::Expression<FieldT>()
{
  this->set_gpu_runnable(true);

  this->template create_field_vector_request<FieldT>( phiTags, phi_ );

  rho_ = this->template create_field_request<FieldT>( rhoTag         );
  t_   = this->template create_field_request<FieldT>( temperatureTag );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
PartialJacobian_Temperature<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  typename Expr::Expression<FieldT>::ValVec&  result = this->get_value_vec();

  const FieldT& rho = rho_->field_ref();
  const FieldT& t   = t_  ->field_ref();

  // calculate \f$x = -\frac{\partial \rho}{\partial T}\f$
  SpatFldPtr<FieldT> x = SpatialFieldStore::get<FieldT>( rho );
  *x <<= rho / ( t );

  /*
   * calculate \f$A_{i,j} = - \phi_i\frac{\partial \rho}{\partial T}\f$
   */
  for(int i=0; i<result.size(); ++i){
    *result[i] <<= ( phi_[i]->field_ref() )*(*x);
  }
}

//--------------------------------------------------------------------

template<typename FieldT>
PartialJacobian_Temperature<FieldT>::
Builder::Builder( const Expr::TagList& resultTags,
                  const Expr::TagList& phiTags,
                  const Expr::Tag&     rhoTag,
                  const Expr::Tag&     temperatureTag )
  : ExpressionBuilder( resultTags ),
    phiTags_       ( phiTags        ),
    rhoTag_        ( rhoTag         ),
    temperatureTag_( temperatureTag )
{
  if( !(resultTags.size()==phiTags.size()) )
  {
    std::ostringstream msg;
    msg << __FILE__ << " : " << __LINE__
        << std::endl
        << "resultTags and phiTags must have the same number of elements."
        << std::endl;
    throw std::runtime_error( msg.str() );
  }
}

// explicit template instantiation
#include <spatialops/structured/FVStaggeredFieldTypes.h>
template class PartialJacobian_Species    <SpatialOps::SVolField>;
template class PartialJacobian_Enthalpy   <SpatialOps::SVolField>;
template class PartialJacobian_Temperature<SpatialOps::SVolField>;
