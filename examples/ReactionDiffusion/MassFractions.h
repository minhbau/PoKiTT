#ifndef MassFractions_Expr_h
#define MassFractions_Expr_h

#include <expression/Expression.h>

namespace pokitt{

/**
 *  \class  MassFractions
 *  \author Nathan Yonkee
 *  \date February, 2015
 *
 *  \brief Calculates primitive variables, mass fractions,
 *  from transported variables, species mass.
 */

template< typename FieldT >
class MassFractions
    : public Expr::Expression<FieldT>
{
  DECLARE_FIELDS( FieldT, rho_ )
  DECLARE_VECTOR_OF_FIELDS( FieldT, rhoYi_ )

  const int nSpec_; //number of species to iterate over

  MassFractions( const Expr::Tag& rhoTag,
                 const Expr::TagList& rhoYiTags);
public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a MassFractions expression
     *  @param yiTags the tag for the mass fractions
     *  @param rhoTag density
     *  @param rhoYiTags tags for transported mass, ordering must be consistent with Cantera input
     *  @param nghost the number of ghost cells to compute in
     */
    Builder( const Expr::TagList& yiTags,
             const Expr::Tag& rhoTag,
             const Expr::TagList& rhoYiTags,
             const int nghost = DEFAULT_NUMBER_OF_GHOSTS );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag rhoTag_;
    const Expr::TagList rhoYiTags_;
  };

  ~MassFractions();
  void evaluate();
};



// ###################################################################
//
//                          Implementation
//
// ###################################################################



template< typename FieldT >
MassFractions<FieldT>::
MassFractions( const Expr::Tag& rhoTag,
               const Expr::TagList& rhoYiTags )
  : Expr::Expression<FieldT>(),
    nSpec_( rhoYiTags.size()+1 )
{
  this->set_gpu_runnable( true );

  rho_ = this->template create_field_request<FieldT>( rhoTag );
  this->template create_field_vector_request<FieldT>( rhoYiTags, rhoYi_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
MassFractions<FieldT>::~MassFractions(){}

//--------------------------------------------------------------------

template< typename FieldT >
void
MassFractions<FieldT>::
evaluate()
{
  using namespace SpatialOps;

  std::vector< FieldT* >& yi = this->get_value_vec();

  const FieldT& rho = rho_->field_ref();

  *yi[nSpec_-1] <<= 1.0;
  for( size_t n = 0; n < (nSpec_-1); ++n ){
    const FieldT& rhoYi = rhoYi_[n]->field_ref();
    *yi[n] <<= rhoYi / rho;
    *yi[nSpec_-1] <<= *yi[nSpec_-1] - *yi[n];
  }

}

//--------------------------------------------------------------------

template< typename FieldT >
MassFractions<FieldT>::
Builder::Builder( const Expr::TagList& resultTags,
                  const Expr::Tag& rhoTag,
                  const Expr::TagList& rhoYiTags,
                  const int nghost )
: ExpressionBuilder( resultTags, nghost ),
  rhoTag_( rhoTag ),
  rhoYiTags_( rhoYiTags )
{}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
MassFractions<FieldT>::
Builder::build() const
{
  return new MassFractions<FieldT>( rhoTag_, rhoYiTags_ );
}

} // namespace pokitt

#endif // MassFractions_Expr_h
