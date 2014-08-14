#ifndef MixtureMolWeightExpr_h
#define MixtureMolWeightExpr_h

#include <expression/ExprLib.h>
        // defines field types
typedef SpatialOps::SVolField   VolField;
//====================================================================

/**
 *  @class MixtureMolWeightExpr
 *  @date January, 2011
 *  @author Naveen Punati
 *  @brief Caculates mixture molecular weight given the invidual species mass fractions and molecular weights
 *         Used for Initial conditions and BC's
 */
class MixtureMolWeightExpr : public Expr::Expression< VolField >
{
public:
  typedef std::vector<const VolField*> ConstSpecT;
  struct Builder : public Expr::ExpressionBuilder
  {
    Builder( const Expr::Tag& result,
             const Expr::Tag& yiTag,
             const std::vector<double>& specMw )
      : ExpressionBuilder(result),
        yiT_( yiTag ),
        specMw_(specMw)
    {}

    Expr::ExpressionBase* build() const{ return new MixtureMolWeightExpr(yiT_,specMw_ ); }

  private:
    const Expr::Tag yiT_;
    const std::vector<double> specMw_;

  };

  void advertise_dependents( Expr::ExprDeps& exprDeps );
  void bind_fields( const Expr::FieldManagerList& fml );
  void evaluate();

private:

  MixtureMolWeightExpr( const Expr::Tag& yiTag,
                        const std::vector<double>& specMw );

  ~MixtureMolWeightExpr();

  const std::vector<double> specMw_;
  const int nspec_;
  Expr::TagList specTags_;

  ConstSpecT yi_;
};

//====================================================================

MixtureMolWeightExpr::
MixtureMolWeightExpr( const Expr::Tag& yiTag,
                      const std::vector<double>& specMw )
  : Expr::Expression<VolField>(),
    specMw_(specMw),
    nspec_(specMw_.size())
{
  specTags_.clear();
  for( int i=0; i<nspec_; ++i ){
    std::ostringstream name;
    name << yiTag.name() << "_" << i;
    specTags_.push_back( Expr::Tag(name.str(),yiTag.context()) );
  }
}

//--------------------------------------------------------------------

MixtureMolWeightExpr::~MixtureMolWeightExpr()
{}

//--------------------------------------------------------------------

void
MixtureMolWeightExpr::
advertise_dependents( Expr::ExprDeps& exprDeps )
{
  for( Expr::TagList::const_iterator itag=specTags_.begin(); itag!=specTags_.end(); ++itag ){
    exprDeps.requires_expression( *itag );
  }
}

//--------------------------------------------------------------------

void
MixtureMolWeightExpr::
bind_fields( const Expr::FieldManagerList& fml )
{
  const Expr::FieldMgrSelector<VolField>::type& fm = fml.field_manager<VolField>();

  yi_.clear();
  for( Expr::TagList::const_iterator itag=specTags_.begin(); itag!=specTags_.end(); ++itag ){
    yi_.push_back( &fm.field_ref( *itag ) );
  }
}

//--------------------------------------------------------------------

void
MixtureMolWeightExpr::
evaluate()
{
  using namespace SpatialOps;
  VolField& result = this->value();

  /// \todo implement in till to avoid double grid loop...
  result  <<= 0.0;
  for( size_t i=0; i<nspec_; ++i ){
    result <<= result + *yi_[i] / specMw_[i];
  }
  result <<= 1.0 / result;
}

#endif // MixtureMolWeightExpr_h
