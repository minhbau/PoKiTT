/*
 * The MIT License
 *
 * Copyright (c) 2015 The University of Utah
 *
  * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifndef MixtureFractionExpr_h
#define MixtureFractionExpr_h

#include <expression/ExprLib.h>

#include <pokitt/MixtureFraction.h>  // mixture fraction support


namespace pokitt{


/**
 * \class SpeciesToMixtureFraction
 * \author James C. Sutherland
 * \brief Calculate the mixture fraction given the species compositions.
 */
template<typename FieldT>
class SpeciesToMixtureFraction
  : public Expr::Expression<FieldT>
{
  DECLARE_VECTOR_OF_FIELDS( FieldT, species_ )

  const MixtureFraction mixfrac_;

  SpeciesToMixtureFraction( const Expr::TagList& specTags,
                       const std::vector<double>& fuelMassFracs,
                       const std::vector<double>& oxidMassFracs )
  : Expr::Expression<FieldT>(),
    mixfrac_( oxidMassFracs, fuelMassFracs, true )
  {
    this->set_gpu_runnable(true);
    this->template create_field_vector_request<FieldT>( specTags, species_ );
  }

  ~SpeciesToMixtureFraction(){}

public:

  class Builder : public Expr::ExpressionBuilder
  {
    const Expr::TagList specTags_;
    const std::vector<double> fuel_, oxid_;
  public:
    Builder( const Expr::Tag& result,
             const Expr::TagList& specTags,
             const std::vector<double>& fuelMassFracs,
             const std::vector<double>& oxidMassFracs,
             const SpatialOps::GhostData nghost = DEFAULT_NUMBER_OF_GHOSTS )
    : ExpressionBuilder(result, nghost),
      specTags_( specTags ),
      fuel_( fuelMassFracs ),
      oxid_( oxidMassFracs )
    {}

    Expr::ExpressionBase* build() const{
      return new SpeciesToMixtureFraction<FieldT>( specTags_, fuel_, oxid_ );
    }
  };

  void evaluate(){
    mixfrac_.species_to_mixfrac( species_, this->value() );
  }

};

//--------------------------------------------------------------------


/**
 * \class MixtureFractionToSpecies
 * \author James C. Sutherland
 * \brief Calculate the species compositions (either unreacted or reacted with infinitely fast chemistry) given the mixture fraction.
 */
template<typename FieldT>
class MixtureFractionToSpecies
  : public Expr::Expression<FieldT>
{
  DECLARE_FIELD( FieldT, mixFracField_ )

  const MixtureFraction mixfrac_;
  const ReactedState state_;

  MixtureFractionToSpecies( const Expr::Tag& mixFracTag,
                            const std::vector<double>& fuelMassFracs,
                            const std::vector<double>& oxidMassFracs,
                            const ReactedState state )
  : Expr::Expression<FieldT>(),
    mixfrac_( oxidMassFracs, fuelMassFracs, true ),
    state_( state )
  {
    this->set_gpu_runnable(true);
    this->template create_field_request<FieldT>( mixFracTag, mixFracField_ );
  }

  ~MixtureFractionToSpecies(){}

public:

  enum ReactedState{
    REACTED_COMPOSITION,
    UNREACTED_COMPOSITION
  };

  class Builder : public Expr::ExpressionBuilder
  {
    const Expr::Tag mixFracTag_;
    const std::vector<double> fuel_, oxid_;
    const ReactedState state_;
  public:
    /**
     * @brief Construct a MixtureFractionToSpecies builder.
     * @param specTags the species mass fractions that we compute in this expression
     * @param mixFracTag the mixture fraction
     * @param fuelMassFracs fuel stream mass fractions
     * @param oxidMassFracs oxidizer stream mass fractions
     * @param state specifies whether you want the unreacted or reacted species compositions.
     *        REACTED_COMPOSITION gives the Burke-Schumann solution (infinitely fast chemistry).
     *        UNREACTED_COMPOSITION gives the unreacted simple mixing solution.
     * @param nghost the number of ghost cells
     */
    Builder( const Expr::TagList& specTags,
             const Expr::Tag& mixFracTag,
             const std::vector<double>& fuelMassFracs,
             const std::vector<double>& oxidMassFracs,
             const ReactedState state,
             const SpatialOps::GhostData nghost = DEFAULT_NUMBER_OF_GHOSTS )
    : Expr::ExpressionBuilder(specTags,nghost),
      mixFracTag_( mixFracTag ),
      fuel_( fuelMassFracs ),
      oxid_( oxidMassFracs ),
      state_( state )
    {}

    Expr::ExpressionBase* build() const{
      return new MixtureFractionToSpecies<FieldT>( mixFracTag_, fuel_, oxid_, state_ );
    }
  };

  void evaluate(){
    switch( state_ ){
      case UNREACTED_COMPOSITION: mixfrac_.mixfrac_to_species   ( mixFracField_, this->get_value_vec() ); break;
      case   REACTED_COMPOSITION: mixfrac_.estimate_product_comp( mixFracField_, this->get_value_vec() ); break;
      default: assert( false );
    }
  }

};

//--------------------------------------------------------------------
} // namespace pokitt

#endif // MixtureFractionExpr_h
