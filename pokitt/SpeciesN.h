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

#ifndef SpeciesN_Expr_h
#define SpeciesN_Expr_h

#include <expression/Expression.h>

namespace pokitt{

/**
 *  \class SpeciesN
 *  \brief Sets the requested species to enforce the constraint that \f$\sum y_i=1\f$.
 *
 *  Note that if "species N" must be greater than one or less than zero to enforce
 *  a sum to unity, then an exception will be thrown.
 */
template< typename FieldT >
class SpeciesN : public Expr::Expression<FieldT>
{
  DECLARE_VECTOR_OF_FIELDS( FieldT, species_ )

  SpeciesN( const Expr::TagList& speciesTags ) : Expr::Expression<FieldT>()
  {
    this->set_gpu_runnable(true);
    this->template create_field_vector_request<FieldT>( speciesTags, species_ );
  }

public:

  class Builder : public Expr::ExpressionBuilder
  {
    Expr::TagList speciesTags_;
  public:
    /**
     *  @brief Build a SpeciesN expression
     *  @param specNTag the tag for species "N" - the one to compute here
     *  @param speciesTags the list of tags for all species
     *  @param nghost number of ghost cells
     */
    Builder( const Expr::Tag& specNTag,
             const Expr::TagList& speciesTags,
             const int nghost = DEFAULT_NUMBER_OF_GHOSTS )
      : ExpressionBuilder( specNTag, nghost )
    {
      for( size_t i=0; i<speciesTags.size(); ++i ){
        if( speciesTags[i] != specNTag ) speciesTags_.push_back( speciesTags[i] );
      }
    }

    /**
     *  @brief Build a SpeciesN expression
     *  @param speciesTags the list of tags for all species
     *  @param specNum the species index (0-based) for the species to absorb the constraint that \f$\sum y_i = 1\f$.
     *  @param nghost number of ghost cells
     */
    Builder( const Expr::TagList& speciesTags,
             const size_t specNum,
             const int nghost = DEFAULT_NUMBER_OF_GHOSTS )
      : ExpressionBuilder( speciesTags[specNum], nghost ),
        speciesTags_( speciesTags )
    {
      speciesTags_.erase( speciesTags_.begin() + specNum );
    }

    Expr::ExpressionBase* build() const{
      return new SpeciesN<FieldT>( speciesTags_ );
    }
  };

  void evaluate()
  {
    FieldT& specN = this->value();

    specN <<= 1.0 - species_[0]->field_ref();
    for( size_t i=1; i<species_.size(); ++i ){
      specN <<= specN - species_[i]->field_ref();
    }

#   ifndef ENABLE_CUDA
    if( nebo_max( specN ) > 1.0 || nebo_min( specN ) < 0.0 ){
      std::ostringstream msg;
      msg << __FILE__ << " : " << __LINE__ << "\nSpecies fractions are out of bounds!\n";
      throw std::runtime_error( msg.str() );
    }
#   endif
  }

};

} // namespace pokitt

#endif // SpeciesN_Expr_h
