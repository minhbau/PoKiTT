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

#ifndef ThermalConductivity_Expr_h
#define ThermalConductivity_Expr_h

#include <expression/Expression.h>

#include <pokitt/CanteraObjects.h> //include cantera wrapper

namespace pokitt{

/**
 *  \class  ThermalConductivity
 *  \author Nate Yonkee
 *  \date July, 2014
 *
 *  \brief Calculates the thermal conductivity using a mixing rule [W/m/K].
 *
 *  This class calculates the thermal conductivity using the following
 *  mixing rule,
 *
 * \f[
 * \lambda = 0.5 \left( \sum_n x_n \lambda_n
 * + \frac{1}{\sum_n x_n/\lambda_n}\right)
 * \f]
 *
 * Where \f$ \lambda_n \f$ is the thermal conductivity of pure n and
 * \f$ x_n \f$ is the mole fraction of species n.
 *
 * Units are W/m/K.
 *
 */

template< typename FieldT >
class ThermalConductivity
    : public Expr::Expression<FieldT>
{
  DECLARE_FIELDS( FieldT, temperature_, mmw_ )
  DECLARE_VECTOR_OF_FIELDS( FieldT, massFracs_ )

  int nSpec_; //number of species to iterate over
  std::vector< std::vector<double> > tCondCoefs_;
  int modelType_; // type of model used by Cantera to estimate pure viscosity
  std::vector<double> molecularWeights_; // molecular weights
  std::vector<double> molecularWeightsInv_; // inverse of molecular weights (diving by MW is expensive)

  ThermalConductivity( const Expr::Tag& temperatureTag,
                       const Expr::TagList& massFracTags,
                       const Expr::Tag& mmwTag );
public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a ThermalConductivity expression
     *  @param resultTag the tag for the mixture averaged thermal conductivity
     *  @param temperatureTag temperature
     *  @param massFracTag tag for mass fraction of each species, ordering must be consistent with Cantera input
     *  @param mmwTag tag for mixture molecular weight
     *  @param nghost the number of ghost cells to compute in
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::Tag& temperatureTag,
             const Expr::TagList& massFracTags,
             const Expr::Tag& mmwTag,
             const int nghost = DEFAULT_NUMBER_OF_GHOSTS );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag temperatureTag_;
    const Expr::TagList massFracTags_;
    const Expr::Tag mmwTag_;
  };

  ~ThermalConductivity();
  void evaluate();
};



// ###################################################################
//
//                          Implementation
//
// ###################################################################



template< typename FieldT >
ThermalConductivity<FieldT>::
ThermalConductivity( const Expr::Tag& temperatureTag,
                     const Expr::TagList& massFracTags,
                     const Expr::Tag& mmwTag )
  : Expr::Expression<FieldT>()
{
  this->set_gpu_runnable( true );

  temperature_ = this->template create_field_request<FieldT>( temperatureTag );
  mmw_         = this->template create_field_request<FieldT>( mmwTag         );

  this->template create_field_vector_request<FieldT>( massFracTags, massFracs_ );

  Cantera::MixTransport* trans = dynamic_cast<Cantera::MixTransport*>( CanteraObjects::get_transport() ); // cast gas transport object as mix transport
  nSpec_ = trans->thermo().nSpecies();

  tCondCoefs_ = trans->getConductivityCoefficients();
  modelType_ = trans->model();

  molecularWeights_ = trans->thermo().molecularWeights();
  molecularWeightsInv_.resize(nSpec_);
  for( size_t n=0; n<nSpec_; ++n)
    molecularWeightsInv_[n] = 1 / molecularWeights_[n];

  CanteraObjects::restore_transport( trans );
}

//--------------------------------------------------------------------

template< typename FieldT >
ThermalConductivity<FieldT>::~ThermalConductivity(){}

//--------------------------------------------------------------------

template< typename FieldT >
void
ThermalConductivity<FieldT>::
evaluate()
{
  using namespace SpatialOps;

  FieldT& mixTCond = this->value();

  const FieldT& temp = temperature_->field_ref();
  const FieldT& mmw  = mmw_        ->field_ref();

  // pre-compute powers of temperature used in polynomial evaluations
  SpatFldPtr<FieldT> sqrtTPtr;  // may be used later on
  SpatFldPtr<FieldT> logtPtr = SpatialFieldStore::get<FieldT>(temp); // log(t)

  FieldT& logt = *logtPtr;
  logt <<= log( temp );

  SpatFldPtr<FieldT> speciesTCondPtr = SpatialFieldStore::get<FieldT>(temp); // temporary to store the thermal conductivity for an individual species
  SpatFldPtr<FieldT> sumPtr          = SpatialFieldStore::get<FieldT>(temp); // for mixing rule
  SpatFldPtr<FieldT> inverseSumPtr   = SpatialFieldStore::get<FieldT>(temp); // 1/sum()

  FieldT& speciesTCond = *speciesTCondPtr;
  FieldT& sum          = *sumPtr;
  FieldT& inverseSum   = *inverseSumPtr;

  sum        <<= 0.0; // set sum to 0 before loop
  inverseSum <<= 0.0; // set inverse sum to 0 before loop

  if( modelType_ == Cantera::cMixtureAveraged ) { // as opposed to CK mode
    sqrtTPtr = SpatialFieldStore::get<FieldT>(temp);
    *sqrtTPtr <<= sqrt( temp );
  }

  for( size_t n = 0; n < nSpec_; ++n){
    const FieldT& yi = massFracs_[n]->field_ref();
    const std::vector<double>& tCondCoefs = tCondCoefs_[n];
    if( modelType_ == Cantera::cMixtureAveraged )
      speciesTCond <<= *sqrtTPtr * ( tCondCoefs[0] + logt * ( tCondCoefs[1] + logt * ( tCondCoefs[2] + logt * ( tCondCoefs[3] + logt * tCondCoefs[4] ))) );
    else
      speciesTCond <<=         exp ( tCondCoefs[0] + logt * ( tCondCoefs[1] + logt * ( tCondCoefs[2] + logt *   tCondCoefs[3]                         )) );
    sum        <<= sum        + yi * speciesTCond * molecularWeightsInv_[n];
    inverseSum <<= inverseSum + yi / speciesTCond * molecularWeightsInv_[n];
  }

  mixTCond <<= 0.5 * ( sum * mmw + 1 / (inverseSum * mmw) ); // mixing rule
}

//--------------------------------------------------------------------

template< typename FieldT >
ThermalConductivity<FieldT>::
Builder::Builder( const Expr::Tag& resultTag,
                  const Expr::Tag& temperatureTag,
                  const Expr::TagList& massFracTags,
                  const Expr::Tag& mmwTag,
                  const int nghost )
: ExpressionBuilder( resultTag, nghost ),
  temperatureTag_( temperatureTag ),
  massFracTags_( massFracTags ),
  mmwTag_( mmwTag )
{}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
ThermalConductivity<FieldT>::
Builder::build() const
{
  return new ThermalConductivity<FieldT>( temperatureTag_, massFracTags_, mmwTag_ );
}

} // namespace pokitt

#endif // ThermalConductivity_Expr_h
