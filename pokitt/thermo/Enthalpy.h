#ifndef SpeciesEnthalpy_Expr_h
#define SpeciesEnthalpy_Expr_h

//#define TIMINGS

#include <expression/Expression.h>
#ifdef TIMINGS
#include <boost/timer.hpp>
#endif

#include <cantera/kernel/ct_defs.h> // contains value of gas constant
#include <cantera/kernel/SpeciesThermoInterpType.h> // contains definitions for which polynomial is being used
#include <cantera/IdealGasMix.h>

#include <pokitt/CanteraObjects.h> //include cantera wrapper

namespace Cantera_CXX{ class IdealGasMix; } //location of polynomial

/**
 *  \class  SpeciesEnthalpy
 *  \author Nate Yonkee
 *  \date July, 2014
 *
 *  \brief Calculates the enthalpy of each species
 *  using either NASA7 or Shomate polynomials with 2 temperature ranges.
 *  Units of J/kg
 *
 *  Five coefficients \f$(a_0,\dots,a_4)\f$ are used to represent
 * \f$ h(T)\f$ as a polynomial in \f$ T \f$ e.g. NASA7:
 *
 * \f[
 * \frac{h(T)}{R} = a_0 T + a_1 T^2/2 + a_2 T^3/3 + a_3 T^4/4 + a_4 T^5/5 + a_5
 * \f]
 *
 */
template< typename FieldT >
class SpeciesEnthalpy
 : public Expr::Expression<FieldT>
{
  typedef std::vector<FieldT*> SpecT;
  const Expr::Tag tTag_;
  const FieldT* t_;
  Cantera_CXX::IdealGasMix* const gasMix_;
  const int nSpec_; //number of species
  bool nasaFlag_; //flag to specify if any NASA polynomials are present
  bool shomateFlag_; //flag to specify if any Shomate polynomials are present


  SpeciesEnthalpy( const Expr::Tag& tTag );
public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a SpeciesEnthalpy expression
     *  @param resultTags the list of tags for each species' mass enthalpy, order must be consistent with Cantera's input file
     *  @param tTag tag for temperature
     */
    Builder( const Expr::TagList& resultTags,
             const Expr::Tag& tTag );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag tTag_;
  };

  ~SpeciesEnthalpy();
  void advertise_dependents( Expr::ExprDeps& exprDeps );
  void bind_fields( const Expr::FieldManagerList& fml );
  void evaluate();
};



// ###################################################################
//
//                          Implementation
//
// ###################################################################



template< typename FieldT >
SpeciesEnthalpy<FieldT>::
SpeciesEnthalpy( const Expr::Tag& tTag )
  : Expr::Expression<FieldT>(),
    tTag_( tTag ),
    gasMix_( CanteraObjects::get_gasmix() ),
    nSpec_( gasMix_->nSpecies() )
{
  this->set_gpu_runnable( true );

  int polyType;
  const Cantera::SpeciesThermo& spThermo = gasMix_->speciesThermo();

  for( size_t n=0; n<nSpec_; ++n){
    polyType = spThermo.reportType(n);
    if( polyType == NASA2 )
      nasaFlag_ = true;
    else if( polyType == SHOMATE2 )
      shomateFlag_ = true;
    else if( polyType != SIMPLE ){
      std::cout << "Error in SpeciesEnthalpy" << std::endl
          <<" Thermo type not supported, type = " << polyType << std::endl
          <<" See <cantera/kernel/SpeciesThermoInterpType.h> for type definition" << std::endl;
      throw std::runtime_error("Problems in SpeciesEnthalpy expression.");
    }
  }
}

//--------------------------------------------------------------------

template< typename FieldT >
SpeciesEnthalpy<FieldT>::
~SpeciesEnthalpy()
{
  CanteraObjects::restore_gasmix(gasMix_);
}

//--------------------------------------------------------------------

template< typename FieldT >
void
SpeciesEnthalpy<FieldT>::
advertise_dependents( Expr::ExprDeps& exprDeps )
{
  exprDeps.requires_expression( tTag_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
SpeciesEnthalpy<FieldT>::
bind_fields( const Expr::FieldManagerList& fml )
{
  t_ = &fml.template field_ref< FieldT >( tTag_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
SpeciesEnthalpy<FieldT>::
evaluate()
{
#ifdef TIMINGS
  boost::timer timer;
#endif
  using namespace SpatialOps;
  using namespace Cantera;
  SpecT& enthalpies = this->get_value_vec();

  const FieldT& temp = *t_;

  SpatFldPtr<FieldT> t2;
  SpatFldPtr<FieldT> t3;
  SpatFldPtr<FieldT> t4;
  SpatFldPtr<FieldT> t5;
  SpatFldPtr<FieldT> recipT;

  // pre-compute powers of t for polynomial evaluations
  if( nasaFlag_ == true || shomateFlag_ == true){
    t2 = SpatialFieldStore::get<FieldT>(*t_); // t^2
    t3 = SpatialFieldStore::get<FieldT>(*t_); // t^3
    t4 = SpatialFieldStore::get<FieldT>(*t_); // t^4

    *t2 <<= temp * temp;
    *t3 <<= *t2 * temp;
    *t4 <<= *t3 * temp;
    if( nasaFlag_ == true ){
      t5 = SpatialFieldStore::get<FieldT>(*t_); // t^5
      *t5 <<= *t4 * temp;
    }
    if( shomateFlag_ == true ){
      recipT = SpatialFieldStore::get<FieldT>(*t_); // = 1/t
      *recipT <<= 1 / temp;
    }
  }

  std::vector<double> c(15,0); //vector of Cantera's coefficients
  int polyType; // type of polynomial
  double minT; // minimum temperature where polynomial is valid
  double maxT; // maximum temperature where polynomial is valid
  double refPressure;
  const Cantera::SpeciesThermo& spThermo = gasMix_->speciesThermo();

  const std::vector<double>& molecularWeights = gasMix_->molecularWeights();

  for( size_t n=0; n<nSpec_; ++n ){
    spThermo.reportParams(n, polyType, &c[0], minT, maxT, refPressure);
    std::vector<double>::iterator ic = c.begin() + 1;
    std::vector<double>::iterator icend = c.end();
    switch (polyType) {
    case SIMPLE: // constant heat capacity
      *enthalpies[n] <<= ( c[1] + c[3] * (temp-c[0]) )
                         / molecularWeights[n];
      break;
    case NASA2:
      /* polynomials are applicable in two temperature ranges - high and low
       * If the temperature is out of range, the value is set to the value at the min or max temp
       */
      for( ; ic != icend; ++ic)
        *ic *= GasConstant / molecularWeights[n]; // dimensionalize the coefficients
      *enthalpies[n] <<= cond( temp <= c[0] && temp >= minT, c[ 6] + c[1] * temp + c[2]/2 * *t2 + c[ 3]/3 * *t3 + c[ 4]/4 * *t4 + c[ 5]/5 * *t5) // if low temp
                             ( temp >  c[0] && temp <= maxT, c[13] + c[8] * temp + c[9]/2 * *t2 + c[10]/3 * *t3 + c[11]/4 * *t4 + c[12]/5 * *t5)  // else if high temp
                             ( temp <  minT                 , c[1] * temp + c[2] * minT * (temp-minT/2) + c[ 3] * minT * minT * (temp - 2*minT/3) + c[ 4] * pow(minT,3) * (temp - 3*minT/4) + c[ 5] * pow(minT,4) * (temp - 4*minT/5) + c[ 6]) // else if out of bounds - low
                             (                                c[8] * temp + c[9] * maxT * (temp-maxT/2) + c[10] * maxT * maxT * (temp - 2*maxT/3) + c[11] * pow(maxT,3) * (temp - 3*maxT/4) + c[12] * pow(maxT,4) * (temp - 4*maxT/5) + c[13]); // else out of bounds - high
    break;
   case SHOMATE2:
     for( ; ic != icend; ++ic)
        *ic *= 1e6 / molecularWeights[n]; // scale the coefficients to keep units consistent
      double minTScaled = minT/1000;
      double maxTScaled = maxT/1000;
      *enthalpies[n] <<= cond( temp <= c[0] && temp >= minT, c[ 6] + c[1] * temp*1e-3 + c[2]/2 * *t2*1e-6 + c[ 3]/3 * *t3*1e-9 + c[ 4]/4 * *t4*1e-12 - c[ 5] * *recipT*1e3 ) // if low temp
                             ( temp >  c[0] && temp <= maxT, c[13] + c[8] * temp*1e-3 + c[9]/2 * *t2*1e-6 + c[10]/3 * *t3*1e-9 + c[11]/4 * *t4*1e-12 - c[12] * *recipT*1e3 )  // else if high range
                             ( temp <  minT, c[1] * temp*1e-3 + c[2] * minTScaled * ( temp*1e-3 - minTScaled/2 ) + c[ 3] * minTScaled * minTScaled * ( temp*1e-3 - 2*minTScaled/3 ) + c[ 4] * pow(minTScaled,3) * ( temp*1e-3 - 3*minTScaled/4 ) - c[ 5] * pow(minTScaled,-1) * ( -temp*1e-3/minTScaled + 2 ) + c[ 6] ) // else if out of bounds - low
                             (               c[8] * temp*1e-3 + c[9] * maxTScaled * ( temp*1e-3 - maxTScaled/2 ) + c[10] * maxTScaled * maxTScaled * ( temp*1e-3 - 2*maxTScaled/3 ) + c[11] * pow(maxTScaled,3) * ( temp*1e-3 - 3*maxTScaled/4 ) - c[12] * pow(maxTScaled,-1) * ( -temp*1e-3/maxTScaled + 2 ) + c[13] ); // else out of bounds - high
    }
  }
#ifdef TIMINGS
    std::cout<<"h time "<<timer.elapsed()<<std::endl;
#endif
}

//--------------------------------------------------------------------

template< typename FieldT >
SpeciesEnthalpy<FieldT>::
Builder::Builder( const Expr::TagList& resultTags,
                  const Expr::Tag& tTag )
  : ExpressionBuilder( resultTags ),
    tTag_( tTag )
{}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
SpeciesEnthalpy<FieldT>::
Builder::build() const
{
  return new SpeciesEnthalpy<FieldT>( tTag_ );
}

#endif // SpeciesEnthalpy_Expr_h
