#ifndef HeatCapacity_Cp_Expr_h
#define HeatCapacity_Cp_Expr_h

#include <expression/Expression.h>

#include <boost/timer.hpp>

#include <pokitt/CanteraObjects.h> //include cantera wrapper

#include <cantera/kernel/ct_defs.h> // contains value of gas constant
#include <cantera/kernel/speciesThermoTypes.h> // contains definitions for which polynomial is being used
#include <cantera/IdealGasMix.h>

namespace Cantera_CXX{ class IdealGasMix; } // location of polynomial coefficients

/**
 *  \class  HeatCapacity_Cp
 *  \author Nate Yonkee
 *  \date May, 2014
 *
 *  \brief Calculates the constant pressure heat capacity of each species
 *  using either NASA7 or Shomate polynomials with 2 temperature ranges.
 *  Units of J/kmole/K
 *
 *  Five coefficients \f$(a_0,\dots,a_4)\f$ are used to represent
 * \f$ c_p^0(T)\f$ as a polynomial in \f$ T \f$
 *
 * NASA7:
 *
 * \f[
 * \frac{c_p(T)}{R} = a_0 + a_1 T + a_2 T^2 + a_3 T^3 + a_4 T^4
 * \f]
 *
 * Shomate:
 *
 * \f[
 * \frac{c_p(T)}{1000} = a_0 + a_1 t + a_2 t^2 + a_3 t^3 + \frac{a_4}{t^2}
 * \f]
 *
 * Where \f[ t = T (K)/1000 \f]
 *
 */
template< typename FieldT >
class HeatCapacity_Cp
 : public Expr::Expression<FieldT>
{
  typedef std::vector<FieldT*> SpecT;
  typedef std::vector<double>  SpecCoefs;
  const Expr::Tag tTag_;
  const FieldT* t_;
  Cantera_CXX::IdealGasMix* const gasMix_;
  const int nSpec_; // number of species
  bool nasaFlag_; // flag to specify if any NASA polynomials are present
  bool shomateFlag_; // flag to specify if any Shomate polynomials are present

  HeatCapacity_Cp( const Expr::Tag& tTag );
public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a HeatCapacity_Cp expression
     *  @param resultTags the list of tags for each species' heat capacity, order must be consistent with Cantera's input file
     *  @param tTag tag for temperature
     */
    Builder( const Expr::TagList& resultTags,
             const Expr::Tag& tTag );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag tTag_;
  };

  ~HeatCapacity_Cp();
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
HeatCapacity_Cp<FieldT>::
HeatCapacity_Cp( const Expr::Tag& tTag )
  : Expr::Expression<FieldT>(),
    tTag_( tTag ),
    gasMix_( CanteraObjects::get_gasmix() ),
    nSpec_( gasMix_->nSpecies() ),
    nasaFlag_( false ),
    shomateFlag_( false )
{
  this->set_gpu_runnable( true );
  const Cantera::SpeciesThermo& spThermo = gasMix_->speciesThermo();

  for(size_t n=0; n<nSpec_; ++n){
    int type = spThermo.reportType(n);
    if( type == NASA2 )
      nasaFlag_ = true;
    else if( type == SHOMATE2 )
      shomateFlag_ = true;
    else if( type != SIMPLE )
      std::cout<<"Thermo type not supported,\nType = "<<type<<" species # "<<n<<std::endl;
  }
}

//--------------------------------------------------------------------

template< typename FieldT >
HeatCapacity_Cp<FieldT>::
~HeatCapacity_Cp()
{
  CanteraObjects::restore_gasmix(gasMix_);
}

//--------------------------------------------------------------------

template< typename FieldT >
void
HeatCapacity_Cp<FieldT>::
advertise_dependents( Expr::ExprDeps& exprDeps )
{
  exprDeps.requires_expression( tTag_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
HeatCapacity_Cp<FieldT>::
bind_fields( const Expr::FieldManagerList& fml )
{
  t_ = &fml.template field_ref< FieldT >( tTag_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
HeatCapacity_Cp<FieldT>::
evaluate()
{
#ifdef TIMINGS
  boost::timer time;
#endif
  using namespace SpatialOps;
  using namespace Cantera;
  SpecT& cps = this->get_value_vec();

  SpatFldPtr<FieldT> t2;
  SpatFldPtr<FieldT> t3;
  SpatFldPtr<FieldT> t4;
  SpatFldPtr<FieldT> recipRecipT;

  // pre-compute powers of t for polynomial evaluations
    if( nasaFlag_ == true | shomateFlag_ == true){
      t2 = SpatialFieldStore::get<FieldT>(*t_);
      t3 = SpatialFieldStore::get<FieldT>(*t_);

      *t2 <<= *t_ * *t_; // t^2
      *t3 <<= *t2 * *t_; // t^3

      if( nasaFlag_ == true ){
        t4 = SpatialFieldStore::get<FieldT>(*t_);
        *t4 <<= *t3 * *t_; // t^4
      }

      if( shomateFlag_ == true ){
        recipRecipT = SpatialFieldStore::get<FieldT>(*t_);
        *recipRecipT <<= 1/ *t2; // = t^-2
      }
    }

  std::vector<double> c(15,0); //vector of Cantera's coefficients
  int polyType;
  double minTemp;
  double maxTemp;
  double refPressure;
  const Cantera::SpeciesThermo& spThermo = gasMix_->speciesThermo();

  for( size_t n=0; n<nSpec_; ++n ){
    spThermo.reportParams(n, polyType, &c[0], minTemp, maxTemp, refPressure);
    if (polyType==CONSTANT_CP | polyType==SIMPLE)
      *cps[n] <<= c[3]; // constant heat capacity
    /* polynomials are applicable in two temperature ranges - high and low
     * If the temperature is out of range, the value is set to the value at the min or max temp
     */
    else if (polyType==NASA2){
      for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic )
        *ic *= GasConstant; // dimensionalize the coefficients
      *cps[n] <<= cond( *t_ <= c[0] && *t_ >= minTemp, c[1] + c[2] * *t_ + c[ 3] * *t2 + c[ 4] * *t3 + c[ 5] * *t4) // if low temp
                      ( *t_ >  c[0] && *t_ <= maxTemp, c[8] + c[9] * *t_ + c[10] * *t2 + c[11] * *t3 + c[12] * *t4)  // else if high temp
                      ( *t_ < minTemp, c[1] + c[2] * minTemp + c[ 3] * minTemp * minTemp + c[ 4] * pow(minTemp,3) + c[ 5] * pow(minTemp,4))  // else if out of bounds - low
                      (                c[8] + c[9] * maxTemp + c[10] * maxTemp * maxTemp + c[11] * pow(maxTemp,3) + c[12] * pow(maxTemp,4)); // else out of bounds - high
    }
    else if (polyType==SHOMATE2){
      for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic )
        *ic *= 1e3; // scale the coefficients
      double minTempScaled = minTemp/1000;
      double maxTempScaled = maxTemp/1000;
      *cps[n] <<= cond( *t_ <= c[0] && *t_ >= minTemp, c[1] + c[2] * *t_*1e-3 + c[ 3] * *t2*1e-6 + c[ 4] * *t3*1e-9 + c[ 5] * *recipRecipT*1e6) // if low temp
                      ( *t_ >  c[0] && *t_ <= maxTemp, c[8] + c[9] * *t_*1e-3 + c[10] * *t2*1e-6 + c[11] * *t3*1e-9 + c[12] * *recipRecipT*1e6)  // else if high temp
                      ( *t_ < minTemp, c[1] + c[2] * minTempScaled + c[ 3] * minTempScaled * minTempScaled + c[ 4] * pow(minTempScaled,3) + c[ 5] * pow(minTempScaled,-2))  // else if out of bounds - low
                      (                c[8] + c[9] * maxTempScaled + c[10] * maxTempScaled * maxTempScaled + c[11] * pow(maxTempScaled,3) + c[12] * pow(maxTempScaled,-2)); // else out of bounds - high
    }
  }
# ifdef TIMINGS
  std::cout << "HeatCapacity_Cp time " << timer.elapsed() << std::endl;
# endif
}

//--------------------------------------------------------------------

template< typename FieldT >
HeatCapacity_Cp<FieldT>::
Builder::Builder( const Expr::TagList& resultTags,
                  const Expr::Tag& tTag )
  : ExpressionBuilder( resultTags ),
    tTag_( tTag )
{}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
HeatCapacity_Cp<FieldT>::
Builder::build() const
{
  return new HeatCapacity_Cp<FieldT>( tTag_ );
}

#endif // HeatCapacity_Cp_Expr_h
