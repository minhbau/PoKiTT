/*
 * The MIT License
 *
 * Copyright (c) 2015-2017 The University of Utah
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

#ifndef ReactionRates_Expr_h
#define ReactionRates_Expr_h

#include <cmath>

#include <expression/Expression.h>

#include <spatialops/structured/MatVecFields.h>
#include <spatialops/structured/MatVecOps.h>

#include <pokitt/CanteraObjects.h> // include cantera wrapper

#include <pokitt/kinetics/ReactionInfo.h>
#include <pokitt/kinetics/AnalyticalJacobian.h>

// these are macros to expand out prop(0) + prop(1) + ...
// this is used to maximize kernel size for the gpu
#define POKITT_SUM2( prop ) (          prop(0) + prop(1))
#define POKITT_SUM3( prop ) (POKITT_SUM2(prop) + prop(2))
#define POKITT_SUM4( prop ) (POKITT_SUM3(prop) + prop(3))
#define POKITT_SUM5( prop ) (POKITT_SUM4(prop) + prop(4))
#define POKITT_SUM6( prop ) (POKITT_SUM5(prop) + prop(5))
#define POKITT_SUM7( prop ) (POKITT_SUM6(prop) + prop(6))
#define POKITT_SUM8( prop ) (POKITT_SUM7(prop) + prop(7))
#define POKITT_SUM9( prop ) (POKITT_SUM8(prop) + prop(8))

#define ARRHENIUS( coef ) ( (coef)[0] * exp( (coef)[1] * logT - (coef)[2] * tRecip ) ) // arrhenius kinetics

#define THD_BDY( i )      ( thdBodies [(i)].thdBdyEff *    yi_[thdBodies [(i)].index]->field_ref() ) // third body effects
#define GIBBS( i )        ( netSpecies[(i)].stoich    * *specG[netSpecies[(i)].index] ) // mass fraction weighted Gibbs energy

namespace pokitt{

/**
 *  \class  ReactionRates
 *  \author Nathan Yonkee
 *  \date July, 2014
 *
 *  \brief Calculates net production rates for each species based on
 *  a kinetics mechanism [kg/m^3/s].
 *
 *  This class calculates net reaction rates for
 *  a chemical kinetics mechanism. Currently this class supports
 *  elementary reversible reactions with arrhenius kinetics. Additionally,
 *  this class supports third body reactions and Lindemann, Troe3, and Troe4
 *  falloff reactions. The Gibbs energy is evaluated using a NASA7 polynomial.
 *
 *  For the case of an elementary, reversible reaction the rate is calculated
 *  as the difference of the forward and reverse rates.
 *
 * \f[
 *  r_{net} = r_{forward}-r_{reverse}
 * \f]
 *
 * Units are kg/m^3/s.
 *
 * The forward and reverse rates of progress follow arrhenius kinetics. For example,
 * the forward rate of progress is calculated as,
 *
 * \f[
 * r_f = k_f \prod_{i=0}^{reactants} C_i
 * \f]
 *
 * The reaction rate constant k is calculated as,
 *
 * \f[
 * k = A \exp(-E/RT)
 * \f]
 *
 * And the reversible rate constant is calculated from the change of reaction
 * in Gibbs free energy.
 *
 * The net rate of production of each species is the difference of the
 * rates of progress of reactions in which it is a product or reactant.
 *
 * This class also supports third body reactions. Third body reactions
 * have rates which are enhanced by the presence of non-participating species,
 *
 * \f[
 * r_{enhance} = M r_{base}
 * \f]
 *
 * Where \f$ r_{base} \f$ is calculated with arrhenius kinetics. M is the sum
 * of the concentrations of third body participants with a weighting factor,
 *
 * \f[
 * M = \sum_{i=0}^{thd bodies} m_i C_i
 * \f]
 *
 * This expression also supports Lindemann and Troe falloff reactions.
 *
 */

template< typename FieldT >
class ReactionRates
    : public Expr::Expression<FieldT>
{
  typedef typename Expr::Expression<FieldT>::ValVec SpecT;
  typedef std::vector< RxnData::SpeciesRxnData > SpecDataVecT; // from CanteraObjects.h

  DECLARE_FIELDS( FieldT, t_, rho_, mmw_ )
  DECLARE_VECTOR_OF_FIELDS( FieldT, yi_ )

  std::vector< ReactionInfo  > rxnInfoVec_; // one for each reaction
  std::vector< const RxnData*> rxnDataVec_; // kinetics data extracted form Cantera

  const int nSpec_; // number of species in the mechanism
  const int nRxns_; // number of reactions in the mechanism
  const double invGasConstant_; // inverse of universal gas constant ( division is expensive)

  std::vector< ThermData > specThermVec_;

  const double standardStatePressure_;

  const Expr::Tag& tempTag_, rhoTag_;
  const Expr::TagList& yiTags_;

  Expr::TagList primTags_;

  const bool computeJacobian_;
  const boost::shared_ptr<ChemicalSourceJacobian>& analyticalJacobian_;

  ReactionRates( const Expr::Tag& tTag,
                 const Expr::Tag& rhoTag,
                 const Expr::TagList& yiTags,
                 const Expr::Tag& mmwTag,
                 const boost::shared_ptr<ChemicalSourceJacobian>& aj );
public:
  class Builder : public Expr::ExpressionBuilder
  {
    const Expr::Tag tTag_;
    const Expr::Tag rhoTag_;
    const Expr::TagList yiTags_;
    const Expr::Tag mmwTag_;
    const boost::shared_ptr<ChemicalSourceJacobian> analyticalJacobian_;
  public:
    /**
     *  @brief Build a ReactionRates expression
     *  @param resultTags tags for the net rate of production of each species
     *  @param tTag temperature
     *  @param rhoTag density
     *  @param yiTags tag for the mass fraction of each species
     *  @param mmwTag tag for mixture molecular weight
     *  @param nghost number of ghost cells
     */
    Builder( const Expr::TagList& resultTags,
             const Expr::Tag& tTag,
             const Expr::Tag& rhoTag,
             const Expr::TagList& yiTags,
             const Expr::Tag& mmwTag,
             const boost::shared_ptr<ChemicalSourceJacobian>& analyticalJacobian = boost::shared_ptr<ChemicalSourceJacobian>(),
             const SpatialOps::GhostData nghost = DEFAULT_NUMBER_OF_GHOSTS )
    : ExpressionBuilder( resultTags, nghost ),
      tTag_( tTag ),
      rhoTag_( rhoTag ),
      yiTags_( yiTags ),
      mmwTag_( mmwTag ),
      analyticalJacobian_( analyticalJacobian )
    {}

    Expr::ExpressionBase* build() const{
      return new ReactionRates<FieldT>( tTag_, rhoTag_, yiTags_, mmwTag_, analyticalJacobian_ );
    }
  };

  ~ReactionRates(){}
  void evaluate();
  void sensitivity( const Expr::Tag& );

  inline bool override_sensitivity() const{ return computeJacobian_; }
};



// ###################################################################
//
//                          Implementation
//
// ###################################################################


template< typename FieldT >
ReactionRates<FieldT>::
ReactionRates( const Expr::Tag& tTag,
               const Expr::Tag& rhoTag,
               const Expr::TagList& yiTags,
               const Expr::Tag& mmwTag,
               const boost::shared_ptr<ChemicalSourceJacobian>& analyticalJacobian )
  : Expr::Expression<FieldT>(),
    nSpec_( CanteraObjects::number_species() ),
    nRxns_( CanteraObjects::number_rxns() ),
    invGasConstant_( 1 / CanteraObjects::gas_constant() ),
    standardStatePressure_( CanteraObjects::reference_pressure() ),
    tempTag_( tTag ),
    rhoTag_( rhoTag ),
    yiTags_( yiTags ),
    analyticalJacobian_( analyticalJacobian ),
    computeJacobian_( analyticalJacobian )
{
  this->set_gpu_runnable( true );

  t_   = this->template create_field_request<FieldT>( tTag   );
  rho_ = this->template create_field_request<FieldT>( rhoTag );
  mmw_ = this->template create_field_request<FieldT>( mmwTag );
  this->template create_field_vector_request<FieldT>( yiTags, yi_ );

  for( int r=0; r<nRxns_; ++r){
    const RxnData& rxnDat = CanteraObjects::rxn_data( r );
    rxnDataVec_.push_back( &rxnDat );
    try{
      rxnInfoVec_.push_back( ReactionInfo( rxnDat ) );
    }
    catch( std::runtime_error& err ){
      std::ostringstream msg;
      msg << err.what() << "\n Error detecting while building ReactionInfo for rxn r = " << r << "\n";
      throw( std::runtime_error( msg.str() ) );
    }
  }

  const double gasConstant = CanteraObjects::gas_constant();
  for( int n=0; n<nSpec_; ++n ){
    ThermData specData = CanteraObjects::species_thermo( n );
    std::vector<double>& c = specData.coefficients;
    switch ( specData.type ) {
    case CONST_POLY: break;
    case NASA_POLY:
      for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic )
        *ic *= gasConstant; // dimensionalize the coefficients
      c[ 2] /= 2; // perform division of coefficients here - minor optimization
      c[ 3] /= 6;
      c[ 4] /= 12;
      c[ 5] /= 20;
      c[ 9] /= 2;
      c[10] /= 6;
      c[11] /= 12;
      c[12] /= 20;
      break;
    default:{
      std::ostringstream msg;
      msg << "Error in " __FILE__ << " : " << __LINE__ << std::endl
          << "The Gibbs polynomial type is not supported, see CanteraObjects.h" << std::endl
          << "species n = " << n << "  polynomial type = " << specData.type << std::endl;
      throw std::runtime_error( msg.str() );
      }
    }
    specThermVec_.push_back( specData );
  }

  primTags_.push_back( rhoTag_ );
  primTags_.push_back( tempTag_ );
  for( int i=0; i<nSpec_-1; ++i ){
    primTags_.push_back( yiTags_[i] );
  }
}

//--------------------------------------------------------------------

template< typename FieldT >
void
ReactionRates<FieldT>::
evaluate()
{
  if( computeJacobian_ ){
    return;
  }
  using namespace SpatialOps;
  SpecT& rRates = this->get_value_vec();
  const FieldT& t   = t_  ->field_ref();
  const FieldT& rho = rho_->field_ref();
  const FieldT& mmw = mmw_->field_ref();

  // for each species: if false add the rate; if true accumulate the rate
  std::vector< bool > isInitialized( rRates.size(), false);

  SpatFldPtr<FieldT> logTPtr = SpatialFieldStore::get<FieldT>( t ); // log(t), used in Gibbs evaluation and arrhenius evaluations
  FieldT& logT = *logTPtr;
  logT <<= log( t );
  /* calculate the delta gibbs energy for each species for use in evaluating reversible rate constants
   * polynomials are applicable in two temperature ranges - high and low
   * Caution: polynomials are used even if temperature is out of range
   * Note: coefficients have been divided by integers during construction of expression
   */
  std::vector< SpatFldPtr<FieldT> > specG; // gibbs energy for each species
  for( int n=0; n<nSpec_; ++n ){
    specG.push_back( SpatialFieldStore::get<FieldT>(t) );
    const ThermoPoly& polyType = specThermVec_[n].type;
    const std::vector<double>& c = specThermVec_[n].coefficients;
    switch ( polyType ){
      case NASA_POLY:
        *specG[n] <<= cond( t <= c[0], c[13] + t * ( c[8] - c[14] - c[8] * logT - t * ( c[9] + t * ( c[10] + t * ( c[11] + t * c[12] )))) )  // if low temp
        (            c[ 6] + t * ( c[1] - c[ 7] - c[1] * logT - t * ( c[2] + t * ( c[ 3] + t * ( c[ 4] + t * c[ 5] )))) ); // else if high temp
        break;
      case CONST_POLY:
        *specG[n] <<=       c[1] + c[3] * ( t    - c[0]        ) // H
        - t * ( c[2] + c[3] * ( logT - log(c[0]) ) ); // -TS
        break;
      default: {
        std::ostringstream msg;
        msg << __FILE__ << " : " << __LINE__
            << "\n Unsupported thermo model for species " << n << " somehow evaded detection\n"
            << "This should have been caught during construction of this expression or CanteraObjects\n"
            << "Check your xml input file to ensure it is not corrupted\n";
        throw std::runtime_error( msg.str() );
      }
    } // switch
  } // species loop

  // First, some temporary fields for values ues in every run of the loop
  SpatFldPtr<FieldT> concPtr     = SpatialFieldStore::get<FieldT>( t ); // molar concentration
  SpatFldPtr<FieldT> logConcPtr  = SpatialFieldStore::get<FieldT>( t ); // log of molar concentration
  SpatFldPtr<FieldT> tRecipPtr   = SpatialFieldStore::get<FieldT>( t ); // 1/t
  FieldT& conc    = *concPtr;
  FieldT& logConc = *logConcPtr;
  FieldT& tRecip  = *tRecipPtr;
  conc    <<= rho / mmw;
  logConc <<= log( conc );
  tRecip  <<= 1/t;

  // temporary fields for forward rate constant and falloff reactions
  SpatFldPtr<FieldT> kPtr        = SpatialFieldStore::get<FieldT>( t ); // forward rate constant
  SpatFldPtr<FieldT> mPtr        = SpatialFieldStore::get<FieldT>( t ); // third body enhancement factor
  SpatFldPtr<FieldT> prPtr       = SpatialFieldStore::get<FieldT>( t ); // normalized pressure factor (falloff reactions)
  SpatFldPtr<FieldT> logFCentPtr = SpatialFieldStore::get<FieldT>( t ); // log10(Fcent) (Troe falloff)
  SpatFldPtr<FieldT> logPrCPtr   = SpatialFieldStore::get<FieldT>( t ); // log10(Pr) + CTROE (Troe falloff)

  FieldT& k        = *kPtr;
  FieldT& m        = *mPtr;
  FieldT& pr       = *prPtr;
  FieldT& logPrC   = *logPrCPtr;
  FieldT& logFCent = *logFCentPtr;

  SpatFldPtr<FieldT> krPtr   = SpatialFieldStore::get<FieldT>( t ); // reverse rate constant
  FieldT& kr = *krPtr;

  for( int r=0; r<nRxns_; ++r ){
    const ReactionInfo& rxnInfo = rxnInfoVec_[r];
    const RxnData& rxnDat = *rxnDataVec_[r];
    const SpecDataVecT& thdBodies  = rxnDat.thdBdySpecies;
    const SpecDataVecT& reactants  = rxnDat.reactants;
    const SpecDataVecT& products   = rxnDat.products;
    const SpecDataVecT& netSpecies = rxnDat.netSpecies;


    const double* kCoefs = rxnDat.kFwdCoefs;
    switch( rxnInfo.kForm ){
      case CONSTANT:   k <<= kCoefs[0];          break;
      case LINEAR:     k <<= kCoefs[0] * t;      break;
      case QUADRATIC:  k <<= kCoefs[0] * t * t;  break;
      case RECIPROCAL: k <<= kCoefs[0] * tRecip; break;
      case ARRHENIUS:  k <<= ARRHENIUS( kCoefs ); break;
    }

    const double baseEff  = rxnDat.thdBdyDefault;
    const double* kPCoefs = rxnDat.kPressureCoefs;
    const double* troe    = rxnDat.troeParams;
    switch( rxnDat.type ){
      case ELEMENTARY: break;
      case THIRD_BODY:
        switch( thdBodies.size() ){
          case 0: k <<= k *   baseEff * conc;                                       break;
          case 1: k <<= k * ( baseEff * conc + rho * (            THD_BDY(0) ) );   break;
          case 2: k <<= k * ( baseEff * conc + rho * ( POKITT_SUM2( THD_BDY  ) ) ); break;
          case 3: k <<= k * ( baseEff * conc + rho * ( POKITT_SUM3( THD_BDY  ) ) ); break;
          case 4: k <<= k * ( baseEff * conc + rho * ( POKITT_SUM4( THD_BDY  ) ) ); break;
          case 5: k <<= k * ( baseEff * conc + rho * ( POKITT_SUM5( THD_BDY  ) ) ); break;
          case 6: k <<= k * ( baseEff * conc + rho * ( POKITT_SUM6( THD_BDY  ) ) ); break;
          case 7: k <<= k * ( baseEff * conc + rho * ( POKITT_SUM7( THD_BDY  ) ) ); break;
          case 8: k <<= k * ( baseEff * conc + rho * ( POKITT_SUM8( THD_BDY  ) ) ); break;
          default:
            m <<= baseEff * conc + rho * ( POKITT_SUM8( THD_BDY  ) );
            for( int i = 8; i != thdBodies.size(); ++i ) // if we have 9 or more third bodies
                m <<= m + rho * THD_BDY(i);
            k <<= k * m;
            break;
        }
        break;
          case LINDEMANN:
            switch( thdBodies.size() ){
              case 0: k <<= k / ( 1 + k / ( ARRHENIUS( kPCoefs ) * ( baseEff * conc                                 ) ) );     break;
              case 1: k <<= k / ( 1 + k / ( ARRHENIUS( kPCoefs ) * ( baseEff * conc + rho * (            THD_BDY(0) ) ) ) );   break;
              case 2: k <<= k / ( 1 + k / ( ARRHENIUS( kPCoefs ) * ( baseEff * conc + rho * ( POKITT_SUM2( THD_BDY  ) ) ) ) ); break;
              case 3: k <<= k / ( 1 + k / ( ARRHENIUS( kPCoefs ) * ( baseEff * conc + rho * ( POKITT_SUM3( THD_BDY  ) ) ) ) ); break;
              case 4: k <<= k / ( 1 + k / ( ARRHENIUS( kPCoefs ) * ( baseEff * conc + rho * ( POKITT_SUM4( THD_BDY  ) ) ) ) ); break;
              case 5: k <<= k / ( 1 + k / ( ARRHENIUS( kPCoefs ) * ( baseEff * conc + rho * ( POKITT_SUM5( THD_BDY  ) ) ) ) ); break;
              case 6: k <<= k / ( 1 + k / ( ARRHENIUS( kPCoefs ) * ( baseEff * conc + rho * ( POKITT_SUM6( THD_BDY  ) ) ) ) ); break;
              case 7: k <<= k / ( 1 + k / ( ARRHENIUS( kPCoefs ) * ( baseEff * conc + rho * ( POKITT_SUM7( THD_BDY  ) ) ) ) ); break;
              case 8: k <<= k / ( 1 + k / ( ARRHENIUS( kPCoefs ) * ( baseEff * conc + rho * ( POKITT_SUM8( THD_BDY  ) ) ) ) ); break;
              default:
                m <<= baseEff * conc + rho * ( POKITT_SUM8( THD_BDY ) );
                for( int i = 8; i != thdBodies.size(); ++i )
                  m <<= m + rho * THD_BDY(i);
                k <<= k / ( 1 + k / ( ARRHENIUS( kPCoefs ) * m ) );
                break;
            } // thd bodies
            break;
              case TROE:
                switch( thdBodies.size() ){
                  case 0: pr <<= ARRHENIUS( kPCoefs ) / k * ( baseEff * conc                                    ); break;
                  case 1: pr <<= ARRHENIUS( kPCoefs ) / k * ( baseEff * conc + rho * (            THD_BDY(0)  ) ); break;
                  case 2: pr <<= ARRHENIUS( kPCoefs ) / k * ( baseEff * conc + rho * ( POKITT_SUM2( THD_BDY ) ) ); break;
                  case 3: pr <<= ARRHENIUS( kPCoefs ) / k * ( baseEff * conc + rho * ( POKITT_SUM3( THD_BDY ) ) ); break;
                  case 4: pr <<= ARRHENIUS( kPCoefs ) / k * ( baseEff * conc + rho * ( POKITT_SUM4( THD_BDY ) ) ); break;
                  case 5: pr <<= ARRHENIUS( kPCoefs ) / k * ( baseEff * conc + rho * ( POKITT_SUM5( THD_BDY ) ) ); break;
                  case 6: pr <<= ARRHENIUS( kPCoefs ) / k * ( baseEff * conc + rho * ( POKITT_SUM6( THD_BDY ) ) ); break;
                  case 7: pr <<= ARRHENIUS( kPCoefs ) / k * ( baseEff * conc + rho * ( POKITT_SUM7( THD_BDY ) ) ); break;
                  case 8: pr <<= ARRHENIUS( kPCoefs ) / k * ( baseEff * conc + rho * ( POKITT_SUM8( THD_BDY ) ) ); break;
                  default:
                    m <<= baseEff * conc + rho * ( POKITT_SUM8( THD_BDY ) );
                    for( int i = 8; i != thdBodies.size(); ++i )
                      m <<= m + rho * THD_BDY(i);
                    pr <<= ARRHENIUS( kPCoefs ) / k * m;
                    break;
                } //thd bodies

                switch( rxnInfo.troeForm ){
                  case T123: logFCent <<= log10( (1-troe[0]) * exp(-t/troe[1]) + troe[0] * exp(-t/troe[2]) + exp(-tRecip * troe[3]) ); break;
                  case T12:  logFCent <<= log10( (1-troe[0]) * exp(-t/troe[1]) + troe[0] * exp(-t/troe[2]) + 0.0                    ); break;
                  case T1:   logFCent <<= log10( (1-troe[0]) * exp(-t/troe[1]) +        0.0                + 0.0                    ); break;
                  case T23:  logFCent <<= log10(            0.0                + troe[0] * exp(-t/troe[2]) + exp(-tRecip * troe[3]) ); break;
                  case T2:   logFCent <<= log10(            0.0                + troe[0] * exp(-t/troe[2]) + 0.0                    ); break;
                  case T13:  logFCent <<= log10( (1-troe[0]) * exp(-t/troe[1]) +        0.0                + exp(-tRecip * troe[3]) ); break;
                  case T3:   logFCent <<= log10(            0.0                +        0.0                + exp(-tRecip * troe[3]) ); break;
                  case NONE:
                  default:{
                    std::ostringstream msg;
                    msg << "Error in " << __FILE__ << " : " << __LINE__ << std::endl
                        <<" Falloff type is TROE, but no terms are flagged for evaluation " << std::endl
                        <<" Reaction # "<< r << std::endl;
                    throw std::runtime_error( msg.str() );
                  }
                }

#        define CTROE ( -0.4 - 0.67 * logFCent ) // macros to make f1 evaluation easier to read
#        define NTROE ( 0.75 - 1.27 * logFCent ) // macros to make k evaluation easier to read
#        define F1 ( logPrC / ( NTROE - 0.14 * logPrC ) )

                logPrC <<= log10( pr ) + CTROE;
                k <<= k * pow(10, logFCent / (1 + square( F1 ) ) ) * pr / (1+pr);
                break;
                  default: {
                    std::ostringstream msg;
                    msg << __FILE__ << " : " << __LINE__
                        << "\n Unidentified reaction type for reaction r = " << r << " somehow evaded detection\n"
                        << "This should have been caught during construction of CanteraObjects\n"
                        << "Check your xml input file to ensure it is not corrupted\n";
                    throw std::runtime_error( msg.str() );
                  }
    } //  switch( rxnDat.type )

    if( rxnDat.reversible ){
      const int netOrder = rxnDat.netOrder;
      switch( netSpecies.size() ){
        case 3: kr <<= k * exp( netOrder * log( standardStatePressure_ * tRecip * invGasConstant_ ) - tRecip * invGasConstant_ * ( POKITT_SUM3( GIBBS ) ) ); break;
        case 2: kr <<= k * exp( netOrder * log( standardStatePressure_ * tRecip * invGasConstant_ ) - tRecip * invGasConstant_ * ( POKITT_SUM2( GIBBS ) ) ); break;
        case 4: kr <<= k * exp( netOrder * log( standardStatePressure_ * tRecip * invGasConstant_ ) - tRecip * invGasConstant_ * ( POKITT_SUM4( GIBBS ) ) ); break;
        case 5: kr <<= k * exp( netOrder * log( standardStatePressure_ * tRecip * invGasConstant_ ) - tRecip * invGasConstant_ * ( POKITT_SUM5( GIBBS ) ) ); break;
      }
    }

#   define C_R(i) ( yi_[reactants[i].index]->field_ref() * rho * reactants[i].invMW ) // concentration of the reactants
#   define C_P(i) ( yi_[products [i].index]->field_ref() * rho * products [i].invMW ) // concentration of products
    switch( rxnInfo.forwardOrder ){
      case ONE:         k <<= k  *         C_R(0);                              break;
      case TWO:         k <<= k  * square( C_R(0) );                            break;
      case ONE_ONE:     k <<= k  *         C_R(0)   *         C_R(1);           break;
      case ONE_ONE_ONE: k <<= k  *         C_R(0)   *         C_R(1)  * C_R(2); break;
      case TWO_ONE:     k <<= k  * square( C_R(0) ) *         C_R(1);           break;
      case ONE_TWO:     k <<= k  *         C_R(0)   * square( C_R(1) );         break;
      default:
        for( int i = 0; i != reactants.size(); ++i){
          switch( reactants[i].stoich ){
            case 1: k <<= k *         C_R(i)  ; break;
            case 2: k <<= k * square( C_R(i) ); break;
            case 3: k <<= k * cube(   C_R(i) ); break;
          }
        }
        break;
    }

    if( rxnDat.reversible ){
      switch( rxnInfo.reverseOrder ){
        case ONE:         kr <<= kr *         C_P(0);                            break;
        case TWO:         kr <<= kr * square( C_P(0) );                          break;
        case ONE_ONE:     kr <<= kr *         C_P(0)         * C_P(1);           break;
        case ONE_ONE_ONE: kr <<= kr *         C_P(0)         * C_P(1)  * C_P(2); break;
        case TWO_ONE:     kr <<= kr * square( C_P(0) )       * C_P(1);           break;
        case ONE_TWO:     kr <<= kr *         C_P(0) * square( C_P(1) );         break;
        default:
          for( int i = 0; i != products.size(); ++i){
            switch( std::abs(products[i].stoich) ){
              case 1: kr <<= kr *         C_P(i)  ; break;
              case 2: kr <<= kr * square( C_P(i) ); break;
              case 3: kr <<= kr * cube(   C_P(i) ); break;
            }
          }
          break;
      }
    }

    SpecDataVecT::const_iterator iNet = netSpecies.begin();
    for( ; iNet != netSpecies.end(); ++iNet ){
      const int index = iNet->index;
      FieldT& ri = *rRates[index];
      if( !isInitialized[index] ){
        isInitialized[index] = true;
        if( rxnDat.reversible ) ri <<=    - ( iNet->stoich * iNet->mw ) * ( k - kr );
        else                    ri <<=    - ( iNet->stoich * iNet->mw ) * ( k      );
      }
      else{
        if( rxnDat.reversible ) ri <<= ri - ( iNet->stoich * iNet->mw ) * ( k - kr );
        else                    ri <<= ri - ( iNet->stoich * iNet->mw ) * ( k      );
      }
    } // sum reaction rates into products and reactants

  } // loop over reactions

  for( int n = 0; n < nSpec_; ++n ){
    if( !isInitialized[n] ) *rRates[n] <<= 0.0; // in case of inerts that have 0.0 rate
  }
}

//--------------------------------------------------------------------

template< typename FieldT >
void
ReactionRates<FieldT>::
sensitivity( const Expr::Tag& tag )
{
  if( computeJacobian_ ){
    const Expr::TagList& rateTags = this->exprNames_;
    if( tag == rhoTag_ ){
      SpecT& rRates = this->get_value_vec();
      const FieldT& T = t_->field_ref();
      const FieldT& rho = rho_->field_ref();
      const FieldT& Mmix = mmw_->field_ref();
      std::vector<const FieldT*> Y;
      for( int i=0; i<nSpec_; ++i ){
        Y.push_back( &( yi_[i]->field_ref() ) );
      }
      std::vector<SpatialOps::SpatFldPtr<FieldT> > drdvPtrs;
      for( int row=0; row<nSpec_; ++row ){
        for( int col=0; col<nSpec_+1; ++col ){
          drdvPtrs.push_back( this->sensitivity_result_ptr( rateTags[row], primTags_[col] ) );
        }
      }
      for( int col=0; col<nSpec_+1; ++col ){
        drdvPtrs.push_back( SpatialOps::SpatialFieldStore::get<FieldT>( *rRates[0] ) );
      }

      SpatialOps::FieldVector<FieldT> rates( rRates );
      SpatialOps::FieldMatrix<FieldT> drdv( drdvPtrs );

      analyticalJacobian_->evaluate_rates_and_jacobian( drdv, rates, T, rho, Y, Mmix );
    }
    else if(std::find(primTags_.begin(), primTags_.end(), tag) == primTags_.end()){
      for( int row=0; row<nSpec_; ++row ){
        this->sensitivity_result( rateTags[row], tag ) <<= 0.0;
      }
    }
  }
}

//--------------------------------------------------------------------

} // namespace pokitt

#undef POKITT_SUM2
#undef POKITT_SUM3
#undef POKITT_SUM4
#undef POKITT_SUM5
#undef POKITT_SUM6
#undef POKITT_SUM7
#undef POKITT_SUM8
#undef POKITT_SUM9

#undef ARRHENIUS
#undef THD_BDY
#undef CTROE
#undef NTROE
#undef F1
#undef GIBBS
#undef C_R
#undef C_P

#endif // ReactionRates_Expr_h
