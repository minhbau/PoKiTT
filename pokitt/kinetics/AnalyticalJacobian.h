/*
 * The MIT License
 *
 * Copyright (c) 2016-2017 The University of Utah
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

#ifndef AnalyticalJacobian_h
#define AnalyticalJacobian_h

/**
 *  \file   AnalyticalJacobian.h
 *  \date   April 11, 2016
 *  \author mike
 */

#include <numeric>
#include <iostream>
#include <vector>

#include <cmath>

#include <pokitt/CanteraObjects.h>
#include <pokitt/kinetics/ReactionInfo.h>

#include <spatialops/Nebo.h>
#include <spatialops/structured/Grid.h>
#include <spatialops/structured/MatVecFields.h>
#include <spatialops/structured/MatVecOps.h>

#include <expression/ExprLib.h>
#include <expression/Expression.h>

namespace pokitt{

#ifdef ENABLE_CUDA
# define LOCATION GPU_INDEX
#else
# define LOCATION CPU_INDEX
#endif

typedef SpatialOps::SVolField FieldT;

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

#define ARRHENIUS( coef ) ( (coef)[0] * exp( (coef)[1] * lnT - (coef)[2] * invT ) )
#define ARRHENIUS_SENS_TEMP( coef ) ( invT * ( (coef)[1] + (coef)[2] * invT ) ) // \frac{1}{k}\pder{k}{T} for a rate constant k and temperature T

#define THD_BDY( i ) ( thdBodies[(i)].thdBdyEff * *YPtr[thdBodies[(i)].index] )
#define THD_BDY_NOY( i ) ( thdBodies[(i)].thdBdyEff )

#define GIBBS( i ) ( netSpecies[(i)].stoich * *specGPtr[netSpecies[(i)].index] )
#define DBDT( i )  ( netSpecies[(i)].stoich * *dBdTSpecPtr[netSpecies[(i)].index] )



/**
 *  \class  ChemicalSourceJacobian
 *  \author Mike Hansen
 *  \date February, 2016
 *
 *  \brief Calculates net production rates for each species as well as the
 *  sensitivities of each species production rates to the primitive variables,
 *  density, temperature, and n-1 mass fractions.
 *
 *  This class supports elementary reversible reactions with modified Arrhenius,
 *  three-body, and Falloff reactions of Lindemann and Troe forms.
 *
 *  NASA7 polynomials and constant heat capacities are supported.
 *
 *  Units of the production rates are kg/m^3/s.
 *
 *  See the writeup for details.
 */

class ChemicalSourceJacobian
{
private:
  void init()
  {
    for( size_t s=0; s<ns; ++s ){
      invMsp.push_back( 1.0 / Msp[s] );
    }

    for( size_t r=0; r<nr; ++r){

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

    for( int n=0; n<ns; ++n ){
      ThermData specData = CanteraObjects::species_thermo( n );
      std::vector<double>& c = specData.coefficients;
      switch ( specData.type ) {
        case CONST_POLY: break;
        case NASA_POLY:
          for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic )
            *ic *= Ru; // dimensionalize the coefficients
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
  }

public:
  /**
     *  @brief Build a ChemicalSourceJacobian object
     */
    ChemicalSourceJacobian()
    : ns( CanteraObjects::number_species() ),
      nr( CanteraObjects::number_rxns() ),
      Msp( CanteraObjects::molecular_weights() ),
      Ru( CanteraObjects::gas_constant() ),
      invRu( 1.0 / CanteraObjects::gas_constant() ),
      pref( CanteraObjects::reference_pressure() )
    {
      init();
    } // constructor



  const size_t ns;
  const size_t nr;

  const std::vector< double > Msp;
  std::vector< double > invMsp;

  const double Ru;
  const double invRu;

  const double pref;

  std::vector< ReactionInfo >  rxnInfoVec_;   // one for each reaction
  std::vector< const RxnData*> rxnDataVec_;   // kinetics data extracted from Cantera
  std::vector< ThermData >     specThermVec_; // species thermo data

  /**
   *  @brief Evaluate species production rates and their sensitivities
   *  @param primitiveSensitivities a FieldMatrix for the sensitivities of prod. rates to primitives
   *  @param productionRates a FieldVector for the production rates
   *  @param T a field containing the temperature
   *  @param rho a field containing the density
   *  @param YPtr a vector of SpatFldPtr containing the mass fraction fields
   *  @param Mmix a field containing the mixture molecular weight
   *
   * Primitive sensitivities has the sensitivities of the ns production rates on its first ns rows.
   * The (ns+1)-th row is left empty.
   *
   * primitiveSensitivities(i,j) will give \pder{w_i}{V_j}.
   *
   * primitiveSensitivities =
   * \pder{w_1}{\rho}     & \pder{w_1}{T}     & \pder{w_1}{Y_1}     & \pder{w_1}{Y_2}     & \cdots & \pder{w_1}{Y_{n_s-1}}     \\
   * \pder{w_2}{\rho}     & \pder{w_2}{T}     & \pder{w_2}{Y_1}     & \pder{w_2}{Y_2}     & \cdots & \pder{w_2}{Y_{n_s-1}}     \\
   * \vdots               & \vdots            & \vdots              & \vdots              & \ddots & \vdots                    \\
   * \pder{w_{n_s}}{\rho} & \pder{w_{n_s}}{T} & \pder{w_{n_s}}{Y_1} & \pder{w_{n_s}}{Y_2} & \cdots & \pder{w_{n_s}}{Y_{n_s-1}} \\
   * 0                    & 0                 & 0                   & 0                   & \cdots & 0
   *
   */
  template< typename FieldT >
  void evaluate_rates_and_jacobian( SpatialOps::FieldMatrix<FieldT>& primitiveSensitivities,
                                    SpatialOps::FieldVector<FieldT>& productionRates,
                                    const FieldT& T,
                                    const FieldT& rho,
                                    const std::vector< const FieldT* >& YPtr,
                                    const FieldT& Mmix ) const;

};


template< typename FieldT >
void
ChemicalSourceJacobian::
evaluate_rates_and_jacobian( SpatialOps::FieldMatrix<FieldT>& primitiveSensitivities,
                             SpatialOps::FieldVector<FieldT>& productionRates,
                             const FieldT& T,
                             const FieldT& rho,
                             const std::vector< const FieldT* >& YPtr,
                             const FieldT& Mmix ) const
{
  // scratch fields - scalars
  SpatialOps::SpatFldPtr<FieldT> invTPtr        = SpatialOps::SpatialFieldStore::get<FieldT>( T ); FieldT& invT        = *invTPtr;        // inverse of temperature
  SpatialOps::SpatFldPtr<FieldT> lnTPtr         = SpatialOps::SpatialFieldStore::get<FieldT>( T ); FieldT& lnT         = *lnTPtr;         // natural logarithm of temperature
  SpatialOps::SpatFldPtr<FieldT> ctPtr          = SpatialOps::SpatialFieldStore::get<FieldT>( T ); FieldT& ct          = *ctPtr;          // molar density
  SpatialOps::SpatFldPtr<FieldT> invMPtr        = SpatialOps::SpatialFieldStore::get<FieldT>( T ); FieldT& invM        = *invMPtr;        // inverse of mixture molecular weight
  SpatialOps::SpatFldPtr<FieldT> kfPtr          = SpatialOps::SpatialFieldStore::get<FieldT>( T ); FieldT& kf          = *kfPtr;          // forward rate constant
  SpatialOps::SpatFldPtr<FieldT> KcPtr          = SpatialOps::SpatialFieldStore::get<FieldT>( T ); FieldT& Kc          = *KcPtr;          // equilibrium constant
  SpatialOps::SpatFldPtr<FieldT> krPtr          = SpatialOps::SpatialFieldStore::get<FieldT>( T ); FieldT& kr          = *krPtr;          // reverse rate constant
  SpatialOps::SpatFldPtr<FieldT> RrPtr          = SpatialOps::SpatialFieldStore::get<FieldT>( T ); FieldT& Rr          = *RrPtr;          // reverse mass action rate
  SpatialOps::SpatFldPtr<FieldT> RnetPtr        = SpatialOps::SpatialFieldStore::get<FieldT>( T ); FieldT& Rnet        = *RnetPtr;        // net mass action rate
  SpatialOps::SpatFldPtr<FieldT> qPtr           = SpatialOps::SpatialFieldStore::get<FieldT>( T ); FieldT& q           = *qPtr;           // rate of progress, q = Rnet * C
  SpatialOps::SpatFldPtr<FieldT> CtbafPtr       = SpatialOps::SpatialFieldStore::get<FieldT>( T ); FieldT& Ctbaf       = *CtbafPtr;       // third body and falloff modifier
  SpatialOps::SpatFldPtr<FieldT> prPtr          = SpatialOps::SpatialFieldStore::get<FieldT>( T ); FieldT& pr          = *prPtr;          // falloff reduced pressure
  SpatialOps::SpatFldPtr<FieldT> fCentPtr       = SpatialOps::SpatialFieldStore::get<FieldT>( T ); FieldT& fCent       = *fCentPtr;       // Fcent (Troe falloff)
  SpatialOps::SpatFldPtr<FieldT> flfConcPtr     = SpatialOps::SpatialFieldStore::get<FieldT>( T ); FieldT& flfConc     = *flfConcPtr;     // falloff third-body concentration
  SpatialOps::SpatFldPtr<FieldT> fTroePtr       = SpatialOps::SpatialFieldStore::get<FieldT>( T ); FieldT& fTroe       = *fTroePtr;       // falloff blending factor for Troe reactions
  SpatialOps::SpatFldPtr<FieldT> gTroePtr       = SpatialOps::SpatialFieldStore::get<FieldT>( T ); FieldT& gTroe       = *gTroePtr;       // piece of the blending factor for Troe reactions
  SpatialOps::SpatFldPtr<FieldT> dRnetdrhoPtr   = SpatialOps::SpatialFieldStore::get<FieldT>( T ); FieldT& dRnetdrho   = *dRnetdrhoPtr;   // \pder{R_net}{\rho}
  SpatialOps::SpatFldPtr<FieldT> dRnetdTPtr     = SpatialOps::SpatialFieldStore::get<FieldT>( T ); FieldT& dRnetdT     = *dRnetdTPtr;     // \pder{R_net}{T}
  SpatialOps::SpatFldPtr<FieldT> dKcdToverKcPtr = SpatialOps::SpatialFieldStore::get<FieldT>( T ); FieldT& dKcdToverKc = *dKcdToverKcPtr; // \pder{K_c}{T}
  SpatialOps::SpatFldPtr<FieldT> dCtbafdrhoPtr  = SpatialOps::SpatialFieldStore::get<FieldT>( T ); FieldT& dCtbafdrho  = *dCtbafdrhoPtr;  // \pder{Ctbaf}{\rho}
  SpatialOps::SpatFldPtr<FieldT> dCtbafdTPtr    = SpatialOps::SpatialFieldStore::get<FieldT>( T ); FieldT& dCtbafdT    = *dCtbafdTPtr;    // \pder{Ctbaf}{T}
  SpatialOps::SpatFldPtr<FieldT> dqdrhoPtr      = SpatialOps::SpatialFieldStore::get<FieldT>( T ); FieldT& dqdrho      = *dqdrhoPtr;      // \pder{q}{\rho}
  SpatialOps::SpatFldPtr<FieldT> dqdTPtr        = SpatialOps::SpatialFieldStore::get<FieldT>( T ); FieldT& dqdT        = *dqdTPtr;        // \pder{q}{T}
  SpatialOps::SpatFldPtr<FieldT> dfTroedrhoPtr  = SpatialOps::SpatialFieldStore::get<FieldT>( T ); FieldT& dfTroedrho  = *dfTroedrhoPtr;  // \pder{fTroe}{\rho}
  SpatialOps::SpatFldPtr<FieldT> dfTroedTPtr    = SpatialOps::SpatialFieldStore::get<FieldT>( T ); FieldT& dfTroedT    = *dfTroedTPtr;    // \pder{fTroe}{T}
  SpatialOps::SpatFldPtr<FieldT> dfCentdTPtr    = SpatialOps::SpatialFieldStore::get<FieldT>( T ); FieldT& dfCentdT    = *dfCentdTPtr;    // \pder{fCent}{T}
  SpatialOps::SpatFldPtr<FieldT> aTroePtr       = SpatialOps::SpatialFieldStore::get<FieldT>( T ); FieldT& aTroe       = *aTroePtr;       // Troe A, part of the blending factor
  SpatialOps::SpatFldPtr<FieldT> bTroePtr       = SpatialOps::SpatialFieldStore::get<FieldT>( T ); FieldT& bTroe       = *bTroePtr;       // Troe B, part of the blending factor
  SpatialOps::SpatFldPtr<FieldT> nsTmpPtr       = SpatialOps::SpatialFieldStore::get<FieldT>( T ); FieldT& nsTmp       = *nsTmpPtr;       // temporary for nth-species calculations

  // scratch fields - vectors
  std::vector< SpatialOps::SpatFldPtr<FieldT> > specGPtr;    // gibbs energy for each species, g_i = h_i - Ts_i, a vector of ns quantities
  std::vector< SpatialOps::SpatFldPtr<FieldT> > dBdTSpecPtr; // \pder{B_i}{T}                , used in dKc/dT  , a vector of ns quantities
  std::vector< SpatialOps::SpatFldPtr<FieldT> > dCtbafdYPtr; // \pder{Ctbaf}{Y_i}                              , a vector of ns quantities
  std::vector< SpatialOps::SpatFldPtr<FieldT> > dfTroedYPtr; // \pder{fTroe}{Y_i}                              , a vector of ns quantities
  std::vector< SpatialOps::SpatFldPtr<FieldT> > dRnetdYPtr;  // \pder{Rnet}{Y_i}                               , a vector of ns-1 quantities
  std::vector< SpatialOps::SpatFldPtr<FieldT> > dqdYPtr;     // \pder{q}{Y_i}                                  , a vector of ns-1 quantities

  for( size_t i=0; i<ns-1; ++i ){ // vectors with (ns-1) x 1 values
    dRnetdYPtr .push_back( SpatialOps::SpatialFieldStore::get<FieldT>( T ) ); *dRnetdYPtr[i] <<= 0.0;
    dqdYPtr    .push_back( SpatialOps::SpatialFieldStore::get<FieldT>( T ) ); *dqdYPtr[i]    <<= 0.0;
  }
  for( size_t i=0; i<ns; ++i ){   // vectors with ns x 1 values
    dCtbafdYPtr  .push_back( SpatialOps::SpatialFieldStore::get<FieldT>( T ) ); *dCtbafdYPtr[i] <<= 0.0;
    dfTroedYPtr  .push_back( SpatialOps::SpatialFieldStore::get<FieldT>( T ) ); *dfTroedYPtr[i] <<= 0.0;
    specGPtr     .push_back( SpatialOps::SpatialFieldStore::get<FieldT>( T ) );
    dBdTSpecPtr  .push_back( SpatialOps::SpatialFieldStore::get<FieldT>( T ) );
  }

  invT <<= 1 / T;
  lnT  <<= log( T );

  invM <<= 1 / Mmix;
  ct <<= rho * invM;

  for( int n=0; n<ns; ++n ){
    const ThermoPoly& polyType = specThermVec_[n].type;
    const std::vector<double>& c = specThermVec_[n].coefficients;
    switch ( polyType ){
      case NASA_POLY:
        *specGPtr[n]    <<= cond( T <= c[0], c[13] + T * ( c[8] - c[14] - c[8] * lnT - T * ( c[9] + T * ( c[10] + T * ( c[11] + T * c[12] )))) )  // if low temp
                                (            c[ 6] + T * ( c[1] - c[ 7] - c[1] * lnT - T * ( c[2] + T * ( c[ 3] + T * ( c[ 4] + T * c[ 5] )))) ); // else if high temp

        *dBdTSpecPtr[n] <<= cond( T <= c[0], invRu * ( ( c[8] - Ru ) * invT + c[9] + T * ( 2 * c[10] + T * ( 3 * c[11] + T * 4 * c[12] ) ) + c[13] * invT * invT ) )  // if low temp
                                (            invRu * ( ( c[1] - Ru ) * invT + c[2] + T * ( 2 * c[ 3] + T * ( 3 * c[ 4] + T * 4 * c[ 5] ) ) + c[ 6] * invT * invT ) ); // else if high temp

        break;
      case CONST_POLY:
        *specGPtr[n]    <<= c[1] + c[3] * ( T - c[0] ) - T * ( c[2] + c[3] * ( lnT - log(c[0]) ) );
        *dBdTSpecPtr[n] <<= invT * ( Msp[n] * invRu * ( c[3] - invT * ( c[3]*c[0] - c[1] ) ) - 1 );
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

  std::vector< bool > isInitialized( ns, false );
  const double rxnOrderTol = ReactionInfo::rxnOrderTol;

  for( size_t r=0; r<nr; ++r ){
    const ReactionInfo& rxnInfo    = rxnInfoVec_[r];
    const RxnData& rxnDat          = *rxnDataVec_[r];
    const SpecDataVecT& thdBodies  = rxnDat.thdBdySpecies;
    const SpecDataVecT& reactants  = rxnDat.reactants;
    const SpecDataVecT& products   = rxnDat.products;
    const SpecDataVecT& netSpecies = rxnDat.netSpecies;

    Rnet      <<= 0;
    dRnetdrho <<= 0;
    dRnetdT   <<= 0;

    dqdrho <<= 0;
    dqdT   <<= 0;

    dCtbafdrho <<= 0.0;
    dCtbafdT   <<= 0.0;

    for( size_t i=0; i<ns-1; ++i ){
      *dqdYPtr[i]     <<= 0.0;
      *dRnetdYPtr[i]  <<= 0.0;
      *dCtbafdYPtr[i] <<= 0.0;
    }

    const double* kCoefs = rxnDat.kFwdCoefs;
    switch( rxnInfo.kForm ){
      case CONSTANT:   kf <<= kCoefs[0];           break;
      case LINEAR:     kf <<= kCoefs[0] * T;       break;
      case QUADRATIC:  kf <<= kCoefs[0] * T * T;   break;
      case RECIPROCAL: kf <<= kCoefs[0] * invT;    break;
      case ARRHENIUS:  kf <<= ARRHENIUS( kCoefs ); break;
    }

#   define C_R(i) ( *YPtr[reactants[i].index] * rho * reactants[i].invMW ) // reactant concentration
#   define C_P(i) ( *YPtr[products [i].index] * rho * products [i].invMW ) // product concentration
    switch( rxnInfo.forwardOrder ){
      case ONE:         Rnet <<= kf  *         C_R(0);                              break;
      case TWO:         Rnet <<= kf  * square( C_R(0) );                            break;
      case ONE_ONE:     Rnet <<= kf  *         C_R(0)   *         C_R(1);           break;
      case ONE_ONE_ONE: Rnet <<= kf  *         C_R(0)   *         C_R(1)  * C_R(2); break;
      case TWO_ONE:     Rnet <<= kf  * square( C_R(0) ) *         C_R(1);           break;
      case ONE_TWO:     Rnet <<= kf  *         C_R(0)   * square( C_R(1) );         break;
      default:
        Rnet <<= kf;
        for( int i = 0; i != reactants.size(); ++i){
          const double rStoich = reactants[i].stoich;
          if     ( fabs( rStoich - 1 ) < rxnOrderTol ) Rnet <<= Rnet *         C_R(i)           ;
          else if( fabs( rStoich - 2 ) < rxnOrderTol ) Rnet <<= Rnet * square( C_R(i)          );
          else if( fabs( rStoich - 3 ) < rxnOrderTol ) Rnet <<= Rnet * cube  ( C_R(i)          );
          else                                         Rnet <<= Rnet * pow   ( C_R(i), rStoich );
        }
        break;
    }

    dRnetdrho <<= Rnet / ct * invM * rxnInfo.sumReactantStoich;
    dRnetdT   <<= Rnet * ARRHENIUS_SENS_TEMP( kCoefs );

    bool nsIsReactant    = false;
    int nsReactantIdx = -1;

    for( size_t sridx=0; sridx<reactants.size(); ++sridx ){
      size_t s = reactants[sridx].index;
      // sridx is the species index in the local reactant list
      // s     is the species index in the global species list

      if( s != ns-1 ){ // to ignore the last species on the first pass
        switch( rxnInfo.forwardOrder ){
          case ONE:
            *dRnetdYPtr[s] <<= kf * rho * invMsp[s];
            break;
          case TWO:
            *dRnetdYPtr[s] <<= kf * rho * invMsp[s] * 2 * C_R(sridx);
            break;
          case ONE_ONE:
            switch( sridx ){
              case 0: *dRnetdYPtr[s] <<= kf * rho * invMsp[s] * C_R(1); break;
              case 1: *dRnetdYPtr[s] <<= kf * rho * invMsp[s] * C_R(0); break;
            }
            break;
              case ONE_ONE_ONE:
                switch( sridx ){
                  case 0: *dRnetdYPtr[s] <<= kf * rho * invMsp[s] * C_R(1) * C_R(2); break;
                  case 1: *dRnetdYPtr[s] <<= kf * rho * invMsp[s] * C_R(0) * C_R(2); break;
                  case 2: *dRnetdYPtr[s] <<= kf * rho * invMsp[s] * C_R(0) * C_R(1); break;
                }
                break;
                  case TWO_ONE:
                    switch( sridx ){
                      case 0: *dRnetdYPtr[s] <<= kf * rho * invMsp[s] * 2 * C_R(0) * C_R(1); break;
                      case 1: *dRnetdYPtr[s] <<= kf * rho * invMsp[s] * square( C_R(0) );    break;
                    }
                    break;
                      case ONE_TWO:
                        switch( sridx ){
                          case 0: *dRnetdYPtr[s] <<= kf * rho * invMsp[s] * square( C_R(1) );    break;
                          case 1: *dRnetdYPtr[s] <<= kf * rho * invMsp[s] * 2 * C_R(1) * C_R(0); break;
                        }
                        break;
                          default:
                            *dRnetdYPtr[s] <<= kf * rho * invMsp[s];
                            const double pStoich = std::abs( products[sridx].stoich );
                            // we don't need to do anything if pStoich == 1 (within the specified tolerance)
                            if     ( fabs( pStoich - 2 ) < rxnOrderTol ) *dRnetdYPtr[s] <<= *dRnetdYPtr[s] * 2       *         C_R(sridx)             ;
                            else if( fabs( pStoich - 3 ) < rxnOrderTol ) *dRnetdYPtr[s] <<= *dRnetdYPtr[s] * 3       * square( C_R(sridx)            );
                            else                                         *dRnetdYPtr[s] <<= *dRnetdYPtr[s] * pStoich * pow   ( C_R(sridx), pStoich-1 );
                            for( int i = 0; i != reactants.size(); ++i){
                              if( !( i==sridx ) ){
                                const double rStoich = reactants[i].stoich;
                                if     ( fabs( rStoich - 1 ) < rxnOrderTol ) *dRnetdYPtr[s] <<= *dRnetdYPtr[s] *         C_R(i)           ;
                                else if( fabs( rStoich - 2 ) < rxnOrderTol ) *dRnetdYPtr[s] <<= *dRnetdYPtr[s] * square( C_R(i)          );
                                else if( fabs( rStoich - 3 ) < rxnOrderTol ) *dRnetdYPtr[s] <<= *dRnetdYPtr[s] * cube  ( C_R(i)          );
                                else                                         *dRnetdYPtr[s] <<= *dRnetdYPtr[s] * pow   ( C_R(i), rStoich );
                              }
                            }
        }
      }
      else{ // nth species is a reactant, so save off its index to prevent a loop below
        nsIsReactant  = true;
        nsReactantIdx = sridx;
      }
    }
    // compute the ns correction term
    if( nsIsReactant ){
      switch( rxnInfo.forwardOrder ){
        case ONE:
          nsTmp <<= kf * rho * invMsp[ns-1];
          break;
        case TWO:
          nsTmp <<= kf * rho * invMsp[ns-1] * 2 * C_R(nsReactantIdx);
          break;
        case ONE_ONE:
          switch( nsReactantIdx ){
            case 0: nsTmp <<= kf * rho * invMsp[ns-1] * C_R(1); break;
            case 1: nsTmp <<= kf * rho * invMsp[ns-1] * C_R(0); break;
          }
          break;
            case ONE_ONE_ONE:
              switch( nsReactantIdx ){
                case 0: nsTmp <<= kf * rho * invMsp[ns-1] * C_R(1) * C_R(2); break;
                case 1: nsTmp <<= kf * rho * invMsp[ns-1] * C_R(0) * C_R(2); break;
                case 2: nsTmp <<= kf * rho * invMsp[ns-1] * C_R(0) * C_R(1); break;
              }
              break;
                case TWO_ONE:
                  switch( nsReactantIdx ){
                    case 0: nsTmp <<= kf * rho * invMsp[ns-1] * 2 * C_R(0) * C_R(1); break;
                    case 1: nsTmp <<= kf * rho * invMsp[ns-1] * square( C_R(0) );    break;
                  }
                  break;
                    case ONE_TWO:
                      switch( nsReactantIdx ){
                        case 0: nsTmp <<= kf * rho * invMsp[ns-1] * square( C_R(1) );    break;
                        case 1: nsTmp <<= kf * rho * invMsp[ns-1] * 2 * C_R(1) * C_R(0); break;
                      }
                      break;
                        default:
                          nsTmp <<= kf * rho * invMsp[ns-1];
                          const double rStoich = reactants[nsReactantIdx].stoich;
                          // we don't need to do anything if rStoich == 1 (within the specified tolerance)
                          if     ( fabs( rStoich - 2 ) < rxnOrderTol ) nsTmp <<= nsTmp * 2       *         C_R(nsReactantIdx)           ;
                          else if( fabs( rStoich - 3 ) < rxnOrderTol ) nsTmp <<= nsTmp * 3       * square( C_R(nsReactantIdx)          );
                          else                                         nsTmp <<= nsTmp * rStoich * pow   ( C_R(nsReactantIdx), rStoich );
                          for( int i = 0; i != reactants.size(); ++i){
                            if( !( i==nsReactantIdx ) ){
                              if     ( fabs( rStoich - 1 ) < rxnOrderTol ) nsTmp <<= nsTmp *         C_R(i)           ;
                              else if( fabs( rStoich - 2 ) < rxnOrderTol ) nsTmp <<= nsTmp * square( C_R(i)          );
                              else if( fabs( rStoich - 3 ) < rxnOrderTol ) nsTmp <<= nsTmp * cube  ( C_R(i)          );
                              else                                         nsTmp <<= nsTmp * pow   ( C_R(i), rStoich );
                            }
                          }
      }
      // subtract from all species sensitivities
      for( size_t s=0; s<ns-1; ++s ){
        *dRnetdYPtr[s] <<= *dRnetdYPtr[s] - nsTmp;
      }
    }

    if( rxnDat.reversible ){
      const double netOrder = rxnDat.netOrder;
      switch( netSpecies.size() ){
        case 3:
          Kc          <<= exp( - ( netOrder * log( pref * invT * invRu ) - invT * invRu * ( POKITT_SUM3( GIBBS ) ) ) );
          dKcdToverKc <<= - POKITT_SUM3( DBDT );
          break;
        case 2:
          Kc          <<= exp( - ( netOrder * log( pref * invT * invRu ) - invT * invRu * ( POKITT_SUM2( GIBBS ) ) ) );
          dKcdToverKc <<= - POKITT_SUM2( DBDT );
          break;
        case 4:
          Kc          <<= exp( - ( netOrder * log( pref * invT * invRu ) - invT * invRu * ( POKITT_SUM4( GIBBS ) ) ) );
          dKcdToverKc <<= - POKITT_SUM4( DBDT );
          break;
        case 5:
          Kc          <<= exp( - ( netOrder * log( pref * invT * invRu ) - invT * invRu * ( POKITT_SUM5( GIBBS ) ) ) );
          dKcdToverKc <<= - POKITT_SUM5( DBDT );
          break;
      }
      kr <<= kf / Kc;

      switch( rxnInfo.reverseOrder ){
        case ONE:         Rr <<= kr *         C_P(0);                            break;
        case TWO:         Rr <<= kr * square( C_P(0) );                          break;
        case ONE_ONE:     Rr <<= kr *         C_P(0)         * C_P(1);           break;
        case ONE_ONE_ONE: Rr <<= kr *         C_P(0)         * C_P(1)  * C_P(2); break;
        case TWO_ONE:     Rr <<= kr * square( C_P(0) )       * C_P(1);           break;
        case ONE_TWO:     Rr <<= kr *         C_P(0) * square( C_P(1) );         break;
        default:
          for( int i = 0; i != products.size(); ++i){
            const double pStoich = std::abs( products[i].stoich );
            if     ( fabs( pStoich - 1 ) < ReactionInfo::rxnOrderTol ) Rr <<= kr *         C_P(i)           ;
            else if( fabs( pStoich - 2 ) < ReactionInfo::rxnOrderTol ) Rr <<= kr * square( C_P(i)          );
            else if( fabs( pStoich - 3 ) < ReactionInfo::rxnOrderTol ) Rr <<= kr * cube  ( C_P(i)          );
            else                                                       Rr <<= kr * pow   ( C_P(i), pStoich );
          }
          break;
      }

      Rnet      <<= Rnet - Rr;
      dRnetdrho <<= dRnetdrho - Rr / ct * invM * rxnInfo.sumProductStoich;
      dRnetdT   <<= dRnetdT   - Rr * ( ARRHENIUS_SENS_TEMP( kCoefs ) - dKcdToverKc );

      nsIsReactant  = false;
      nsReactantIdx = -1;

      for( size_t sridx=0; sridx<products.size(); ++sridx ){
        size_t s = products[sridx].index;
        // sridx is the species index in the product list
        // s     is the species index in the species list

        if( s != ns-1 ){ // to ignore the last species on the first pass
          switch( rxnInfo.reverseOrder ){
            case ONE:
              *dRnetdYPtr[s] <<= *dRnetdYPtr[s] - kr * rho * invMsp[s];
              break;
            case TWO:
              *dRnetdYPtr[s] <<= *dRnetdYPtr[s] - kr * rho * invMsp[s] * 2 * C_P(sridx);
              break;
            case ONE_ONE:
              switch( sridx ){
                case 0: *dRnetdYPtr[s] <<= *dRnetdYPtr[s] - kr * rho * invMsp[s] * C_P(1); break;
                case 1: *dRnetdYPtr[s] <<= *dRnetdYPtr[s] - kr * rho * invMsp[s] * C_P(0); break;
              }
              break;
                case ONE_ONE_ONE:
                  switch( sridx ){
                    case 0: *dRnetdYPtr[s] <<= *dRnetdYPtr[s] - kr * rho * invMsp[s] * C_P(1) * C_P(2); break;
                    case 1: *dRnetdYPtr[s] <<= *dRnetdYPtr[s] - kr * rho * invMsp[s] * C_P(0) * C_P(2); break;
                    case 2: *dRnetdYPtr[s] <<= *dRnetdYPtr[s] - kr * rho * invMsp[s] * C_P(0) * C_P(1); break;
                  }
                  break;
                    case TWO_ONE:
                      switch( sridx ){
                        case 0: *dRnetdYPtr[s] <<= *dRnetdYPtr[s] - kr * rho * invMsp[s] * 2 * C_P(0) * C_P(1); break;
                        case 1: *dRnetdYPtr[s] <<= *dRnetdYPtr[s] - kr * rho * invMsp[s] * square( C_P(0) );    break;
                      }
                      break;
                        case ONE_TWO:
                          switch( sridx ){
                            case 0: *dRnetdYPtr[s] <<= *dRnetdYPtr[s] - kr * rho * invMsp[s] * square( C_P(1) );    break;
                            case 1: *dRnetdYPtr[s] <<= *dRnetdYPtr[s] - kr * rho * invMsp[s] * 2 * C_P(1) * C_P(0); break;
                          }
                          break;
                            default:
                              nsTmp <<= - kr * rho * invMsp[s];
                              const double pStoich = std::abs( products[sridx].stoich );
                              // we don't need to do anything if pStoich == 1 (within the specified tolerance)
                              if     ( fabs( pStoich - 2 ) < rxnOrderTol ) nsTmp <<= nsTmp * 2       *         C_P(sridx)             ;
                              else if( fabs( pStoich - 3 ) < rxnOrderTol ) nsTmp <<= nsTmp * 3       * square( C_P(sridx)            );
                              else                                         nsTmp <<= nsTmp * pStoich * pow   ( C_P(sridx), pStoich-1 );

                              for( int i = 0; i != products.size(); ++i){
                                if( !( i==sridx ) ){
                                  if     ( fabs( pStoich - 1 ) < rxnOrderTol ) nsTmp <<= nsTmp * C_P(i)                   ;
                                  else if( fabs( pStoich - 2 ) < rxnOrderTol ) nsTmp <<= nsTmp * square( C_P(i)          );
                                  else if( fabs( pStoich - 3 ) < rxnOrderTol ) nsTmp <<= nsTmp * cube  ( C_P(i)          );
                                  else                                         nsTmp <<= nsTmp * pow   ( C_P(i), pStoich );

                                  break;
                                }
                              }
                              *dRnetdYPtr[s] <<= *dRnetdYPtr[s] - nsTmp;
          }
        }
        else{ // nth species is a reactant, so save off its index to prevent a loop below
          nsIsReactant  = true;
          nsReactantIdx = sridx;
        }
      }
      // compute the ns correction term
      if( nsIsReactant ){
        switch( rxnInfo.reverseOrder ){
          case ONE:
            nsTmp <<= kr * rho * invMsp[ns-1];
            break;
          case TWO:
            nsTmp <<= kr * rho * invMsp[ns-1] * 2 * C_P(nsReactantIdx);
            break;
          case ONE_ONE:
            switch( nsReactantIdx ){
              case 0: nsTmp <<= kr * rho * invMsp[ns-1] * C_P(1); break;
              case 1: nsTmp <<= kr * rho * invMsp[ns-1] * C_P(0); break;
            }
            break;
              case ONE_ONE_ONE:
                switch( nsReactantIdx ){
                  case 0: nsTmp <<= kr * rho * invMsp[ns-1] * C_P(1) * C_P(2); break;
                  case 1: nsTmp <<= kr * rho * invMsp[ns-1] * C_P(0) * C_P(2); break;
                  case 2: nsTmp <<= kr * rho * invMsp[ns-1] * C_P(0) * C_P(1); break;
                }
                break;
                  case TWO_ONE:
                    switch( nsReactantIdx ){
                      case 0: nsTmp <<= kr * rho * invMsp[ns-1] * 2 * C_P(0) * C_P(1); break;
                      case 1: nsTmp <<= kr * rho * invMsp[ns-1] * square( C_P(0) );    break;
                    }
                    break;
                      case ONE_TWO:
                        switch( nsReactantIdx ){
                          case 0: nsTmp <<= kr * rho * invMsp[ns-1] * square( C_P(1) );    break;
                          case 1: nsTmp <<= kr * rho * invMsp[ns-1] * 2 * C_P(1) * C_P(0); break;
                        }
                        break;
                          default:
                            nsTmp <<= kr * rho * invMsp[ns-1];
                            const double pStoich = std::abs( products[nsReactantIdx].stoich );
                            // we don't need to do anything if pStoich == 1 (within the specified tolerance)
                            if     ( fabs( pStoich - 2 ) < rxnOrderTol ) nsTmp <<= nsTmp * 2       *         C_P(nsReactantIdx)             ;
                            else if( fabs( pStoich - 3 ) < rxnOrderTol ) nsTmp <<= nsTmp * 3       * square( C_P(nsReactantIdx)            );
                            else                                         nsTmp <<= nsTmp * pStoich * pow   ( C_P(nsReactantIdx), pStoich-1 );
                            for( int i = 0; i != products.size(); ++i){
                              if( !( i==nsReactantIdx ) ){
                                if     ( fabs( pStoich - 1 ) < rxnOrderTol ) nsTmp <<= nsTmp *         C_P(i)           ;
                                else if( fabs( pStoich - 2 ) < rxnOrderTol ) nsTmp <<= nsTmp * square( C_P(i)          );
                                else if( fabs( pStoich - 3 ) < rxnOrderTol ) nsTmp <<= nsTmp * cube  ( C_P(i)          );
                                else                                         nsTmp <<= nsTmp * pow   ( C_P(i), pStoich );
                              }
                            }
        }
        // subtract correction from all species sensitivities
        for( size_t s=0; s<ns-1; ++s ){
          *dRnetdYPtr[s] <<= *dRnetdYPtr[s] + nsTmp;
        }
      }
    }

    const double baseEff  = rxnDat.thdBdyDefault;
    const double* kPCoefs = rxnDat.kPressureCoefs;
    const double* troe    = rxnDat.troeParams;
    switch( rxnDat.type ){
      case ELEMENTARY:
        Ctbaf      <<= 1.0;
        dCtbafdrho <<= 0.0;
        dCtbafdT   <<= 0.0;
        for( size_t s=0; s<ns-1; ++s ){
          *dCtbafdYPtr[s] <<= 0.0;
        }
        break;
      case THIRD_BODY:

        Ctbaf <<= baseEff * ct;
        for( int i=0; i<thdBodies.size(); ++i )
          Ctbaf <<= Ctbaf + rho * THD_BDY(i);

        dCtbafdrho <<= baseEff * invM;
        for( int i=0; i<thdBodies.size(); ++i )
          dCtbafdrho <<= dCtbafdrho + THD_BDY(i);

        dCtbafdT   <<= 0.0;

        for( size_t s=0; s<ns-1; ++s )
          *dCtbafdYPtr[s] <<= rho * baseEff * ( invMsp[s] - invMsp[ns-1] );

        for( int i=0; i<thdBodies.size(); ++i )
          *dCtbafdYPtr[thdBodies[i].index] <<= *dCtbafdYPtr[thdBodies[i].index] + rho * THD_BDY_NOY(i);



        // stage 2 to correct for ns-th species as a third-body
        for( size_t s=0; s<thdBodies.size(); ++s ){
          size_t idx = thdBodies[s].index;
          if( idx == ns-1 ){
            for( size_t ss=0; ss<ns-1; ++ss ){
              *dCtbafdYPtr[ss] <<= *dCtbafdYPtr[ss] - rho * THD_BDY_NOY(s);
            }
          }
        }
        break;
          case LINDEMANN:
            flfConc <<= baseEff * ct;
            for (size_t i = 0; i < thdBodies.size(); ++i)
              flfConc <<= flfConc + rho * THD_BDY(i);

            pr <<= ARRHENIUS(kPCoefs) / kf * flfConc;

            Ctbaf <<= pr / (1. + pr);

            dCtbafdT <<= Ctbaf / (1. + pr) * (ARRHENIUS_SENS_TEMP(kPCoefs) - ARRHENIUS_SENS_TEMP(kCoefs));

            nsTmp <<= Ctbaf / (1. + pr) / flfConc;
            dCtbafdrho <<= nsTmp * baseEff * invM;

            for (size_t s = 0; s < ns; ++s) {
              *dCtbafdYPtr[s] <<= nsTmp * baseEff * rho * invMsp[s];
            }
            for (size_t i = 0; i < thdBodies.size(); ++i){
              dCtbafdrho <<= dCtbafdrho + nsTmp * THD_BDY(i);
              *dCtbafdYPtr[thdBodies[i].index] <<= *dCtbafdYPtr[thdBodies[i].index] + nsTmp * rho * THD_BDY_NOY(i);
            }
            break;

            // stage 2 to correct for ns-th species as a third-body
            for( size_t s=0; s<thdBodies.size(); ++s ){
              size_t idx = thdBodies[s].index;
              if( idx == ns-1 ){
                for( size_t ss=0; ss<ns-1; ++ss ){
                  *dCtbafdYPtr[ss] <<= *dCtbafdYPtr[ss] - nsTmp * rho * ( THD_BDY_NOY(s) );
                }
              }
            }
            break;

              case TROE:
                switch( rxnInfo.troeForm ){
                  case T123:
                    fCent <<= (1-troe[0]) * exp(-T/troe[1]) + troe[0] * exp(-T/troe[2]) + exp(-invT * troe[3]);
                    dfCentdT <<= (troe[0]-1) / troe[1] * exp(-T/troe[1]) - troe[0] / troe[2] * exp(-T/troe[2]) + exp(-invT * troe[3]) * troe[3] * invT * invT;
                    break;
                  case T12:
                    fCent <<= (1-troe[0]) * exp(-T/troe[1]) + troe[0] * exp(-T/troe[2]);
                    dfCentdT <<= (troe[0]-1) / troe[1] * exp(-T/troe[1]) - troe[0] / troe[2] * exp(-T/troe[2]);
                    break;
                  case T1:
                    fCent <<= (1-troe[0]) * exp(-T/troe[1]);
                    dfCentdT <<= (troe[0]-1) / troe[1] * exp(-T/troe[1]);
                    break;
                  case T23:
                    fCent <<= troe[0] * exp(-T/troe[2]) + exp(-invT * troe[3]);
                    dfCentdT <<= - troe[0] / troe[2] * exp(-T/troe[2]) + exp(-invT * troe[3]) * troe[3] * invT * invT;
                    break;
                  case T2:
                    fCent <<= troe[0] * exp(-T/troe[2]);
                    dfCentdT <<= - troe[0] / troe[2] * exp(-T/troe[2]);
                    break;
                  case T13:
                    fCent <<= (1-troe[0]) * exp(-T/troe[1]) + exp(-invT * troe[3]);
                    dfCentdT <<= (troe[0]-1) / troe[1] * exp(-T/troe[1]) + exp(-invT * troe[3]) * troe[3] * invT * invT;
                    break;
                  case T3:
                    fCent <<= exp(-invT * troe[3]);
                    dfCentdT <<= exp(-invT * troe[3]) * troe[3] * invT * invT;
                    break;
                  case NONE:
                  default:{
                    std::ostringstream msg;
                    msg << "Error in " << __FILE__ << " : " << __LINE__ << std::endl
                        <<" Falloff type is TROE, but no terms are flagged for evaluation " << std::endl
                        <<" Reaction # "<< r << std::endl;
                    throw std::runtime_error( msg.str() );
                  }
                }

                flfConc <<= baseEff * ct;
                for( int i=0; i<thdBodies.size(); ++i )
                  flfConc <<= flfConc + rho * THD_BDY(i);

                pr      <<= ARRHENIUS( kPCoefs ) / kf * flfConc;

                aTroe <<=          log10( pr ) -   0.67 * log10( fCent ) - 0.4;
                bTroe <<= - 0.14 * log10( pr ) - 1.1762 * log10( fCent ) + 0.806;
                gTroe <<= 1 / ( 1 + square( aTroe / bTroe ) );

                fTroe <<= pow( fCent, gTroe );
                Ctbaf <<= fTroe * pr / ( 1 + pr );

                dfTroedT <<= fTroe * ( gTroe / fCent * dfCentdT + log( fCent ) *
                             ( -2.0 * gTroe * gTroe / std::log(10) * aTroe / cube( bTroe ) *
                             ( ( bTroe + 0.14 * aTroe ) * ( ARRHENIUS_SENS_TEMP( kPCoefs ) - ARRHENIUS_SENS_TEMP( kCoefs ) ) -
                             ( 0.67 * bTroe - 1.1762 * aTroe ) * dfCentdT / fCent ) ) );
                dCtbafdT <<= 1 / ( 1 + 1 / pr ) * dfTroedT + fTroe * pr / ( ( 1 + pr ) * ( 1 + pr ) ) * ( ARRHENIUS_SENS_TEMP( kPCoefs ) - ARRHENIUS_SENS_TEMP( kCoefs ) );

                nsTmp <<= 1 / flfConc * ( -2.0 / (1 + 1 / pr ) * fTroe * log( fCent ) * gTroe * gTroe / std::log(10) * aTroe / cube( bTroe ) * ( bTroe + 0.14 * aTroe ) +
                          pr * fTroe / ( ( 1 + pr ) * ( 1 + pr ) ) );
                dCtbafdrho <<= nsTmp * baseEff * invM;
                for( int i=0; i<thdBodies.size(); ++i )
                  dCtbafdrho <<= dCtbafdrho + nsTmp * THD_BDY(i);

                for( size_t s=0; s<ns-1; ++s ){
                  *dCtbafdYPtr[s] <<= nsTmp * baseEff * rho * ( invMsp[s] - invMsp[ns-1] );
                }
                for( int i=0; i<thdBodies.size(); ++i )
                  *dCtbafdYPtr[thdBodies[i].index] <<= *dCtbafdYPtr[thdBodies[i].index] + nsTmp * rho * THD_BDY_NOY(i);


                // stage 2 to correct for ns-th species as a third-body
                for( size_t s=0; s<thdBodies.size(); ++s ){
                  size_t idx = thdBodies[s].index;
                  if( idx == ns-1 ){
                    for( size_t ss=0; ss<ns-1; ++ss ){
                      *dCtbafdYPtr[ss] <<= *dCtbafdYPtr[ss] - nsTmp * rho * ( THD_BDY_NOY(s) );
                    }
                  }
                }

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


    q <<= Rnet * Ctbaf;

    dqdrho <<= dRnetdrho * Ctbaf + dCtbafdrho * Rnet;
    dqdT   <<= dRnetdT   * Ctbaf + dCtbafdT   * Rnet;
    for( size_t s=0; s<ns-1; ++s ){
      *dqdYPtr[s] <<= *dRnetdYPtr[s] * Ctbaf + *dCtbafdYPtr[s] * Rnet;
    }

    // rates of formation
    SpecDataVecT::const_iterator iNet = netSpecies.begin();
    for( ; iNet != netSpecies.end(); ++iNet ){
      const int index = iNet->index;

      if( !isInitialized[index] ){
        isInitialized[index] = true;
        productionRates(index)          <<= - ( iNet->stoich * iNet->mw ) * q;
        primitiveSensitivities(index,0) <<= - ( iNet->stoich * iNet->mw ) * dqdrho;
        primitiveSensitivities(index,1) <<= - ( iNet->stoich * iNet->mw ) * dqdT;

        for( size_t s=0; s<ns-1; ++s )
          primitiveSensitivities(index,2+s) <<= - ( iNet->stoich * iNet->mw ) * *dqdYPtr[s];

      }
      else{
        productionRates(index)          <<= productionRates(index)          - ( iNet->stoich * iNet->mw ) * q;
        primitiveSensitivities(index,0) <<= primitiveSensitivities(index,0) - ( iNet->stoich * iNet->mw ) * dqdrho;
        primitiveSensitivities(index,1) <<= primitiveSensitivities(index,1) - ( iNet->stoich * iNet->mw ) * dqdT;

        for( size_t s=0; s<ns-1; ++s )
          primitiveSensitivities(index,2+s) <<= primitiveSensitivities(index,2+s) - ( iNet->stoich * iNet->mw ) * *dqdYPtr[s];

      }
    }


  } // reaction loop

  for( int i=0; i<ns; ++i ){
    if( !isInitialized[i] ){
      productionRates(i) <<= 0.0;
      primitiveSensitivities(i,0) <<= 0.0;
      primitiveSensitivities(i,1) <<= 0.0;
      for( int j=0; j<ns-1; ++j ){
        primitiveSensitivities(i,2+j) <<= 0.0;
      }
    }
  }
}

}; // namespace pokitt
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

#endif // AnalyticalJacobian_h
