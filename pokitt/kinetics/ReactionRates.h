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

#ifndef ReactionRates_Expr_h
#define ReactionRates_Expr_h

#include <expression/Expression.h>
#include <spatialops/util/TimeLogger.h>
#include <expression/FieldRequest.h>

#include <pokitt/CanteraObjects.h> // include cantera wrapper

#include <cantera/kernel/ct_defs.h> // contains value of gas constant
#include <cantera/kernel/reaction_defs.h> // reaction type definitions
#include <cantera/kernel/speciesThermoTypes.h> // contains definitions for which polynomial is being used

#ifdef ENABLE_CUDA
#include <cuda_profiler_api.h>
#include <cuda_runtime_api.h>
#endif

#define PROP_SUM2( prop ) (        prop(0) + prop(1))
#define PROP_SUM3( prop ) (PROP_SUM2(prop) + prop(2))
#define PROP_SUM4( prop ) (PROP_SUM3(prop) + prop(3))
#define PROP_SUM5( prop ) (PROP_SUM4(prop) + prop(4))
#define PROP_SUM6( prop ) (PROP_SUM5(prop) + prop(5))
#define PROP_SUM7( prop ) (PROP_SUM6(prop) + prop(6))
#define PROP_SUM8( prop ) (PROP_SUM7(prop) + prop(7))
#define PROP_SUM9( prop ) (PROP_SUM8(prop) + prop(8))

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

enum PressureEnhancement{
  NO_ENHANCE,
  THIRD_BODY,
  LINDEMANN,
  TROE
};

enum RateConstantForm{
  CONSTANT,
  LINEAR,
  QUADRATIC,
  RECIPROCAL,
  ARRHENIUS
};

enum TroeTerms{
  NONE,
  T1,
  T2,
  T12,
  T3,
  T13,
  T23,
  T123
};

enum ReactionOrder{
  ONE,
  TWO,
  ONE_ONE,
  ONE_ONE_ONE,
  TWO_ONE,
  ONE_TWO,
  OTHER
};

template< typename FieldT >
struct ReactionInfo{
  ReactionInfo< FieldT >(): reversible( false ), kForm( ARRHENIUS ), pFactor( NO_ENHANCE ), troeForm( NONE ),
		                    kRateCoefs(), defaultThdBdyEff( 0 ), thdBdyInd(), thdBdyEff(), numThdBdy( 0 ), thdYi(),
		                    kPressureCoefs(), troeCoefs(), rInd(), rSto(), rMWInv(), rYi(), forwardOrder( OTHER ),
		                    pInd(), pSto(), pMWInv(), pYi(), reverseOrder( OTHER ), nInd(), nSto(), nMW(), netOrder( 0 ), netNumSpec( 0 )

  {}
  typedef std::vector< boost::shared_ptr< const Expr::FieldRequest<FieldT> > > FieldReqVecT;
  bool reversible;
  RateConstantForm kForm;
  PressureEnhancement pFactor;
  TroeTerms troeForm;
  std::vector< double > kRateCoefs;
  double defaultThdBdyEff;
  std::vector< int > thdBdyInd;
  std::vector< double > thdBdyEff;
  FieldReqVecT thdYi;
  int numThdBdy;
  std::vector< double > kPressureCoefs;
  std::vector< double > troeCoefs;

  std::vector< int > rInd;
  std::vector< int > rSto;
  std::vector< double > rMWInv;
  FieldReqVecT rYi;
  ReactionOrder forwardOrder;

  std::vector< int > pInd;
  std::vector< int > pSto;
  std::vector< double > pMWInv;
  FieldReqVecT pYi;
  ReactionOrder reverseOrder;

  std::vector< int > nInd;
  std::vector< int > nSto;
  std::vector< double > nMW;
  int netOrder;
  int netNumSpec;
};

template< typename FieldT >
class ReactionRates
    : public Expr::Expression<FieldT>
{
  typedef std::vector<FieldT*> SpecT;
  typedef std::vector< boost::shared_ptr< const Expr::FieldRequest<FieldT> > > FieldReqVecT;

  DECLARE_FIELDS( FieldT, t_, rho_, mmw_ )
  DECLARE_VECTOR_OF_FIELDS( FieldT, yi_ )

  std::vector< ReactionInfo< FieldT > > rxnInfo_;

  int nSpec_; // number of species in the mechanism
  int nRxns_; // number of reactions in the mechanism
  const double invGasConstant_; // inverse of universal gas constant ( division is expensive)
  bool first_;

  typedef std::vector<double> PolyVals; // values used for polynomial
  PolyVals minTVec_; // vector of minimum temperatures for polynomial evaluations
  PolyVals maxTVec_; // vector of maximum temperatures for polynomial evaluations
  std::vector< PolyVals > cVec_; // vector of polynomial coefficients
  std::vector<int> polyTypeVec_; // vector of polynomial types

  ReactionRates( const Expr::Tag& tTag,
                 const Expr::Tag& pTag,
                 const Expr::TagList& yiTags,
                 const Expr::Tag& mmwTag );
public:
  class Builder : public Expr::ExpressionBuilder
  {
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
             const int nghost = DEFAULT_NUMBER_OF_GHOSTS );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag tTag_;
    const Expr::Tag rhoTag_;
    const Expr::TagList yiTags_;
    const Expr::Tag mmwTag_;
  };

  ~ReactionRates(){}
  void evaluate();
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
               const Expr::Tag& mmwTag )
  : Expr::Expression<FieldT>(),
    invGasConstant_( 1 / Cantera::GasConstant )
{
  first_ = true;
  this->set_gpu_runnable( true );

  t_   = this->template create_field_request<FieldT>( tTag   );
  rho_   = this->template create_field_request<FieldT>( rhoTag   );
  mmw_ = this->template create_field_request<FieldT>( mmwTag );
  this->template create_field_vector_request<FieldT>( yiTags, yi_ );

  using namespace Cantera;
  Cantera_CXX::IdealGasMix* gasMix = CanteraObjects::get_gasmix();
  const std::vector<ReactionData>& canteraDataVec = gasMix->getReactionData(); // contains kinetics data for each reaction
  nSpec_ = gasMix->nSpecies();
  nRxns_ = canteraDataVec.size();
  std::vector< double > molecularWeights = gasMix->molecularWeights();
  std::vector< double > molecularWeightsInv(nSpec_, 0.0);
  for( size_t n=0; n<nSpec_; ++n ){
	  molecularWeightsInv[n] = 1 / molecularWeights[n];
  }

  rxnInfo_.resize( nRxns_, ReactionInfo<FieldT>() );
  for( size_t r=0; r<nRxns_; ++r){
    const Cantera::ReactionData& canteraData = canteraDataVec[r]; // ReactionData for reaction r
    ReactionInfo<FieldT>& rxnInfo = rxnInfo_[r];

    rxnInfo.kRateCoefs       = canteraData.rateCoeffParameters;
    rxnInfo.kPressureCoefs   = canteraData.auxRateCoeffParameters;
    rxnInfo.defaultThdBdyEff = canteraData.default_3b_eff;
    rxnInfo.troeCoefs        = canteraData.falloffParameters;
    rxnInfo.reversible       = canteraData.reversible;

    /* here we are reading Cantera's product and reactant stoichiometric coefficients and storing them as ints
     * This speeds up evaluating the product (multiplication) of concentrations when evaluating reaction rates of progress
     * only elementary reactions are supported
     */
    std::vector< double >::const_iterator iCSto;
    std::vector< int >::const_iterator iCInd = canteraData.reactants.begin();
    for( iCSto = canteraData.rstoich.begin(); iCSto!=canteraData.rstoich.end(); ++iCSto, ++iCInd){
      rxnInfo.rInd.push_back( *iCInd );
      if(      fabs( *iCSto-1 ) < 1e-2 ) rxnInfo.rSto.push_back( 1 );
      else if( fabs( *iCSto-2 ) < 1e-2 ) rxnInfo.rSto.push_back( 2 );
      else if( fabs( *iCSto-3 ) < 1e-2 ) rxnInfo.rSto.push_back( 3 );
      else{
        std::ostringstream msg;
        msg << "Error in " __FILE__ << " : " << __LINE__ << std::endl
            <<" Non-integer reactant stoichiometric coefficient" << std::endl
            <<" Reaction # "<< r <<", stoichiometric coefficient = " << *iCSto << std::endl;
        throw std::runtime_error( msg.str() );
      }
    }
    if( rxnInfo.rInd.size() < 1 || rxnInfo.rInd.size() > 3 ){
      std::ostringstream msg;
      msg << "Error in " __FILE__ << " : " << __LINE__ << std::endl
          <<" Number of reactants must be 1 <= n <= 3 " << std::endl
          <<" Reaction # "<< r <<", number of reactants n = " << rxnInfo.rInd.size() << std::endl;
      throw std::runtime_error( msg.str() );
    }
    switch( rxnInfo.rSto.size() ){
    case 1:{
      if( rxnInfo.rSto[0] == 1 ) rxnInfo.forwardOrder = ONE;
      else if( rxnInfo.rSto[0] == 2 ) rxnInfo.forwardOrder = TWO;
      else {
        rxnInfo.forwardOrder = OTHER;
        std::cout << r << " OTHER forward reaction\n";
      }
      break;
    }
    case 2:{
      if( rxnInfo.rSto[0] == 1 && rxnInfo.rSto[1] == 1 ) rxnInfo.forwardOrder = ONE_ONE;
      else if( rxnInfo.rSto[0] == 1 && rxnInfo.rSto[1] == 2 ) rxnInfo.forwardOrder = ONE_TWO;
      else if( rxnInfo.rSto[0] == 2 && rxnInfo.rSto[1] == 1 ) rxnInfo.forwardOrder = TWO_ONE;
      else{
        rxnInfo.forwardOrder = OTHER;
        std::cout << r << " OTHER forward reaction\n";
      }
      break;
    }
    case 3:{
      if(rxnInfo.rSto[0] == 1 && rxnInfo.rSto[1] == 1 && rxnInfo.rSto[2] == 1 ){
        rxnInfo.forwardOrder = ONE_ONE_ONE;
      }
      else{
        rxnInfo.forwardOrder = OTHER;
        std::cout << r << " OTHER forward reaction\n";
      }
      break;
    }
    default:{
      rxnInfo.reverseOrder = OTHER;
      std::cout << r << " OTHER forward reaction\n";
    }
    }

    iCInd = canteraData.products.begin();
    for( iCSto = canteraData.pstoich.begin(); iCSto!=canteraData.pstoich.end(); ++iCSto, ++iCInd ){
      rxnInfo.pInd.push_back( *iCInd );
      if(      fabs( *iCSto - 1 ) < 1e-2 ) rxnInfo.pSto.push_back( -1 );
      else if( fabs( *iCSto - 2 ) < 1e-2 ) rxnInfo.pSto.push_back( -2 );
      else if( fabs( *iCSto - 3 ) < 1e-2 ) rxnInfo.pSto.push_back( -3 );
      else{
        std::ostringstream msg;
        msg << "Error in " __FILE__ << " : " << __LINE__ << std::endl
            <<" Non-integer product stoichiometric coefficient" << std::endl
            <<" Reaction # "<< r <<", stoichiometric coefficient = " << *iCSto << std::endl;
        throw std::runtime_error( msg.str() );
      }
    }
    if( rxnInfo.pInd.size() < 1 || rxnInfo.pInd.size() > 3 ){
      std::ostringstream msg;
      msg << "Error in " __FILE__ << " : " << __LINE__ << std::endl
          <<" Number of products must be 1 <= n <= 3 " << std::endl
          <<" Reaction # "<< r <<", number of products n = " << rxnInfo.pInd.size() << std::endl;
      throw std::runtime_error( msg.str() );
    }
    switch( rxnInfo.pSto.size() ){
    case 1:{
      if( rxnInfo.pSto[0] == -1 ) rxnInfo.reverseOrder = ONE;
      else if( rxnInfo.pSto[0] == -2 ) rxnInfo.reverseOrder = TWO;
      else {
        rxnInfo.reverseOrder = OTHER;
        std::cout << r << " OTHER reverse reaction\n";
      }
      break;
    }
    case 2:{
      if( rxnInfo.pSto[0] == -1 && rxnInfo.pSto[1] == -1 ) rxnInfo.reverseOrder = ONE_ONE;
      else if( rxnInfo.pSto[0] == -1 && rxnInfo.pSto[1] == -2 ) rxnInfo.reverseOrder = ONE_TWO;
      else if( rxnInfo.pSto[0] == -2 && rxnInfo.pSto[1] == -1 ) rxnInfo.reverseOrder = TWO_ONE;
      else{
        rxnInfo.reverseOrder = OTHER;
        std::cout << r << " OTHER reaction\n";
      }
      break;
    }
    case 3:{
      if(rxnInfo.pSto[0] == -1 && rxnInfo.pSto[1] == -1 && rxnInfo.pSto[2] == -1 ){
        rxnInfo.reverseOrder = ONE_ONE_ONE;
      }
      else{
        rxnInfo.reverseOrder = OTHER;
        std::cout << r << " OTHER reaction\n";
      }
      break;
    }
    default:{
      rxnInfo.reverseOrder = OTHER;
      std::cout << r << " OTHER reaction\n";
    }
    }

    rxnInfo.netOrder = 0; // difference of forward and reverse rxn orders, used for equilibrium K
    std::map< int, int > netSpecies; // start with a map to avoid duplicating participating species which aren't reacting
    std::vector< int >::iterator iInd = rxnInfo.rInd.begin();
    std::vector< int >::iterator iSto = rxnInfo.rSto.begin();
    for(  ; iSto != rxnInfo.rSto.end(); ++iSto, ++iInd ){
      netSpecies.insert( std::make_pair( *iInd, *iSto ) );
      rxnInfo.netOrder += *iSto;
    }

    iInd = rxnInfo.pInd.begin();
    for( iSto = rxnInfo.pSto.begin(); iSto != rxnInfo.pSto.end(); ++iSto, ++iInd ){
      if( netSpecies.find( *iInd ) != netSpecies.end() )
        netSpecies[ *iInd ] += *iSto;
      else
        netSpecies.insert( std::make_pair( *iInd, *iSto ) );
      rxnInfo.netOrder += *iSto;
    }

    std::map< int, int >::iterator iNet; // now we keep the net species
    for( iNet = netSpecies.begin(); iNet != netSpecies.end(); ++iNet ){
      if( iNet->second != 0 ){
        rxnInfo.nInd.push_back( iNet->first );
        rxnInfo.nSto.push_back( iNet->second );
      }
    }

    rxnInfo.netNumSpec = rxnInfo.nInd.size(); // doesn't catch reversibles
    if( rxnInfo.netNumSpec < 2 || rxnInfo.netNumSpec > 5 ){
      std::ostringstream msg;
      msg << "Error in " __FILE__ << " : " << __LINE__ << std::endl
          <<" The number of reacting species must be 2 <= n <= 5 " << std::endl
          <<" Reaction # "<< r <<", number of active species n = " << rxnInfo.netNumSpec << std::endl;
      throw std::runtime_error( msg.str() );
    }

    for( iInd = rxnInfo.rInd.begin(); iInd != rxnInfo.rInd.end(); ++iInd){
      rxnInfo.rMWInv.push_back( molecularWeightsInv[*iInd] );
      rxnInfo.rYi.push_back( yi_[*iInd] );
    }

    for( iInd = rxnInfo.pInd.begin(); iInd != rxnInfo.pInd.end(); ++iInd){
      rxnInfo.pMWInv.push_back( molecularWeightsInv[*iInd] );
      rxnInfo.pYi.push_back( yi_[*iInd] );
    }

    for( iInd = rxnInfo.nInd.begin(); iInd != rxnInfo.nInd.end(); ++iInd){
      rxnInfo.nMW.push_back( molecularWeights[*iInd] );
    }

    /*
     * Here we are checking if the rate constant is Arrhenius in form
     * If it is, we need to perform an exponential evaluation (expensive)
     * If not, then we determine which power of temperature to use
     */
    const std::vector<double>& arrhenius = rxnInfo.kRateCoefs;
    if( fabs( arrhenius[2] ) < 1e-6 ){ // i.e. 0 activation energy
      if( fabs( arrhenius[1] ) < 1e-6 )
        rxnInfo.kForm = CONSTANT;
      else if( fabs( arrhenius[1] - 1 ) < 1e-6 )
        rxnInfo.kForm = LINEAR;
      else if( fabs( arrhenius[1] - 2 ) < 1e-6 )
        rxnInfo.kForm = QUADRATIC;
      else if( fabs( arrhenius[1] + 1 ) < 1e-6 )
        rxnInfo.kForm = RECIPROCAL;
      else
        rxnInfo.kForm = ARRHENIUS; // the exponent for temperature is not and integer
    }
    else
      rxnInfo.kForm = ARRHENIUS;

    /* Here we are converting the efficiency map into a vector of pairs.
     * This should perform slightly better and is the same type as
     * other containers in the struct
     */
    std::map<int, double>::const_iterator iEff;
    std::map<int, double>::const_iterator iEffE = canteraData.thirdBodyEfficiencies.end();
    for( iEff  = canteraData.thirdBodyEfficiencies.begin(); iEff!=iEffE ; ++iEff){
        rxnInfo.thdBdyInd.push_back( iEff->first );
        rxnInfo.thdBdyEff.push_back( iEff->second );
        rxnInfo.thdBdyEff.back() = molecularWeightsInv[iEff->first] * ( rxnInfo.thdBdyEff.back() - rxnInfo.defaultThdBdyEff );
        rxnInfo.thdYi.push_back( yi_[iEff->first]);
    }
    rxnInfo.numThdBdy = rxnInfo.thdBdyInd.size();


    switch( canteraData.reactionType ){
    case ELEMENTARY_RXN: rxnInfo.pFactor = NO_ENHANCE; break;
    case THREE_BODY_RXN: rxnInfo.pFactor = THIRD_BODY; break;
    case FALLOFF_RXN:
      switch( canteraData.falloffType ){
      case SIMPLE_FALLOFF: rxnInfo.pFactor = LINDEMANN; break;
      case TROE3_FALLOFF:
        if( fabs( canteraData.falloffParameters[1] ) < 1e-8 && fabs( canteraData.falloffParameters[2] ) < 1e-8 ){
          rxnInfo.pFactor  = LINDEMANN;
          rxnInfo.troeForm = NONE;
        }
        else{
        rxnInfo.pFactor = TROE;
        if(      fabs( canteraData.falloffParameters[1] ) < 1e-8 )
          rxnInfo.troeForm = T2;
        else if( fabs( canteraData.falloffParameters[2] ) < 1e-8 )
          rxnInfo.troeForm = T1;
        else
          rxnInfo.troeForm = T12;
        }
        break;
      case TROE4_FALLOFF:
        rxnInfo.pFactor = TROE;
        if( fabs( canteraData.falloffParameters[1] ) < 1e-8 && fabs( canteraData.falloffParameters[2] ) < 1e-8 )
          rxnInfo.troeForm = T3;
        else if( fabs( canteraData.falloffParameters[1] ) < 1e-8 )
          rxnInfo.troeForm = T23;
        else if( fabs( canteraData.falloffParameters[2] ) < 1e-8 )
          rxnInfo.troeForm = T13;
        else
          rxnInfo.troeForm = T123;
        break;
      }
    }

  }

  // collect polynomial coefficients for evaluating the gibbs energy function
  const Cantera::SpeciesThermo& spThermo = gasMix->speciesThermo();
  std::vector<double> c(15,0); //vector of Cantera's coefficients
  int polyType;
  double minT;
  double maxT;
  double refPressure;

  for( size_t n=0; n<nSpec_; ++n ){
    spThermo.reportParams(n, polyType, &c[0], minT, maxT, refPressure);
    polyTypeVec_.push_back(polyType); // vector of polynomial types
    minTVec_.push_back(minT); // vector of minimum temperatures for polynomial evaluations
    maxTVec_.push_back(maxT); // vector of maximum temperatures for polynomial evaluations
    switch (polyType) {
    case SIMPLE: break;
    case NASA2:
      for( std::vector<double>::iterator ic = c.begin() + 1; ic!=c.end(); ++ic)
        *ic *= Cantera::GasConstant; // dimensionalize the coefficients
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
          <<" Thermo type not supported, type = " << polyType << std::endl
          <<" See <cantera/kernel/SpeciesThermoInterpType.h> for type definition" << std::endl;
      throw std::runtime_error( msg.str() );
    }
    }
    cVec_.push_back(c); // vector of polynomial coefficients
  }

  CanteraObjects::restore_gasmix(gasMix);
}

//--------------------------------------------------------------------

template< typename FieldT >
void
ReactionRates<FieldT>::
evaluate()
{
  using namespace SpatialOps;
#ifdef ENABLE_CUDA
  if( !first_ )cudaProfilerStart();
#endif

  SpecT& rRates = this->get_value_vec();
  const FieldT& t   = t_->field_ref();
  const FieldT& rho = rho_->field_ref();
  const FieldT& mmw = mmw_->field_ref();

  std::vector< bool > initialized ( nSpec_, false); // avoids needing to assign ri <<= 0

  // logs and reciprocals to save time during the loop
  SpatFldPtr<FieldT> concPtr     = SpatialFieldStore::get<FieldT>( t ); // molar concentration
  SpatFldPtr<FieldT> logConcPtr  = SpatialFieldStore::get<FieldT>( t ); // log of molar concentration
  SpatFldPtr<FieldT> logTPtr     = SpatialFieldStore::get<FieldT>( t ); // log(t)
  SpatFldPtr<FieldT> tRecipPtr   = SpatialFieldStore::get<FieldT>( t ); // 1/t

  FieldT& conc    = *concPtr;
  FieldT& logConc = *logConcPtr;
  FieldT& logT   = *logTPtr;
  FieldT& tRecip  = *tRecipPtr;

  conc    <<= rho / mmw;
  logConc <<= log( conc );
  logT    <<= log(t);
  tRecip  <<= 1/t;

  // calculate the delta gibbs energy for each species for use in evaluating reversible rate constants
  std::vector< SpatFldPtr<FieldT> > specG; // gibbs energy for each species
  for( int n=0; n<nSpec_; ++n ){
    specG.push_back( SpatialFieldStore::get<FieldT>(t) );
    const int polyType = polyTypeVec_[n];
    const std::vector<double>& c = cVec_[n];
    /* polynomials are applicable in two temperature ranges - high and low
     * Caution: polynomials are used even if temperature is out of range
     * Note: coefficients have been divided by integers during construction of expression
     */
    switch ( polyType ){
    case NASA2:
      *specG[n] <<= cond( t <= c[0], c[ 6] + t * ( c[1] - c[ 7] - c[1] * logT - t * ( c[2] + t * ( c[ 3] + t * ( c[ 4] + t * c[ 5] )))) )  // if low temp
                        (            c[13] + t * ( c[8] - c[14] - c[8] * logT - t * ( c[9] + t * ( c[10] + t * ( c[11] + t * c[12] )))) );  // else if high temp
      break;
    case SIMPLE:
      *specG[n] <<=       c[1] + c[3] * ( t    - c[0]        ) // H
                  - t * ( c[2] + c[3] * ( logT - log(c[0]) ) ); // -TS
      break;
    } // switch
  } // species loop

  // temporary fields for forward rate constant and falloff reactions
  SpatFldPtr<FieldT> kPtr        = SpatialFieldStore::get<FieldT>( t ); // forward rate constant
  SpatFldPtr<FieldT> mPtr        = SpatialFieldStore::get<FieldT>( t ); // third body enhancement factor
  SpatFldPtr<FieldT> prPtr       = SpatialFieldStore::get<FieldT>( t ); // normalized pressure factor
  SpatFldPtr<FieldT> logFCentPtr = SpatialFieldStore::get<FieldT>( t ); // log10(Fcent)
  SpatFldPtr<FieldT> logPrCPtr   = SpatialFieldStore::get<FieldT>( t ); // log10(Pr) + CTROE

  FieldT& k  = *kPtr;
  FieldT& m  = *mPtr;
  FieldT& pr = *prPtr;
  FieldT& logPrC = *logPrCPtr;
  FieldT& logFCent = *logFCentPtr;

  // fields for the reverse rates
  SpatFldPtr<FieldT> krPtr   = SpatialFieldStore::get<FieldT>( t ); // equilibrium constant
  FieldT& kr = *krPtr;

  for( int r=0; r<nRxns_; ++r ){
    const ReactionInfo< FieldT >& rxnInfo = rxnInfo_[r];

    const std::vector<double>& rateCoefs = rxnInfo.kRateCoefs;
#define ARRHENIUS(A,b,E) ((A) * exp( (b) * logT - (E) * tRecip )) // convenience definition
    switch( rxnInfo.kForm ){
    case CONSTANT:   k <<= rateCoefs[0];          break;
    case LINEAR:     k <<= rateCoefs[0] * t;      break;
    case QUADRATIC:  k <<= rateCoefs[0] * t * t;  break;
    case RECIPROCAL: k <<= rateCoefs[0] * tRecip; break;
    case ARRHENIUS:  k <<= ARRHENIUS( rateCoefs[0], rateCoefs[1], rateCoefs[2] ); break;
    }

    const double baseEff             = rxnInfo.defaultThdBdyEff;
    const std::vector< double >& thdEff = rxnInfo.thdBdyEff;
    const FieldReqVecT           thdYi  = rxnInfo.thdYi;
#   define THD_BDY(i) ( thdEff[(i)] * thdYi[(i)]->field_ref() ) // represents non-default third bodies
    switch( rxnInfo.pFactor ){
    case NO_ENHANCE: break;
    case THIRD_BODY:
      switch( rxnInfo.numThdBdy ){
      case 0:
        k <<= k * baseEff * conc;
        break;
      case 1:
        k <<= k * ( baseEff * conc + rho * ( THD_BDY(0) ) );
        break;
      case 2:
        k <<= k * ( baseEff * conc + rho * ( PROP_SUM2( THD_BDY ) ) );
        break;
      case 3:
        k <<= k * ( baseEff * conc + rho * ( PROP_SUM3( THD_BDY ) ) );
        break;
      case 4:
        k <<= k * ( baseEff * conc + rho * ( PROP_SUM4( THD_BDY ) ) );
        break;
      case 5:
        k <<= k * ( baseEff * conc + rho * ( PROP_SUM5( THD_BDY ) ) );
        break;
      case 6:
        k <<= k * ( baseEff * conc + rho * ( PROP_SUM6( THD_BDY ) ) );
        break;
      case 7:
        k <<= k * ( baseEff * conc + rho * ( PROP_SUM7( THD_BDY ) ) );
        break;
      case 8:
        k <<= k * ( baseEff * conc + rho * ( PROP_SUM8( THD_BDY ) ) );
        break;
      default:{
        m <<= baseEff * conc;
        for( int i = 0; i != rxnInfo.numThdBdy; ++i )
          m <<= m + THD_BDY(i) * rho;
        k <<= k * m;
        break;
        }
      }
    break;
    case LINDEMANN:{
      const double A = rxnInfo.kPressureCoefs[0];
      const double b = rxnInfo.kPressureCoefs[1];
      const double E = rxnInfo.kPressureCoefs[2];
        switch( rxnInfo.numThdBdy ){
        case 0:
          k <<= k / ( 1 + k / ( ARRHENIUS(A,b,E) *   baseEff * conc ) );
          break;
        case 1:
          k <<= k / ( 1 + k / ( ARRHENIUS(A,b,E) * ( baseEff * conc + rho * ( THD_BDY(0) ) ) ) );
          break;
        case 2:
          k <<= k / ( 1 + k / ( ARRHENIUS(A,b,E) * ( baseEff * conc + rho * ( PROP_SUM2( THD_BDY )  ) ) ) );
          break;
        case 3:
          k <<= k / ( 1 + k / ( ARRHENIUS(A,b,E) * ( baseEff * conc + rho * ( PROP_SUM3( THD_BDY )  ) ) ) );
          break;
        case 4:
          k <<= k / ( 1 + k / ( ARRHENIUS(A,b,E) * ( baseEff * conc + rho * ( PROP_SUM4( THD_BDY ) ) ) ) );
          break;
        case 5:
          k <<= k / ( 1 + k / ( ARRHENIUS(A,b,E) * ( baseEff * conc + rho * ( PROP_SUM5( THD_BDY )  ) ) ) );
          break;
        case 6:
          k <<= k / ( 1 + k / ( ARRHENIUS(A,b,E) * ( baseEff * conc + rho * ( PROP_SUM6( THD_BDY ) ) ) ) );
          break;
        case 7:
          k <<= k / ( 1 + k / ( ARRHENIUS(A,b,E) * ( baseEff * conc + rho * ( PROP_SUM7( THD_BDY ) ) ) ) );
          break;
        case 8:
          k <<= k / ( 1 + k / ( ARRHENIUS(A,b,E) * ( baseEff * conc + rho * ( PROP_SUM8( THD_BDY ) ) ) ) );
          break;
        default:{
          m <<= baseEff * conc;
          for( int i = 0; i != rxnInfo.numThdBdy; ++i )
            m <<= m + THD_BDY(i) * rho;
          k <<= k / ( 1 + k / ( ARRHENIUS(A,b,E) * m ) );
          break;
          }
        }
        break;
    }
      case TROE:{
        const double A = rxnInfo.kPressureCoefs[0];
        const double b = rxnInfo.kPressureCoefs[1];
        const double E = rxnInfo.kPressureCoefs[2];
        switch( rxnInfo.numThdBdy ){
              case 0:
                pr <<= ARRHENIUS( A, b, E ) / k * ( baseEff * conc );
                break;
              case 1:
                pr <<= ARRHENIUS( A, b, E ) / k * ( baseEff * conc + rho * THD_BDY(0) );
                break;
              case 2:
                pr <<= ARRHENIUS( A, b, E ) / k * ( baseEff * conc + rho * ( PROP_SUM2( THD_BDY ) ) );
                break;
              case 3:
                pr <<= ARRHENIUS( A, b, E ) / k * ( baseEff * conc + rho * ( PROP_SUM3( THD_BDY ) ) );
                break;
              case 4:
                pr <<= ARRHENIUS( A, b, E ) / k * ( baseEff * conc + rho * ( PROP_SUM4( THD_BDY ) ) );
                break;
              case 5:
                pr <<= ARRHENIUS( A, b, E ) / k * ( baseEff * conc + rho * ( PROP_SUM5( THD_BDY ) ) );
                break;
              case 6:
                pr <<= ARRHENIUS( A, b, E ) / k * ( baseEff * conc + rho * ( PROP_SUM6( THD_BDY ) ) );
                break;
              case 7:
                pr <<= ARRHENIUS( A, b, E ) / k * ( baseEff * conc + rho * ( PROP_SUM7( THD_BDY ) ) );
                break;
              case 8:
                pr <<= ARRHENIUS( A, b, E ) / k * ( baseEff * conc + rho * ( PROP_SUM8( THD_BDY ) ) );
                break;
              default:{
                m <<= baseEff * conc;
                for( int i = 0; i != rxnInfo.numThdBdy; ++i )
                  m <<= m + THD_BDY(i) * rho;
                pr <<= ARRHENIUS( A, b, E ) / k * m;
                break;
                }
              }
#undef ARRHENIUS
#   undef THD_BDY
        const std::vector<double>&      troe = rxnInfo.troeCoefs;
        switch( rxnInfo.troeForm ){
        case T1:   logFCent <<= log10( (1-troe[0]) * exp(-t/troe[1])                                                      ); break;
        case T2:   logFCent <<= log10(                                 troe[0] * exp(-t/troe[2])                          ); break;
        case T12:  logFCent <<= log10( (1-troe[0]) * exp(-t/troe[1]) + troe[0] * exp(-t/troe[2])                          ); break;
        case T3:   logFCent <<= log10(                                                             exp(-tRecip * troe[3]) ); break;
        case T13:  logFCent <<= log10( (1-troe[0]) * exp(-t/troe[1])                             + exp(-tRecip * troe[3]) ); break;
        case T23:  logFCent <<= log10(                                 troe[0] * exp(-t/troe[2]) + exp(-tRecip * troe[3]) ); break;
        case T123: logFCent <<= log10( (1-troe[0]) * exp(-t/troe[1]) + troe[0] * exp(-t/troe[2]) + exp(-tRecip * troe[3]) ); break;
        case NONE:{
          std::ostringstream msg;
          msg << "Error in " << __FILE__ << " : " << __LINE__ << std::endl
              <<" Falloff type is TROE, but no terms are flagged for evaluation " << std::endl
              <<" Reaction # "<< r << std::endl;
          throw std::runtime_error( msg.str() );
        }
        }

#     define CTROE ( -0.4 - 0.67 * logFCent ) // macros to make f1 evaluation easier to read
#     define NTROE ( 0.75 - 1.27 * logFCent ) // macros to make k evaluation easier to read
#     define F1 ( logPrC / ( NTROE - 0.14 * logPrC ) )

        logPrC <<= log10( pr ) + CTROE;
        k <<= k * pow(10, logFCent / (1 + square( F1 ) ) ) * pr / (1+pr);

#undef CTROE
#undef NTROE
#undef F1
        break;
      }
    }

    const std::vector< int >& nInd = rxnInfo.nInd;
    const std::vector< int >& nSto = rxnInfo.nSto;
#   define GIBBS(i) ( nSto[(i)] * *specG[ nInd[(i)] ] ) // more legible code
    if( rxnInfo.reversible ){
      const int netOrder = rxnInfo.netOrder;
      switch( rxnInfo.netNumSpec ){
      case 3:
        kr <<= k * exp( netOrder * logConc - tRecip * invGasConstant_ * ( PROP_SUM3( GIBBS )                      ) );
        break;
      case 2:
        kr <<= k * exp( netOrder * logConc - tRecip * invGasConstant_ * ( PROP_SUM2( GIBBS )                                  ) );
        break;
      case 4:
        kr <<= k * exp( netOrder * logConc - tRecip * invGasConstant_ * ( PROP_SUM4( GIBBS )            ) );
        break;
      case 5:
        kr <<= k * exp( netOrder * logConc - tRecip * invGasConstant_ * ( PROP_SUM5( GIBBS ) ) );
        break;
      }
    }
#   undef GIBBS



    const std::vector< double >& rMWInv = rxnInfo.rMWInv;
    const FieldReqVecT& rYi = rxnInfo.rYi;
#define C_R(i) ( rYi[(i)]->field_ref() * rho * rMWInv[(i)] )
    switch( rxnInfo.forwardOrder ){
    case ONE:
      k <<= k  *  C_R(0);
      break;
    case TWO:
      k <<= k  * square( C_R(0) );
      break;
    case ONE_ONE:
      k <<= k  * C_R(0) * C_R(1);
      break;
    case ONE_ONE_ONE:
      k <<= k  * C_R(0) * C_R(1) * C_R(2);
      break;
    case TWO_ONE:
      k <<= k  * square( C_R(0) ) * C_R(1) ;
      break;
    case ONE_TWO:
      k <<= k  * C_R(0) * square( C_R(1) );
      break;
    default:{
      std::vector< int >::const_iterator iSto  = rxnInfo.rSto.begin();
      std::vector< int >::const_iterator iStoE = rxnInfo.rSto.end();
      std::vector< double >::const_iterator iMWInv = rxnInfo.rMWInv.begin();
      typename FieldReqVecT::const_iterator iYi = rYi.begin();
      for( ; iSto != iStoE; ++iSto, ++iYi, ++iMWInv ){
        switch( *iSto ){
        case 1: k <<= k * (*iYi)->field_ref() * *iMWInv * rho; break;
        case 2: k <<= k * square( (*iYi)->field_ref() * *iMWInv * rho ); break;
        case 3: k <<= k * cube( (*iYi)->field_ref() * *iMWInv * rho ); break;
        }
      }
      break;
    }
    }
#undef C_R

    const std::vector< double >& pMWInv = rxnInfo.pMWInv;
    const FieldReqVecT& pYi = rxnInfo.pYi;
#define C_P(i) ( pYi[(i)]->field_ref() * rho * pMWInv[(i)] )
    if( rxnInfo.reversible ){
      switch( rxnInfo.reverseOrder ){
      case ONE:
        kr <<= kr  * C_P(0);  break;
      case TWO:
        kr <<= kr  * square( C_P(0) ); break;
      case ONE_ONE:
        kr <<= kr  *  C_P(0) * C_P(1); break;
      case ONE_ONE_ONE:
        kr <<= kr  * C_P(0) * C_P(1) * C_P(2); break;
      case TWO_ONE:
        kr <<= kr  * square( C_P(0) ) * C_P(1) ; break;
      case ONE_TWO:
        kr <<= kr  * C_P(0) * square( C_P(1) ); break;
      default:{
        std::vector< int >::const_iterator iSto  = rxnInfo.pSto.begin();
        std::vector< int >::const_iterator iStoE = rxnInfo.pSto.end();
        std::vector< double >::const_iterator iMWInv = rxnInfo.pMWInv.begin();
        typename FieldReqVecT::const_iterator iYi = pYi.begin();
        for( ; iSto != iStoE; ++iSto, ++iYi, ++iMWInv ){
          switch( *iSto ){
          case 1: kr <<= kr * (*iYi)->field_ref() * *iMWInv * rho; break;
          case 2: kr <<= kr * square( (*iYi)->field_ref() * *iMWInv * rho ); break;
          case 3: kr <<= kr * cube( (*iYi)->field_ref() * *iMWInv * rho ); break;
          }
        }
        break;
      }
    }
    }
#undef C_P
    std::vector< int >::const_iterator iInd = rxnInfo.nInd.begin();
    std::vector< int >::const_iterator iIndE = rxnInfo.nInd.end();
     std::vector< double >::const_iterator iMW = rxnInfo.nMW.begin();
     std::vector< int >::const_iterator       iSto = rxnInfo.nSto.begin();
    for( ; iInd != iIndE; ++iInd, ++iMW, ++iSto ){
      FieldT& ri = *rRates[*iInd];
      if( !initialized[*iInd] ){
        initialized[*iInd] = true;
        if( rxnInfo.reversible ) ri <<=    - ( *iSto * *iMW ) * ( k - kr );
        else                     ri <<=    - ( *iSto * *iMW ) * ( k      );
      }
      else{
        if( rxnInfo.reversible ) ri <<= ri - ( *iSto * *iMW ) * ( k - kr );
        else                     ri <<= ri - ( *iSto * *iMW ) * ( k      );
      }
    }
  } // loop over reactions

    for( int n = 0; n < nSpec_; ++n ){
       if( !initialized[n] )*rRates[n] <<= 0.0; //in case of inerts that have 0.0 rxn rates
    }
#ifdef ENABLE_CUDA
  if( !first_ ) cudaProfilerStop();
#endif
  first_ = false;
}

//--------------------------------------------------------------------

template< typename FieldT >
ReactionRates<FieldT>::
Builder::Builder( const Expr::TagList& resultTags,
                  const Expr::Tag& tTag,
                  const Expr::Tag& rhoTag,
                  const Expr::TagList& yiTags,
                  const Expr::Tag& mmwTag,
                  const int nghost )
: ExpressionBuilder( resultTags, nghost ),
  tTag_( tTag ),
  rhoTag_( rhoTag ),
  yiTags_( yiTags ),
  mmwTag_( mmwTag )
{}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
ReactionRates<FieldT>::
Builder::build() const
{
  return new ReactionRates<FieldT>( tTag_, rhoTag_, yiTags_, mmwTag_ );
}

} // namespace pokitt

#endif // ReactionRates_Expr_h
