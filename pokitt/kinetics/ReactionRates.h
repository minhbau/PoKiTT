#ifndef ReactionRates_Expr_h
#define ReactionRates_Expr_h

#include <expression/Expression.h>

#include <pokitt/CanteraObjects.h> // include cantera wrapper
#include <pokitt/thermo/Temperature.h>

#include <cantera/kernel/ct_defs.h> // contains value of gas constant
#include <cantera/kernel/reaction_defs.h> // reaction type definitions
#include <cantera/kernel/speciesThermoTypes.h> // contains definitions for which polynomial is being used


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
  typedef std::vector<FieldT*> SpecT;

  const Expr::Tag tTag_;
  const Expr::TagList tPowerTags_;
  const Expr::Tag pTag_;
  const Expr::Tag mmwTag_;
  const Expr::TagList massFracTags_;

  const FieldT* t_;
  std::vector<const FieldT*> tPowers_;
  const FieldT* p_;
  const FieldT* mmw_; // mixture molecular weight

  std::vector<const FieldT*> massFracs_;

  std::vector<Cantera::ReactionData> rxnDataVec_;

  int nSpec_; // number of species in the mechanism
  int nRxns_; // number of reactions in the mechanism

  std::vector<int> dOrder_; // difference between the forward reaction order and backwards reaction order
  std::vector<std::vector<int> > rStoich_; // reactant stoich coeffs, stores ints instead of the doubles from Cantera
  std::vector<std::vector<int> > pStoich_; // product stoich coeffs, stores ints instead of the doubles from Cantera

  std::vector< double > molecularWeights_; // molecular weights
  std::vector< double > molecularWeightsInv_; // inverse of molecular weights (dividing by MW is expensive)

  typedef std::vector<double> PolyVals; // values used for polynomial
  PolyVals minTVec_; // vector of minimum temperatures for polynomial evaluations
  PolyVals maxTVec_; // vector of maximum temperatures for polynomial evaluations
  std::vector< PolyVals > cVec_; // vector of polynomial coefficients
  std::vector<int> polyTypeVec_; // vector of polynomial types

  ReactionRates( const Expr::Tag& tTag,
                 const Expr::Tag& pTag,
                 const Expr::TagList& massFracTags,
                 const Expr::Tag& mmwTag );
public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a ReactionRates expression
     *  @param resultTags tags for the net rate of production of each species
     *  @param tTag temperature
     *  @param pTag pressure
     *  @param massFracTags tag for the mass fraction of each species
     *  @param mmwTag tag for mixture molecular weight
     */
    Builder( const Expr::TagList& resultTags,
             const Expr::Tag& tTag,
             const Expr::Tag& pTag,
             const Expr::TagList& massFracTags,
             const Expr::Tag& mmwTag );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag tTag_;
    const Expr::Tag pTag_;
    const Expr::TagList massFracTags_;
    const Expr::Tag mmwTag_;
  };

  ~ReactionRates(){}
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
ReactionRates<FieldT>::
ReactionRates( const Expr::Tag& tTag,
               const Expr::Tag& pTag,
               const Expr::TagList& massFracTags,
               const Expr::Tag& mmwTag )
  : Expr::Expression<FieldT>(),
    tTag_( tTag ),
    tPowerTags_( Temperature<FieldT>::temperature_powers_tags() ), // temperature powers are auto-generated
    pTag_( pTag ),
    massFracTags_( massFracTags ),
    mmwTag_( mmwTag )
{
  this->set_gpu_runnable( true );
  Cantera_CXX::IdealGasMix* const gasMix = CanteraObjects::get_gasmix();

  nSpec_ = gasMix->nSpecies();

  rxnDataVec_ = gasMix->getReactionData(); // contains kinetics data for each reaction
  nRxns_ = rxnDataVec_.size();
  molecularWeights_ = gasMix->molecularWeights();
  molecularWeightsInv_.resize(nSpec_);
  for( size_t n=0; n<nSpec_; ++n )
    molecularWeightsInv_[n] = 1 / molecularWeights_[n];
  for( size_t r=0; r<nRxns_; ++r){
    const Cantera::ReactionData& rxnData = rxnDataVec_[r]; // ReactionData for reaction r

    /* here we are reading Cantera's product and reactant stoichiometric coefficients and storing them as ints
     * This speeds up evaluating the product (multiplication) of concentrations when evaluating reaction rates of progress
     * only elementary reactions are supported
     */
    std::vector<int> pStoich;
    std::vector<int> rStoich;
    const std::vector<double>& canteraRStoich = rxnData.rstoich;
    const std::vector<double>& canteraPStoich  = rxnData.pstoich;

    // dOrder is the difference between forward and reverse order, this is used in evaluating reverse rate constant
    int dOrder=0;
    std::vector<double>::const_iterator iStoich;
    for( iStoich = canteraRStoich.begin(); iStoich!=canteraRStoich.end(); ++iStoich ){
      if( fabs(*iStoich-1) <1e-2 )
        rStoich.push_back(1);
      else if( fabs(*iStoich-2) < 1e-2 )
        rStoich.push_back(2);
      else if( fabs(*iStoich-3) < 1e-2 )
        rStoich.push_back(3);
      else{
        std::cout << "Error in ReactionRates" << std::endl
            <<" Non-integer reactant stoichiometric coefficient" << std::endl
            <<" Reaction # "<< r <<", stoichiometric coefficient = " << *iStoich << std::endl;
        throw std::runtime_error("Problems in ReactionRates expression.");
      }
      dOrder += *iStoich;
    }
    for( iStoich = canteraPStoich.begin(); iStoich!=canteraPStoich.end(); ++iStoich ){
      if( fabs(*iStoich-1) <1e-2 )
        pStoich.push_back(1);
      else if( fabs(*iStoich-2) < 1e-2 )
        pStoich.push_back(2);
      else if( fabs(*iStoich-3) < 1e-2 )
        pStoich.push_back(3);
      else{
        std::cout << "Error in ReactionRates" << std::endl
            <<" Non-integer product stoichiometric coefficient" << std::endl
            <<" Reaction # "<< r <<", stoichiometric coefficient = " << *iStoich << std::endl;
        throw std::runtime_error("Problems in ReactionRates expression.");
      }
      dOrder -= *iStoich;
    }

    rStoich_.push_back(rStoich); // store stoich coeffs for each reaction
    pStoich_.push_back(pStoich);
    dOrder_.push_back(dOrder); // store dOrder for each reaction
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
      c[2]  = c[2]/2; // perform division of coefficients here - minor optimization
      c[3]  = c[3]/6;
      c[4]  = c[4]/12;
      c[5]  = c[5]/20;
      c[9]  = c[9]/2;
      c[10] = c[10]/6;
      c[11] = c[11]/12;
      c[12] = c[12]/20;
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
advertise_dependents( Expr::ExprDeps& exprDeps )
{
  exprDeps.requires_expression( tTag_         );
  exprDeps.requires_expression( tPowerTags_   );
  exprDeps.requires_expression( pTag_         );
  exprDeps.requires_expression( massFracTags_ );
  exprDeps.requires_expression( mmwTag_       );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
ReactionRates<FieldT>::
bind_fields( const Expr::FieldManagerList& fml )
{
  const typename Expr::FieldMgrSelector<FieldT>::type& fm = fml.field_manager<FieldT>();
  t_ = &fm.field_ref( tTag_ );

  tPowers_.clear();
  BOOST_FOREACH( const Expr::Tag& tag, tPowerTags_ ){
    tPowers_.push_back( &fm.field_ref(tag) );
  }

  p_ = &fm.field_ref( pTag_ );
  mmw_ = &fm.field_ref( mmwTag_ );

  massFracs_.clear();
  BOOST_FOREACH( const Expr::Tag& tag, massFracTags_ ){
    massFracs_.push_back( &fm.field_ref(tag) );
  }
}

//--------------------------------------------------------------------

template< typename FieldT >
void
ReactionRates<FieldT>::
evaluate()
{
  using namespace SpatialOps;
  using namespace Cantera;

  SpecT& rRates = this->get_value_vec();
  const FieldT& t = *t_;
  const FieldT& p = *p_;

  SpatFldPtr<FieldT> concPtr    = SpatialFieldStore::get<FieldT>(t); // total mass concentration [kg/m^3]
  SpatFldPtr<FieldT> conc2Ptr   = SpatialFieldStore::get<FieldT>(t); // total mass concentration squared
  SpatFldPtr<FieldT> logConcPtr = SpatialFieldStore::get<FieldT>(t); // log of molar concentration

  SpatFldPtr<FieldT> kPtr       = SpatialFieldStore::get<FieldT>(t); // forward rate constant
  SpatFldPtr<FieldT> krPtr      = SpatialFieldStore::get<FieldT>(t); // reverse rate constant

  SpatFldPtr<FieldT> logTPtr    = SpatialFieldStore::get<FieldT>(t); // log(t)

  SpatFldPtr<FieldT> dgPtr      = SpatialFieldStore::get<FieldT>(t); // delta gibbs energy for a reaction
  std::vector< SpatFldPtr<FieldT> > gPtrvec; // gibbs energy for each species
  for( size_t n=0; n<nSpec_; ++n)
    gPtrvec.push_back( SpatialFieldStore::get<FieldT>(t) );

  FieldT& conc    = *concPtr;
  FieldT& conc2   = *conc2Ptr;
  FieldT& logConc = *logConcPtr;

  FieldT& k  = *kPtr;
  FieldT& kr = *krPtr;

  FieldT& logT   = *logTPtr;
  FieldT& dg = *dgPtr;

  const FieldT& tRecip = *tPowers_[4];
  logT    <<= log(t);
  conc    <<= tRecip * p / GasConstant; // molar concentration
  logConc <<= log(conc);
  conc    <<= conc * *mmw_; // mass concentration

  { // calculate the delta gibbs energy for each species for use in evaluating reversible rate constants
    const FieldT& t2 = *tPowers_[0]; // t^2
    const FieldT& t3 = *tPowers_[1]; // t^3
    const FieldT& t4 = *tPowers_[2]; // t^4
    const FieldT& t5 = *tPowers_[3]; // t^5

    SpatFldPtr<FieldT> tlogTPtr = SpatialFieldStore::get<FieldT>(t); // t * log(t)
    FieldT&            tlogT    = *tlogTPtr;
    tlogT <<= t * logT;

    for( size_t n=0; n<nSpec_; ++n ){
      const int polyType = polyTypeVec_[n];
      const std::vector<double>& c = cVec_[n];
      const double minT = minTVec_[n];
      const double maxT = maxTVec_[n];
      switch ( polyType ){
        case SIMPLE:
          *gPtrvec[n] <<=        c[1] + c[3] * ( t    - c[0]        ) // H
                         - t * ( c[2] + c[3] * ( logT - log(c[0]) ) ); // -TS
          break;
        case NASA2:{
          /* polynomials are applicable in two temperature ranges - high and low
           * Caution: polynomials are used even if temperature is out of range
           * Note: coefficients have been divided by integers during construction of expression
           */
          *gPtrvec[n] <<= cond( t <= c[0], c[ 6] + (c[1] - c[ 7]) * t - c[2] * t2 - c[ 3] * t3 - c[ 4] * t4 - c[ 5] * t5 - c[1] * tlogT ) // if low temp
                              (            c[13] + (c[8] - c[14]) * t - c[9] * t2 - c[10] * t3 - c[11] * t4 - c[12] * t5 - c[8] * tlogT ); // else if high temp
          break;
        }
      } // switch
    } // species loop
  }

  for( size_t n=0; n<nSpec_; ++n )
    *rRates[n] <<= 0.0; // set net rates of production to 0.0 before starting summation over reactions

  for( size_t r=0; r<nRxns_; ++r ){
    const Cantera::ReactionData& rxnData = rxnDataVec_[r]; // ReactionData for reaction r

    const std::vector<double>& arrhenius = rxnData.rateCoeffParameters;
    // evaluation of forward rate constant
    if( fabs( arrhenius[1] ) < 1e-6 && fabs( arrhenius[2] ) < 1e-6 ) // rate constant does not depend on temp
      k <<= arrhenius[0];
    else
      k <<= arrhenius[0] * exp( arrhenius[1] * logT - arrhenius[2] * tRecip );

    if( rxnData.reactionType == THREE_BODY_RXN || rxnData.reactionType == FALLOFF_RXN ){ // third body
      SpatFldPtr<FieldT> mPtr  = SpatialFieldStore::get<FieldT>(t); // third body enhancement factor
      FieldT& m = *mPtr;

      m <<= rxnData.default_3b_eff * conc / *mmw_; // evaluate default third body enhancement
      for( std::map<int, double>::const_iterator iEff = rxnData.thirdBodyEfficiencies.begin(); iEff!= rxnData.thirdBodyEfficiencies.end(); ++iEff)
        m <<= m + *massFracs_[iEff->first] * conc * ( iEff->second - rxnData.default_3b_eff ) * molecularWeightsInv_[iEff->first]; // correct for non-default species;

      const int fallType = rxnData.falloffType;
      if( fallType == 0 ) // no falloff
        k <<= k * m;
      else{
        SpatFldPtr<FieldT> prPtr = SpatialFieldStore::get<FieldT>(t); // reduced pressure for falloff
        FieldT&            pr    = *prPtr;

        const std::vector<double>& auxParam = rxnData.auxRateCoeffParameters;
        pr <<= auxParam[0] * exp( auxParam[1] * logT - auxParam[2] * tRecip ) * m / k;

        switch( fallType ){
          case SIMPLE_FALLOFF: //Lindemann
            k <<= k * pr / (1 + pr);
            break;
          // Troe
          case TROE3_FALLOFF:  // fall through
          case TROE4_FALLOFF:{
            SpatFldPtr<FieldT> logPrPtr    = SpatialFieldStore::get<FieldT>(t); // log10(Pr)
            SpatFldPtr<FieldT> fCentPtr    = SpatialFieldStore::get<FieldT>(t); // "Fcent" factor for falloff evaluation
            SpatFldPtr<FieldT> logFCentPtr = SpatialFieldStore::get<FieldT>(t); // log10(Fcent)
            SpatFldPtr<FieldT> f1Ptr       = SpatialFieldStore::get<FieldT>(t); // "f1" factor for falloff evaluation

            FieldT& logPr    = *logPrPtr;
            FieldT& fCent    = *fCentPtr;
            FieldT& logFCent = *logFCentPtr;
            FieldT& f1       = *f1Ptr;

            const std::vector<double>& troe = rxnData.falloffParameters;

            if( fallType == TROE3_FALLOFF ) // Troe3 uses 3 parameters and has two terms
              logFCent <<= log10( cond( fabs(troe[1]) > 1e-8, (1-troe[0]) * exp(-t/troe[1]) ) (0.0) // if troe[1] is 0, this term is set to 0
                                + cond( fabs(troe[2]) > 1e-8,    troe[0]  * exp(-t/troe[2]) ) (0.0) ); // if troe[2] is 0, this term is set to 0
            else // Troe4 uses 4 parameters and has three terms
              logFCent <<= log10( cond( fabs(troe[1]) > 1e-8, (1-troe[0]) * exp( -t      / troe[1] ) ) (0.0) // if troe[1] is 0, this term is set to 0
                                + cond( fabs(troe[2]) > 1e-8,    troe[0]  * exp( -t      / troe[2] ) ) (0.0) // if troe[2] is 0, this term is set to 0
                                +                                    1.0  * exp( -tRecip * troe[3] )        );

            logPr <<= log10(pr);
#           define CTROE -0.4 - 0.67 * logFCent // macros to make f1 evaluation easier to read
#           define NTROE 0.75 - 1.27 * logFCent
            f1 <<= ( logPr + CTROE ) / ( NTROE - 0.14 * (logPr +CTROE) ); // calculate the "f1" factor

            fCent <<= pow(10, logFCent / (1 + f1*f1) ); // calculate the "Fcent" factor
            k <<= fCent * k * pr / (1+pr);

            break;
          } // Troe
          default:{
            std::ostringstream msg;
            msg << "Error in " __FILE__ << " : " << __LINE__ << std::endl
                <<" Falloff type not supported, type = " << fallType << std::endl
                <<" See <cantera/kernel/reaction_defs.h> for type definition" << std::endl;
            throw std::runtime_error( msg.str() );
          }
        } // switch (fallType)
      } // falloff
    } // third body

    const std::vector<int>& rStoich = rStoich_[r]; // integer stoich. coefficients stored when expression was built
    const std::vector<int>& pStoich = pStoich_[r];

    const std::vector<int>& reactants = rxnData.reactants; // species index for each reactant or product
    const std::vector<int>& products  = rxnData.products;

    std::vector<int>::const_iterator iStoich; // stoichiometric coefficient
    std::vector<int>::const_iterator iSpec; // species index

    if( rxnData.reversible == true ){
      dg <<= 0.0; // set value to 0 before summation

      iStoich = rStoich.begin();
      for( iSpec = reactants.begin(); iSpec!=reactants.end(); ++iSpec, ++iStoich)
        dg <<= dg - *iStoich * *gPtrvec[*iSpec];
      iStoich = pStoich.begin();
      for( iSpec = products.begin(); iSpec!=products.end(); ++iSpec, ++iStoich)
        dg <<= dg + *iStoich * *gPtrvec[*iSpec];

      kr <<= k * exp( dg * tRecip / GasConstant + dOrder_[r] * logConc ); // reversible rate constant
    }
    else
      kr <<= 0.0; // 0 if not reversible

    iStoich = rStoich.begin();
      for( iSpec=reactants.begin(); iSpec!=reactants.end(); ++iSpec, ++iStoich){ // calculate r_f = k*C_i^stoich
        switch ( *iStoich ){
          case 1:
            k <<= k * *massFracs_[*iSpec] * conc * molecularWeightsInv_[*iSpec];
            break;
          case 2:
            k <<= k * *massFracs_[*iSpec] * *massFracs_[*iSpec] * conc * molecularWeightsInv_[*iSpec] * conc * molecularWeightsInv_[*iSpec];
            break;
          case 3:
            k <<= k * *massFracs_[*iSpec] * conc * molecularWeightsInv_[*iSpec] * *massFracs_[*iSpec] * conc * molecularWeightsInv_[*iSpec] * *massFracs_[*iSpec] * conc * molecularWeightsInv_[*iSpec];
            break;
          default:{
            std::ostringstream msg;
            msg << "Error in " << __FILE__ << " : " << __LINE__ << std::endl
                <<" Non-integer reactant stoichiometric coefficient" << std::endl
                <<" Reaction # "<< r <<", stoichiometric coefficient = " << *iStoich << std::endl;
            throw std::runtime_error( msg.str() );
          }
        }
      }

      if( rxnData.reversible == true){
        iStoich = pStoich.begin();
          for( iSpec=products.begin(); iSpec!=products.end(); ++iSpec, ++iStoich ){ // calculate r_r = kr*C_i^stoich
            switch ( *iStoich ){
              case 1:
                kr <<= kr * *massFracs_[*iSpec] * conc * molecularWeightsInv_[*iSpec];
                break;
              case 2:
                kr <<= kr * *massFracs_[*iSpec] * *massFracs_[*iSpec] * conc * molecularWeightsInv_[*iSpec] * conc * molecularWeightsInv_[*iSpec];
                break;
              case 3:
                kr <<= kr * *massFracs_[*iSpec] * conc * molecularWeightsInv_[*iSpec] * *massFracs_[*iSpec] * conc * molecularWeightsInv_[*iSpec] * *massFracs_[*iSpec] * conc * molecularWeightsInv_[*iSpec];
                break;
              default:{
                std::ostringstream msg;
                msg << "Error in " << __FILE__ << " : " << __LINE__ << std::endl
                    <<" Non-integer product stoichiometric coefficient" << std::endl
                    <<" Reaction # "<< r <<", stoichiometric coefficient = " << *iStoich << std::endl;
                throw std::runtime_error( msg.str() );
              }
            }
          }
      }


    iStoich = rStoich.begin();
    for( iSpec = reactants.begin(); iSpec!=reactants.end(); ++iSpec, ++iStoich )
      *rRates[*iSpec] <<= *rRates[*iSpec] - *iStoich * (k - kr) * molecularWeights_[*iSpec]; // sum rate of progress into species rates of production

    iStoich = pStoich.begin();
    for( iSpec = products.begin(); iSpec!=products.end(); ++iSpec, ++iStoich )
      *rRates[*iSpec] <<= *rRates[*iSpec] + *iStoich * (k - kr) * molecularWeights_[*iSpec];

  }
}

//--------------------------------------------------------------------

template< typename FieldT >
ReactionRates<FieldT>::
Builder::Builder( const Expr::TagList& resultTags,
                  const Expr::Tag& tTag,
                  const Expr::Tag& pTag,
                  const Expr::TagList& massFracTags,
                  const Expr::Tag& mmwTag )
: ExpressionBuilder( resultTags ),
  tTag_( tTag ),
  pTag_( pTag ),
  massFracTags_( massFracTags ),
  mmwTag_( mmwTag )
{}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
ReactionRates<FieldT>::
Builder::build() const
{
  return new ReactionRates<FieldT>( tTag_, pTag_, massFracTags_, mmwTag_ );
}

} // namespace pokitt

#endif // ReactionRates_Expr_h
