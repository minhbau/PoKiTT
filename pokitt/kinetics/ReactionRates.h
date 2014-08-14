#ifndef ReactionRates_Expr_h
#define ReactionRates_Expr_h

//#define TIMINGS

#include <expression/Expression.h>
#ifdef TIMINGS
#include <boost/timer.hpp>
#endif

#include <pokitt/CanteraObjects.h> // include cantera wrapper

#include <cantera/kernel/ct_defs.h> // contains value of gas constant
#include <cantera/kernel/reaction_defs.h>
#include <cantera/kernel/speciesThermoTypes.h> // contains definitions for which polynomial is being used
#include <cantera/IdealGasMix.h>

namespace Cantera_CXX{ class IdealGasMix; } // location of polynomial coefficients for gibbs energy

/**
 *  \class  ReactionRates
 *  \author Nate Yonkee
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
  const Expr::Tag pTag_;
  Expr::TagList massFracTags_;
  const Expr::Tag mmwTag_;
  const FieldT* t_;
  const FieldT* p_;
  const FieldT* mmw_; // mixture molecular weight
  std::vector<const FieldT*> massFracs_;
  Cantera_CXX::IdealGasMix* const gasMix_;
  const double pressure_;
  const int nSpec_; // number of species in the mechanism
  const int nRxns_; // number of reactions in the mechanism
  std::vector<int> dOrder_; // difference between the forward reaction order and backwards reaction order
  std::vector<std::vector<int> > rstoich_; // reactant stoich coeffs, stores ints instead of the doubles from Cantera
  std::vector<std::vector<int> > pstoich_; // product stoich coeffs, stores ints instead of the doubles from Cantera


  ReactionRates( const Expr::Tag& tTag,
                 const Expr::Tag& pTag,
                 const Expr::Tag& massFracTag,
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
     *  @param massFracTag tag for the mass fraction of each species
     *  @param mmwTag tag for mixture molecular weight
     */
    Builder( const Expr::TagList& resultTags,
             const Expr::Tag& tTag,
             const Expr::Tag& pTag,
             const Expr::Tag& massFracTag,
             const Expr::Tag& mmwTag );

    Expr::ExpressionBase* build() const;

  private:
    const Expr::Tag tTag_;
    const Expr::Tag pTag_;
    const Expr::Tag massFracTag_;
    const Expr::Tag mmwTag_;
  };

  ~ReactionRates();
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
    const Expr::Tag& massFracTag,
    const Expr::Tag& mmwTag )
    : Expr::Expression<FieldT>(),
      tTag_( tTag ),
      pTag_( pTag ),
      mmwTag_( mmwTag ),
      gasMix_( CanteraObjects::get_gasmix() ),
      pressure_( gasMix_->pressure()),
      nSpec_( gasMix_->nSpecies() ),
      nRxns_( gasMix_->getReactionData().size())
      {
  this->set_gpu_runnable( true );

  massFracTags_.clear();
  for( size_t n=0; n<nSpec_; ++n ){
    std::ostringstream name;
    name << massFracTag.name() << "_" << n;
    massFracTags_.push_back( Expr::Tag(name.str(),massFracTag.context()) );
  }

  const std::vector<Cantera::ReactionData>& rxnDataVec = gasMix_->getReactionData(); // contains kinetics data for each reaction
  for( size_t r=0; r<nRxns_; ++r){
    const Cantera::ReactionData& rxnData = rxnDataVec[r]; // ReactionData for reaction r

    /* here we are reading Cantera's product and reactant stoichiometric coefficients and storing them as ints
     * This speeds up evaluating the product (multiplication) of concentrations when evaluating reaction rates of progress
     */
    std::vector<int> pstoich;
    std::vector<int> rstoich;
    const std::vector<double>& reactantStoich = rxnData.rstoich;
    const std::vector<double>& productStoich = rxnData.pstoich;

    // dOrder is the difference between forward and reverse order, this is used in evaluating reverse rate constant
    double dOrder=0.0;

    std::vector<double>::const_iterator sit;
    for( sit = reactantStoich.begin(); sit!=reactantStoich.end(); ++sit){
      if( fabs(*sit-1) <1e-4 )
        rstoich.push_back(1);
      else if( fabs(*sit-2) < 1e-4 )
        rstoich.push_back(2);
      else if( fabs(*sit-3) < 1e-4 )
        rstoich.push_back(3);
      else{
        std::cout << "Error in ReactionRates" << std::endl
            <<" Non-integer reactant stoichiometric coefficient" << std::endl
            <<" Reaction # "<< r <<", stoichiometric coefficient = " << *sit << std::endl;
        throw std::runtime_error("Problems in ReactionRates expression.");
      }
      dOrder+=*sit;
    }

    for( sit = productStoich.begin(); sit!=productStoich.end(); ++sit){
      if( fabs(*sit-1) <1e-4 )
        pstoich.push_back(1);
      else if( fabs(*sit-2) < 1e-4 )
        pstoich.push_back(2);
      else if( fabs(*sit-3) < 1e-4 )
        pstoich.push_back(3);
      else{
        std::cout << "Error in ReactionRates" << std::endl
            <<" Non-integer product stoichiometric coefficient" << std::endl
            <<" Reaction # "<< r <<", stoichiometric coefficient = " << *sit << std::endl;
        throw std::runtime_error("Problems in ReactionRates expression.");
      }
      dOrder-=*sit;
    }

    rstoich_.push_back(rstoich); // store stoich coeffs for each reaction
    pstoich_.push_back(pstoich);
    dOrder_.push_back(dOrder); // store dOrder for each reaction
  }
      }



//--------------------------------------------------------------------

template< typename FieldT >
ReactionRates<FieldT>::
~ReactionRates()
{
  CanteraObjects::restore_gasmix(gasMix_);
}

//--------------------------------------------------------------------

template< typename FieldT >
void
ReactionRates<FieldT>::
advertise_dependents( Expr::ExprDeps& exprDeps )
{
  exprDeps.requires_expression( tTag_ );
  exprDeps.requires_expression( pTag_ );
  exprDeps.requires_expression( massFracTags_ );
  exprDeps.requires_expression( mmwTag_ );
}

//--------------------------------------------------------------------

template< typename FieldT >
void
ReactionRates<FieldT>::
bind_fields( const Expr::FieldManagerList& fml )
{
  const typename Expr::FieldMgrSelector<FieldT>::type& fm = fml.field_manager<FieldT>();
  t_ = &fm.field_ref( tTag_ );
  p_ = &fm.field_ref( pTag_ );
  mmw_ = &fm.field_ref( mmwTag_ );
  for (size_t n=0; n<nSpec_; ++n) massFracs_.push_back(&fm.field_ref( massFracTags_[n] ));
}

//--------------------------------------------------------------------

template< typename FieldT >
void
ReactionRates<FieldT>::
evaluate()
{
#ifdef TIMINGS
  boost::timer timer;
#endif
  using namespace SpatialOps;
  using namespace Cantera;

  SpecT& rRates = this->get_value_vec();
  const FieldT& t = *t_;
  const FieldT& p = *p_;

  int mtype;
# ifdef ENABLE_CUDA
  mtype = GPU_INDEX;
# else
  mtype = CPU_INDEX;
# endif

  SpatFldPtr<FieldT> concPtr  = SpatialFieldStore::get<FieldT>(t); // total mass concentration [kg/m^3]
  SpatFldPtr<FieldT> conc2Ptr = SpatialFieldStore::get<FieldT>(t); // total mass concentration squared
  SpatFldPtr<FieldT> logConcPtr  = SpatialFieldStore::get<FieldT>(t); // log of molar concentration

  SpatFldPtr<FieldT> kPtr  = SpatialFieldStore::get<FieldT>(t,mtype); // forward rate constant
  SpatFldPtr<FieldT> krPtr  = SpatialFieldStore::get<FieldT>(t); // reverse rate constant

  SpatFldPtr<FieldT> trecipPtr  = SpatialFieldStore::get<FieldT>(t); // t^-1
  SpatFldPtr<FieldT> logtPtr  = SpatialFieldStore::get<FieldT>(t); // log(t)

  SpatFldPtr<FieldT> dgPtr = SpatialFieldStore::get<FieldT>(t); // delta gibbs energy for a reaction
  std::vector<SpatFldPtr<FieldT> > gPtrvec; // gibbs energy for each species
  for( size_t n=0; n<nSpec_; ++n)
    gPtrvec.push_back( SpatialFieldStore::get<FieldT>(t) );

  FieldT& conc = *concPtr;
  FieldT& conc2 = *conc2Ptr;
  FieldT& logConc = *logConcPtr;

  FieldT& k = *kPtr;
  FieldT& kr = *krPtr;

  FieldT& trecip = *trecipPtr;
  FieldT& logt = *logtPtr;

  FieldT& dg = *dgPtr;

  trecip <<= 1/ t;
  logt <<= log(t);
  conc <<= trecip * p / GasConstant; // molar concentration
  logConc <<= log(conc);
  conc <<= conc * *mmw_; // mass concentration
  conc2 <<= conc * conc;

  { // add scope resolution so that fields of powers of t go out of scope after evaluation of gibbs energy polynomial
    // pre-compute powers of t used in gibbs energy polynomial for each species
    SpatFldPtr<FieldT> t2Ptr  = SpatialFieldStore::get<FieldT>(t); // t^2
    SpatFldPtr<FieldT> t3Ptr  = SpatialFieldStore::get<FieldT>(t); // t^3
    SpatFldPtr<FieldT> t4Ptr  = SpatialFieldStore::get<FieldT>(t); // t^4
    SpatFldPtr<FieldT> t5Ptr  = SpatialFieldStore::get<FieldT>(t); // t^5
    SpatFldPtr<FieldT> tlogtPtr  = SpatialFieldStore::get<FieldT>(t); // t * log(t)

    FieldT& t2 = *t2Ptr;
    FieldT& t3 = *t3Ptr;
    FieldT& t4 = *t4Ptr;
    FieldT& t5 = *t5Ptr;
    FieldT& tlogt = *tlogtPtr;

    t2 <<= t * t;
    t3 <<= t2 * t;
    t4 <<= t3 * t;
    t5 <<= t4 * t;
    tlogt <<= t * logt;

    const Cantera::SpeciesThermo& spThermo = gasMix_->speciesThermo();
    std::vector<double> c(15,0); // vector of Cantera's polynomial coefficients
    int polyType; // type of polynomial
    double minT; // minimum temperature where polynomial is valid
    double maxT; // maximum temperature where polynomial is valid
    double refPressure;

    for( size_t n=0; n<nSpec_; ++n ){
      spThermo.reportParams(n, polyType, &c[0], minT, maxT, refPressure);
      if ( polyType==SIMPLE){
        *gPtrvec[n] <<=        c[1] + c[3] * ( t - c[0]           ) // H
                       - t * ( c[2] + c[3] * ( logt - log(c[0]) ) ); // -TS
      }
      /* polynomials are applicable in two temperature ranges - high and low
        * If the temperature is out of range, the value is set to the value at the min or max temp
        */
      else if ( polyType==NASA2 ){
        std::vector<double>::iterator ic = c.begin() + 1;
        std::vector<double>::iterator icend = c.end();
        for( ; ic != icend; ++ic)
          *ic *= GasConstant; // dimensionalize the coefficients
        *gPtrvec[n] <<= cond( t <= c[0], c[ 6] + (c[1] - c[ 7]) * t - c[2]/2 * t2 - c[ 3]/6 * t3 - c[ 4]/12 * t4 - c[ 5]/20 * t5 - c[1] * tlogt ) // if low temp
                            (            c[13] + (c[8] - c[14]) * t - c[9]/2 * t2 - c[10]/6 * t3 - c[11]/12 * t4 - c[12]/20 * t5 - c[8] * tlogt ); // else if high temp
      }
      else{
        std::cout << "Error in ReactionRates" << std::endl
             <<" Thermo type not supported, type = " << polyType << std::endl
             <<" See <cantera/kernel/SpeciesThermoInterpType.h> for type definition" << std::endl;
        throw std::runtime_error("Problems in ReactionRates expression.");
      }
    }
  } // powers of t now go out of scope and free the associated memory

  for( size_t n=0; n<nSpec_; ++n )
    *rRates[n] <<= 0.0; // set net rates of production to 0.0 before starting summation over reactions

  const std::vector<Cantera::ReactionData>& rxnDataVec = gasMix_->getReactionData(); // vector of ReactionData for every reaction

  const std::vector<double>& molecularWeights = gasMix_->molecularWeights();

  std::vector<double> molecularWeightsInv(nSpec_);
  for( size_t n=0; n<nSpec_; ++n)
    molecularWeightsInv[n] = 1 / molecularWeights[n];

  for( size_t r=0; r<nRxns_; ++r ){

    const Cantera::ReactionData& rxnData = rxnDataVec[r]; // ReactionData for reaction r

    const std::vector<double>& arrhenius = rxnData.rateCoeffParameters;

    if( fabs( arrhenius[1] ) < 1e-6 && fabs( arrhenius[2] ) < 1e-6 ) // evaluation of forward rate constant
      k <<= arrhenius[0];
    else
      k <<= arrhenius[0] * exp( arrhenius[1] * logt - arrhenius[2] * trecip );

    if( rxnData.reactionType == THREE_BODY_RXN || rxnData.reactionType == FALLOFF_RXN ){ // third body
      SpatFldPtr<FieldT> mPtr  = SpatialFieldStore::get<FieldT>(t); // third body enhancement factor
      FieldT& m = *mPtr;

      m <<= rxnData.default_3b_eff * conc / *mmw_; // evaluate default third body enhancement
      for( std::map<int, double>::const_iterator miter = rxnData.thirdBodyEfficiencies.begin(); miter!= rxnData.thirdBodyEfficiencies.end(); ++miter)
        m <<= m + *massFracs_[miter->first] * conc * ( miter->second - rxnData.default_3b_eff ) * molecularWeightsInv[miter->first]; // correct for non-default species;

      int fallType = rxnData.falloffType;
      if( fallType == 0 ) // no falloff
        k <<= k * m;
      else{
        SpatFldPtr<FieldT> PrPtr  = SpatialFieldStore::get<FieldT>(t); // reduced pressure for falloff
        SpatFldPtr<FieldT> logPrPtr  = SpatialFieldStore::get<FieldT>(t); // log10(Pr)
        SpatFldPtr<FieldT> FcentPtr  = SpatialFieldStore::get<FieldT>(t); // "Fcent" factor for falloff evaluation
        SpatFldPtr<FieldT> logFcentPtr  = SpatialFieldStore::get<FieldT>(t); // log10(Fcent)
        SpatFldPtr<FieldT> f1Ptr  = SpatialFieldStore::get<FieldT>(t); // "f1" factor for falloff evaluation

        FieldT& Pr = *PrPtr;
        FieldT& logPr = *logPrPtr;
        FieldT& Fcent = *FcentPtr;
        FieldT& logFcent = *logFcentPtr;
        FieldT& f1 = *f1Ptr;

        const std::vector<double>& auxparam = rxnData.auxRateCoeffParameters;

        Pr <<= auxparam[0] * exp( auxparam[1] * logt - auxparam[2] * trecip ) * m / k;

        if( fallType == SIMPLE_FALLOFF ) //Lindemann
          k <<= k * Pr / (1 + Pr);

        else if ( fallType == TROE3_FALLOFF || fallType == TROE4_FALLOFF ){ // Troe
          const std::vector<double>& troe = rxnData.falloffParameters;

          if( fallType == TROE3_FALLOFF ) // Troe3 uses 3 parameters and has two terms
            logFcent <<= log10( cond( fabs(troe[1]) > 1e-8, (1-troe[0]) * exp(-t/troe[1]) ) (0.0) // if troe[1] is 0, this term is set to 0
                              + cond( fabs(troe[2]) > 1e-8,    troe[0]  * exp(-t/troe[2]) ) (0.0) ); // if troe[2] is 0, this term is set to 0
          else // Troe4 uses 4 parameters and has three terms
            logFcent <<= log10( cond( fabs(troe[1]) > 1e-8, (1-troe[0]) * exp( -t      / troe[1] ) ) (0.0) // if troe[1] is 0, this term is set to 0
                              + cond( fabs(troe[2]) > 1e-8,    troe[0]  * exp( -t      / troe[2] ) ) (0.0) // if troe[2] is 0, this term is set to 0
                              +                                    1.0  * exp( -trecip * troe[3] )        );

          logPr <<= log10(Pr);
#define CTROE -0.4 - 0.67 * logFcent // macros to make f1 evaluation easier to read
#define NTROE 0.75 - 1.27 * logFcent
          f1 <<= (logPr+CTROE) / ( NTROE - 0.14 * (logPr+CTROE) ); // calculate the "f1" factor

          Fcent <<= pow(10, logFcent / (1 + f1*f1) ); // calculate the "Fcent" factor
          k <<= Fcent * k * Pr / (1+Pr);
        } // Troe
        else{
          std::cout << "Error in ReactionRates" << std::endl
               <<" Falloff type not supported, type = " << fallType << std::endl
               <<" See <cantera/kernel/reaction_defs.h> for type definition" << std::endl;
          throw std::runtime_error("Problems in ReactionRates expression.");
        }
      } // falloff
    } // third body

    const std::vector<int>& rstoich = rstoich_[r]; // integer stoich. coefficients stored when expression was built
    const std::vector<int>& pstoich = pstoich_[r];

    const std::vector<int>& reactants = rxnData.reactants; // species index for each reactant or product
    const std::vector<int>& products = rxnData.products;

    std::vector<int>::const_iterator sit;
    std::vector<int>::const_iterator rit;

    if( rxnData.reversible == true){
      dg <<= 0.0; // set value to 0 before summation

      sit = rstoich.begin();
      for( rit = reactants.begin(); rit!=reactants.end(); ++rit, ++sit)
        dg <<= dg - *sit * *gPtrvec[*rit];
      sit = pstoich.begin();
      for( rit = products.begin(); rit!=products.end(); ++rit, ++sit)
        dg <<= dg + *sit * *gPtrvec[*rit];

      kr <<= k * exp( dg * trecip / GasConstant + dOrder_[r] * logConc ); // reversible rate constant
    }
    else
      kr <<= 0.0; // 0 if not reversible

    sit = rstoich.begin();
      for( rit=reactants.begin(); rit!=reactants.end(); ++rit, ++sit){
        switch (*sit){
        case 1:
//          k <<= k * *massFracs_[*rit] * conc / molecularWeights[*rit];
          k <<= k * *massFracs_[*rit] * conc * molecularWeightsInv[*rit];
          break;
        case 2:
//          k <<= k * *massFracs_[*rit] * *massFracs_[*rit] * conc / molecularWeights[*rit] * conc / molecularWeights[*rit];
          k <<= k * *massFracs_[*rit] * *massFracs_[*rit] * conc * molecularWeightsInv[*rit] * conc * molecularWeightsInv[*rit];
//          k <<= k * *massFracs_[*rit] * *massFracs_[*rit] * conc2 *  molecularWeightsInv2[*rit];
          break;
        case 3:
          k <<= k * *massFracs_[*rit] * conc * molecularWeightsInv[*rit] * *massFracs_[*rit] * conc * molecularWeightsInv[*rit] * *massFracs_[*rit] * conc * molecularWeightsInv[*rit];
          break;
        default:
          std::cout << "Error in ReactionRates" << std::endl
              <<" Non-integer reactant stoichiometric coefficient" << std::endl
              <<" Reaction # "<< r <<", stoichiometric coefficient = " << *sit << std::endl;
          throw std::runtime_error("Problems in ReactionRates expression.");
          break;
        }
      }

      if( rxnData.reversible == true){
        sit = pstoich.begin();
          for( rit=products.begin(); rit!=products.end(); ++rit, ++sit){
            switch (*sit){
            case 1:
//              kr <<= kr * *massFracs_[*rit] * conc / molecularWeights[*rit];
              kr <<= kr * *massFracs_[*rit] * conc * molecularWeightsInv[*rit];
              break;
            case 2:
//              kr <<= kr * *massFracs_[*rit] * *massFracs_[*rit] * conc / molecularWeights[*rit] * conc / molecularWeights[*rit];
              kr <<= kr * *massFracs_[*rit] * *massFracs_[*rit] * conc * molecularWeightsInv[*rit] * conc * molecularWeightsInv[*rit];
//              kr <<= kr * *massFracs_[*rit] * *massFracs_[*rit] * conc2 *  molecularWeightsInv2[*rit];
              break;
            case 3:
              kr <<= kr * *massFracs_[*rit] * conc * molecularWeightsInv[*rit] * *massFracs_[*rit] * conc * molecularWeightsInv[*rit] * *massFracs_[*rit] * conc * molecularWeightsInv[*rit];
              break;
            default:
              std::cout << "Error in ReactionRates" << std::endl
                  <<" Non-integer product stoichiometric coefficient" << std::endl
                  <<" Reaction # "<< r <<", stoichiometric coefficient = " << *sit << std::endl;
              throw std::runtime_error("Problems in ReactionRates expression.");
              break;
            }
          }
      }


    sit = rstoich.begin();
    for( rit = reactants.begin(); rit!=reactants.end(); ++rit, ++sit)
      *rRates[*rit] <<= *rRates[*rit] - *sit * (k - kr) * molecularWeights[*rit]; // sum rate of progress into species rates of production

    sit = pstoich.begin();
    for( rit = products.begin(); rit!=products.end(); ++rit, ++sit)
      *rRates[*rit] <<= *rRates[*rit] + *sit * (k - kr) * molecularWeights[*rit];

  }
#ifdef TIMINGS
    std::cout<<"rr time "<<timer.elapsed()<<std::endl;
#endif
}

//--------------------------------------------------------------------

template< typename FieldT >
ReactionRates<FieldT>::
Builder::Builder( const Expr::TagList& resultTags,
                  const Expr::Tag& tTag,
                  const Expr::Tag& pTag,
                  const Expr::Tag& massFracTag,
                  const Expr::Tag& mmwTag )
: ExpressionBuilder( resultTags ),
  tTag_( tTag ),
  pTag_( pTag ),
  massFracTag_( massFracTag ),
  mmwTag_( mmwTag )
{}

//--------------------------------------------------------------------

template< typename FieldT >
Expr::ExpressionBase*
ReactionRates<FieldT>::
Builder::build() const
{
  return new ReactionRates<FieldT>( tTag_, pTag_, massFracTag_, mmwTag_ );
}


#endif // ReactionRates_Expr_h
