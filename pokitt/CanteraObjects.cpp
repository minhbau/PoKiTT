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

#include <stdexcept>
#include <cassert>
#include <iostream>
#include <memory>

#include "CanteraObjects.h"

#include <cantera/thermo.h>
#include <cantera/base/Solution.h>

#include <cantera/base/ct_defs.h>              // value of gas constant
#include <cantera/thermo/speciesThermoTypes.h> // definitions for which polynomial is being used
#include <cantera/kinetics/GasKinetics.h>
#include <cantera/kinetics/Reaction.h>
#include <cantera/transport/MixTransport.h>

#ifdef EXPRESSION_THREADS

#include <boost/thread/mutex.hpp>
class CanteraMutex
{
  const boost::mutex::scoped_lock lock;
  inline boost::mutex& get_mutex() const
  {
    static boost::mutex m; return m;
  }
public:
  CanteraMutex() : lock( get_mutex() ) {}
  ~CanteraMutex(){}
};

#else

struct CanteraMutex{
  inline explicit CanteraMutex(){}
  inline ~CanteraMutex(){}
};

#endif  // EXPRESSION_THREADS

ThermData::ThermData( const Cantera::MultiSpeciesThermo& spThermo, const int i )
  : index( i )
{
  double refPressure; int cType;
  coefficients.resize( 15 );
  spThermo.reportParams(i, cType, &(coefficients[0]), minTemp, maxTemp, refPressure);
  switch( cType ){
    case SIMPLE     :
    case CONSTANT_CP: type = CONST_POLY;   break;
    case NASA2      : type = NASA_POLY;    break;
    default:{
      std::ostringstream msg;
      msg << __FILE__ << " : " << __LINE__
          << "\nThermo type not currently supported,\n Type = " << cType
          << ", species # " << i
          << "\nSee speciesThermoTypes.h in Cantera\n\n";
      throw std::runtime_error( msg.str() );
    }
  }
}

RxnData::SpeciesRxnData::
SpeciesRxnData( int index,
                int stoich,
                double mw,
                double thdBdyEff )
  : index( index ),
    stoich( stoich ),
    mw( mw ),
    invMW( 1/mw ),
    thdBdyEff( thdBdyEff )
{}

RxnData::RxnData( const IdealGasPtr gas,
                  const Cantera::Reaction& rxn,
                  const std::vector<double>& MW )
  : reversible( rxn.reversible )
{
  const std::string rxnType = rxn.type();

  for( int i=0; i<3; ++i ){
    kFwdCoefs     [i] = 0.0;
    kPressureCoefs[i] = 0.0;
  }
  for( int i=0; i<4; ++i ) troeParams[i]=0;

  // ensure that the type of each reaction is supported by PoKiTT.
  if( rxnType.compare( "elementary"        ) != 0 &&
      rxnType.compare( "elementary-legacy" ) != 0 &&
      rxnType.compare( "three-body"        ) != 0 &&
      rxnType.compare( "three-body-legacy" ) != 0 &&
      rxnType.compare( "falloff"           ) != 0  )
  {
    std::ostringstream msg;
    msg << __FILE__ << " : " << __LINE__
        << "\nUnsupported reaction type encountered for \n"
        << "\t" << rxn.equation() << "\n\t(Rxn type: " << rxnType << ")\n";
    throw std::runtime_error( msg.str() );
  }

  if( rxnType.compare( "elementary" ) == 0 ){
    const auto rate = std::dynamic_pointer_cast<Cantera::ArrheniusRate>( const_cast<Cantera::Reaction&>(rxn).rate() );
    kFwdCoefs[0] = rate->preExponentialFactor();
    kFwdCoefs[1] = rate->temperatureExponent();
    kFwdCoefs[2] = rate->activationEnergy_R();
    type = ELEMENTARY;
  }
  else if( rxnType.compare( "elementary-legacy" ) == 0 ){
    const Cantera::ElementaryReaction2& elemRate = dynamic_cast<Cantera::ElementaryReaction2&>( const_cast<Cantera::Reaction&>(rxn) );
    kFwdCoefs[0] = elemRate.rate.preExponentialFactor();
    kFwdCoefs[1] = elemRate.rate.temperatureExponent();
    kFwdCoefs[2] = elemRate.rate.activationEnergy_R();
    type = ELEMENTARY;
  }
  else if( rxnType.compare( "falloff" ) == 0 ){
    const Cantera::FalloffReaction& fallRxn = dynamic_cast<const Cantera::FalloffReaction&>(rxn);
    kFwdCoefs     [0] = fallRxn.high_rate.preExponentialFactor();
    kFwdCoefs     [1] = fallRxn.high_rate.temperatureExponent();
    kFwdCoefs     [2] = fallRxn.high_rate.activationEnergy_R();
    kPressureCoefs[0] = fallRxn.low_rate.preExponentialFactor();
    kPressureCoefs[1] = fallRxn.low_rate.temperatureExponent();
    kPressureCoefs[2] = fallRxn.low_rate.activationEnergy_R();
    fallRxn.falloff->getParameters( troeParams );
    thdBdyDefault = fallRxn.third_body.default_efficiency;
    for( const Cantera::Composition::value_type& vt : fallRxn.third_body.efficiencies ){
      const int spIx = gas->speciesIndex(vt.first);
      thdBdySpecies.push_back( SpeciesRxnData( spIx, 0, MW[spIx], vt.second ) );
    }

    const std::string rxnType = fallRxn.falloff->type();
    if( rxnType.compare("Lindemann") == 0 )
      type = LINDEMANN;
    else if( rxnType.compare("Troe") == 0 or rxnType.compare("SRI") == 0 )
      type = TROE;
    else{
      std::ostringstream msg;
      msg << __FILE__ << " : " << __LINE__
          << "\nFalloff type not supported,\n Type = " << rxnType
          << "\nSee reaction_defs.h in Cantera\n";
      throw std::runtime_error( msg.str() );
    }
  }
  else if( rxnType.compare( "three-body" ) == 0 ){
    Cantera::ThreeBodyReaction3& tbRxn = dynamic_cast<Cantera::ThreeBodyReaction3&>( const_cast<Cantera::Reaction&>(rxn) );
    const auto rate = std::dynamic_pointer_cast<Cantera::ArrheniusRate>( tbRxn.rate() );
    kFwdCoefs[0]  = rate->preExponentialFactor();
    kFwdCoefs[1]  = rate->temperatureExponent();
    kFwdCoefs[2]  = rate->activationEnergy_R();
    thdBdyDefault = tbRxn.thirdBody()->default_efficiency;
    type          = THIRD_BODY;
    for( const Cantera::Composition::value_type& vt : tbRxn.thirdBody()->efficiencies ){
      const int spIx = gas->speciesIndex(vt.first);
      thdBdySpecies.push_back( SpeciesRxnData( spIx, 0, MW[spIx], vt.second ) );
    }
  }
  else if( rxnType.compare( "three-body-legacy" ) == 0 ){
    Cantera::ThreeBodyReaction2& tbRxn = dynamic_cast<Cantera::ThreeBodyReaction2&>( const_cast<Cantera::Reaction&>(rxn) );
    kFwdCoefs[0]  = tbRxn.rate.preExponentialFactor();
    kFwdCoefs[1]  = tbRxn.rate.temperatureExponent();
    kFwdCoefs[2]  = tbRxn.rate.activationEnergy_R();
    thdBdyDefault = tbRxn.third_body.default_efficiency;
    type          = THIRD_BODY;
    for( const Cantera::Composition::value_type& vt : tbRxn.third_body.efficiencies ){
      const int spIx = gas->speciesIndex(vt.first);
      thdBdySpecies.push_back( SpeciesRxnData( spIx, 0, MW[spIx], vt.second ) );
    }
  }

  /* Here we are reading Cantera's product and reactant stoichiometric coefficients and storing them as ints
   * This helps performance because we only support elementary reaction mechanisms
   * We also check if stoichiometry is consistent with elementary reactions
   */
  for( const Cantera::Composition::value_type& reactant : rxn.reactants ){
    const std::string& spNam = reactant.first;
    const size_t spIx = gas->speciesIndex(spNam);
    const double stoich = reactant.second;
    const int istoich = int(stoich);
    if(      fabs( stoich - 1 ) < 1e-10 ) reactants.push_back( SpeciesRxnData( spIx, 1, MW[spIx], thdBdyDefault ) );
    else if( fabs( stoich - 2 ) < 1e-10 ) reactants.push_back( SpeciesRxnData( spIx, 2, MW[spIx], thdBdyDefault ) );
    else if( fabs( stoich - 3 ) < 1e-10 ) reactants.push_back( SpeciesRxnData( spIx, 3, MW[spIx], thdBdyDefault ) );
    else{
      std::ostringstream msg;
      msg << "Error in " __FILE__ << " : " << __LINE__ << std::endl
          <<" Non-integer reactant stoichiometric coefficient (" << stoich << ") on species " << spNam<< std::endl
          <<" in reaction: " << rxn.equation()  << "\n\n";
      throw std::runtime_error( msg.str() );
    }
  }
  if( reactants.size() < 1 || reactants.size() > 3 ){
    std::ostringstream msg;
    msg << "Error in " __FILE__ << " : " << __LINE__ << std::endl
        << "In reaction: " << rxn.equation() << std::endl
        <<" Number of reactants must be 1 <= n <= 3 " << std::endl
        <<" Number of reactants n = " << reactants.size() << std::endl << std::endl;
    throw std::runtime_error( msg.str() );
  }

  for( const Cantera::Composition::value_type& prod : rxn.products ){
    const std::string& spNam = prod.first;
    const size_t spIx = gas->speciesIndex(spNam);
    const double stoich = prod.second;
    const int istoich = int(stoich);
    if(      fabs( stoich - 1 ) < 1e-10 ) products.push_back( SpeciesRxnData( spIx, -1, MW[spIx], thdBdyDefault ) );
    else if( fabs( stoich - 2 ) < 1e-10 ) products.push_back( SpeciesRxnData( spIx, -2, MW[spIx], thdBdyDefault ) );
    else if( fabs( stoich - 3 ) < 1e-10 ) products.push_back( SpeciesRxnData( spIx, -3, MW[spIx], thdBdyDefault ) );
    else{
      std::ostringstream msg;
      msg << "Error in " __FILE__ << " : " << __LINE__ << std::endl
          <<" Non-integer product stoichiometric coefficient" << std::endl
          <<" Stoichiometric coefficient = " << istoich << std::endl;
      throw std::runtime_error( msg.str() );
    }
  }
  if( rxn.products.size() < 1 || rxn.products.size() > 3 ){
    std::ostringstream msg;
    msg << "Error in " __FILE__ << " : " << __LINE__ << std::endl
        <<" Number of products must be 1 <= n <= 3 " << std::endl
        <<" Number of products n = " << products.size() << std::endl;
    throw std::runtime_error( msg.str() );
  }

  netOrder = 0; // difference of forward and reverse rxn orders, used for equilibrium K
  std::map< int, SpeciesRxnData > netMap; // maps to find duplicates, we will extract non-0 entries to a vector later
  for( const SpeciesRxnData& srd : reactants ){
    netOrder += srd.stoich;
    netMap.insert( std::pair<int,SpeciesRxnData>(srd.index, srd) );
  }
  for( const SpeciesRxnData& srd : products ){
    netOrder += srd.stoich;
    if( netMap.find( srd.index ) == netMap.end() ){
      netMap.insert( std::pair<int,SpeciesRxnData>(srd.index, srd) );
    }
    else{
      const SpeciesRxnData& oldData = netMap.find( srd.index )->second;
      netMap.erase( srd.index );
      if( srd.stoich - oldData.stoich != 0 )
        netMap.insert( std::pair<int,SpeciesRxnData > ( srd.index, SpeciesRxnData( srd.index, srd.stoich + oldData.stoich, srd.mw, srd.thdBdyEff) ));
    }
  }

  std::map< int, SpeciesRxnData >::iterator iNet; // now we keep the net species
  for( iNet = netMap.begin(); iNet != netMap.end(); ++iNet ){
    netSpecies.push_back( iNet->second );
  }

  if( netSpecies.size() < 2 || netSpecies.size() > 5 ){
    std::ostringstream msg;
    msg << "Error in " __FILE__ << " : " << __LINE__ << std::endl
        <<" The number of reacting species must be 2 <= n <= 5 " << std::endl
        <<" Number of active species n = " << netSpecies.size() << std::endl;
    throw std::runtime_error( msg.str() );
  }

}

//====================================================================

CanteraObjects::Setup::
Setup( const std::string transName,
       const std::string inputFileName,
       const std::string inpGroup )
  : transportName( transName ),
    inputFile( inputFileName ),
    inputGroup( inpGroup )
{}

//====================================================================

CanteraObjects::CanteraObjects()
  : gasConstant_( Cantera::GasConstant ),
    hasBeenSetup_( false ),
    hasTransport_( true ) // set to false if no transport is found
{}

//--------------------------------------------------------------------

CanteraObjects&
CanteraObjects::self()
{
  CanteraMutex lock;
  static CanteraObjects s;
  return s;
}

//--------------------------------------------------------------------

CanteraObjects::~CanteraObjects()
{
  while( ! soln_.empty() )  soln_.pop();
  hasBeenSetup_ = false;
}

//--------------------------------------------------------------------

void
CanteraObjects::setup_cantera( const Setup& options,
                               const int ncopies )
{
  CanteraObjects& co = CanteraObjects::self();
  assert( !co.hasBeenSetup_ );
  co.options_ = options;

  if( options.transportName.compare("Mix") != 0 ) co.hasTransport_ = false;

  for( int i=0; i<ncopies; ++i ){
    co.build_new_soln();
  }
  co.extract_thermo_data();
  co.extract_kinetics_data();

  if( co.hasTransport_ ) co.extract_mix_transport_data();

  co.hasBeenSetup_ = true;
}

//--------------------------------------------------------------------

void
CanteraObjects::build_new_soln()
{
  try{
    std::shared_ptr<Cantera::Solution> soln = Cantera::newSolution( options_.inputFile, options_.inputGroup, options_.transportName );
    soln_.push( soln );
//    gas_.push( soln->thermo() );
    if( soln->thermo()->nSpecies() > 1 and hasTransport_ ){
      TransPtr tp = soln->transport();
      assert( tp->transportType().compare("CK_Mix")==0 or tp->transportType().compare("Mix")==0 );
//      trans_.push( tp );
    }
    else{
      hasTransport_ = false;
    }
  }
  catch( Cantera::CanteraError& err ){
    throw std::runtime_error("Error initializing cantera.  Check for proper location of input files.\n" + err.getMessage() );
  }
}

//--------------------------------------------------------------------

void
CanteraObjects::extract_thermo_data()
{
  assert( !hasBeenSetup_ );
  const IdealGasPtr gas = soln_.front()->thermo();
  phaseName_        = gas->name();
  numSpecies_       = gas->nSpecies();
  molecularWeights_ = gas->molecularWeights();
  speciesNames_.clear();
  const Cantera::MultiSpeciesThermo& spThermo = gas->speciesThermo();
  for( int i = 0; i != numSpecies_; ++i ){
    thermDataMap_.insert( std::pair< int, ThermData >( i, ThermData( spThermo, i ) ) );
    speciesNames_.push_back( gas->speciesName(i) );
    speciesIndices_.insert( std::pair< std::string, int>( gas->speciesName(i), i ) );
  }
}

//--------------------------------------------------------------------

void
CanteraObjects::extract_kinetics_data()
{
  assert( !hasBeenSetup_ );
  const SolnPtr soln = soln_.front();
  const auto kin = std::dynamic_pointer_cast<const Cantera::GasKinetics>( soln->kinetics() );
  if( kin == nullptr ){
    numRxns_ = 0;
    return;
  }
  const std::vector< std::shared_ptr<Cantera::Reaction> >& rxnVec = kin->getReactionData(); // contains kinetics data for each reaction
  numRxns_ = rxnVec.size();
  for( size_t r=0; r<numRxns_; ++r){
    RxnData rxnDat( soln->thermo(), *rxnVec[r], molecularWeights_);
    try{
      std::vector<RxnData::SpeciesRxnData>::iterator iThd = rxnDat.thdBdySpecies.begin();
      for( ; iThd != rxnDat.thdBdySpecies.end(); ++iThd) // we evaluate M assuming everything is default, this corrects for non-default species
        iThd->thdBdyEff = iThd->invMW * ( iThd->thdBdyEff - rxnDat.thdBdyDefault );
      rxnDataMap_.insert( std::pair< int, RxnData >( r, rxnDat ) );
    }
    catch( std::runtime_error& err ){
      std::ostringstream msg;
      msg << __FILE__ << " : " << __LINE__
          << " \n\nError occured when parsing data for reaction:\n\t "
          << rxnVec[r]->equation() << "\n\n"
          << err.what() << std::endl;
      throw std::runtime_error( msg.str() );
    }
  }
}

//--------------------------------------------------------------------

void
CanteraObjects::extract_mix_transport_data()
{
  assert( !hasBeenSetup_ );
  assert( hasTransport_ );
  const TransPtr trans = soln_.front()->transport();
  const std::string transType = trans->transportType();
  if( transType.compare("CK_Mix") != 0 and transType.compare("Mix") != 0 ){
    std::ostringstream msg;
    msg << __FILE__ << " : " << __LINE__
        << " \n Unsupported model for transport properties detected " << std::endl
        << " Only mixture average transport properties are currently supported" << std::endl
        << " Supplied model = " << transType << std::endl;
    throw std::runtime_error( msg.str() );
  }
  const auto mixTrans = std::dynamic_pointer_cast<Cantera::MixTransport>( trans ); // cast gas transport object as mix transport
  diffusionCoefs_   = mixTrans->getDiffusionPolyCoefficients(); // diffusion coefficient parameters
  viscosityCoefs_   = mixTrans->getViscosityCoefficients();
  thermalCondCoefs_ = mixTrans->getConductivityCoefficients();
}

//--------------------------------------------------------------------

SolnPtr
CanteraObjects::get_cantera_solution()
{
  CanteraMutex lock;
  CanteraObjects& co = CanteraObjects::self();
  assert( co.hasBeenSetup_ );
  if( co.soln_.size() == 0 ) co.build_new_soln();
  SolnPtr soln = co.soln_.front();
  co.soln_.pop();
  return soln;
}

//--------------------------------------------------------------------

IdealGasPtr
CanteraObjects::get_thermo()
{
  // jcs potential bug. Thermo and Transport objects are linked (at least from transport -> thermo).
  // We can get the thermo object from the transport object, but not vice-versa.
  // That means that calls to restore_thermo will get things out of sync on these queues. This could cause big problems.
  CanteraMutex lock;
  CanteraObjects& co = CanteraObjects::self();
  assert( co.hasBeenSetup_ );
  if( co.soln_.size() == 0 ) co.build_new_soln();
  SolnPtr soln = co.soln_.front();
  co.soln_.pop();
  IdealGasPtr gas = soln->thermo();
  co.gasMapper_[gas] = soln;
  return gas;
}

//--------------------------------------------------------------------

TransPtr
CanteraObjects::get_transport()
{
  CanteraMutex lock;
  CanteraObjects& co = CanteraObjects::self();
  assert( co.hasBeenSetup_ );

  if( !co.hasTransport_ ){
    std::ostringstream msg;
    msg << __FILE__ << " : " << __LINE__
        << " \n\nERROR:\n"
        << "\tget_transport() called when no transport data was available\n"
        << "\tThis is likely because the cantera input file lacked transport data for one or more species\n"
        << "\tThe Cantera error message follows\n"
        << co.transportErrorMessage_.str() << std::endl;
    throw std::runtime_error( msg.str() );
  }

  if( co.soln_.size() == 0 ) co.build_new_soln();
  SolnPtr soln = co.soln_.front();
  co.soln_.pop();
  TransPtr trans = soln->transport();
  co.transMapper_[trans] = soln;
  return trans;
}

//--------------------------------------------------------------------

double
CanteraObjects::gas_constant()
{
  const CanteraObjects& co = CanteraObjects::self();
  assert( co.hasBeenSetup_ );
  return co.gasConstant_;
}

//--------------------------------------------------------------------

double
CanteraObjects::reference_pressure()
{
  CanteraMutex lock;
  IdealGasPtr thermo = CanteraObjects::get_thermo();
  const double rp = thermo->refPressure();
  CanteraObjects::restore_thermo(thermo);
  return rp;
}

//--------------------------------------------------------------------

int
CanteraObjects::number_species()
{
  const CanteraObjects& co = CanteraObjects::self();
  assert( co.hasBeenSetup_ );
  return co.numSpecies_;
}

//--------------------------------------------------------------------

int
CanteraObjects::number_rxns()
{
  const CanteraObjects& co = CanteraObjects::self();
  assert( co.hasBeenSetup_ );
  return co.numRxns_;
}

//--------------------------------------------------------------------

const std::vector<double>&
CanteraObjects::molecular_weights()
{
  const CanteraObjects& co = CanteraObjects::self();
  assert( co.hasBeenSetup_ );
  return co.molecularWeights_;
}

//--------------------------------------------------------------------

std::string
CanteraObjects::species_name( const int i )
{
  const CanteraObjects& co = CanteraObjects::self();
  assert( co.hasBeenSetup_ );
  if( i >= number_species() ){
    std::ostringstream msg;
    msg << "Error in " __FILE__ << " : " << __LINE__ << std::endl
        << " species_name() called for a species that doesn't exist " << std::endl
        << " Species # "<< i << std::endl;
    throw std::runtime_error( msg.str() );
  }
  return co.speciesNames_[i];
}

//--------------------------------------------------------------------

const std::vector<std::string>&
CanteraObjects::species_names()
{
  const CanteraObjects& co = CanteraObjects::self();
  assert( co.hasBeenSetup_ );
  return co.speciesNames_;
}

//--------------------------------------------------------------------

int
CanteraObjects::species_index( const std::string& name )
{
  const CanteraObjects& co = CanteraObjects::self();
  assert( co.hasBeenSetup_ );
  if( co.speciesIndices_.find(name) != co.speciesIndices_.end() )
    return co.speciesIndices_.find(name)->second;
  else{
    std::ostringstream msg;
    msg << "Error in " __FILE__ << " : " << __LINE__ << std::endl
        << " species_index called for a species that doesn't exist " << std::endl
        << " Species name "<< name << std::endl;
    throw std::runtime_error( msg.str() );
  }
}

//--------------------------------------------------------------------

const ThermData&
CanteraObjects::species_thermo( const int i )
{
  const CanteraObjects& co = CanteraObjects::self();
  assert( co.hasBeenSetup_ );
  if( co.thermDataMap_.find(i) != co.thermDataMap_.end() )
    return co.thermDataMap_.find(i)->second;
  else{
    std::ostringstream msg;
    msg << "Error in " __FILE__ << " : " << __LINE__ << std::endl
        <<" species_thermo called for a species that doesn't exist " << std::endl
        <<" Species # "<< i << std::endl;
    throw std::runtime_error( msg.str() );
  }
}

//--------------------------------------------------------------------

const RxnData&
CanteraObjects::rxn_data( const int r )
{
  const CanteraObjects& co = CanteraObjects::self();
  assert( co.hasBeenSetup_ );
  if( co.rxnDataMap_.find(r) != co.rxnDataMap_.end() )
    return co.rxnDataMap_.find(r)->second;
  else{
    std::ostringstream msg;
    msg << "Error in " __FILE__ << " : " << __LINE__ << std::endl
        <<" rxn_data called for a reaction that doesn't exist " << std::endl
        <<" Reaction # "<< r << std::endl;
    throw std::runtime_error( msg.str() );
  }
}

//--------------------------------------------------------------------

const std::vector< std::vector< double > >&
CanteraObjects::diffusion_coefs()
{
  const CanteraObjects& co = CanteraObjects::self();
  assert( co.hasBeenSetup_ );
  return co.diffusionCoefs_;
}

//--------------------------------------------------------------------

const std::vector< std::vector< double > >&
CanteraObjects::viscosity_coefs()
{
  const CanteraObjects& co = CanteraObjects::self();
  assert( co.hasBeenSetup_ );
  return co.viscosityCoefs_;
}

//--------------------------------------------------------------------

const std::vector< std::vector< double > >&
CanteraObjects::thermal_cond_coefs()
{
  const CanteraObjects& co = CanteraObjects::self();
  assert( co.hasBeenSetup_ );
  return co.thermalCondCoefs_;
}

//--------------------------------------------------------------------

const std::string&
CanteraObjects::phase_name()
{
  const CanteraObjects& co = CanteraObjects::self();
  assert( co.hasBeenSetup_ );
  return co.phaseName_;
}

//--------------------------------------------------------------------

void
CanteraObjects::restore_soln( SolnPtr soln )
{
  CanteraMutex lock;
  CanteraObjects& co = CanteraObjects::self();
  assert( co.hasBeenSetup_ );
  co.soln_.push( soln );
}

//--------------------------------------------------------------------

void
CanteraObjects::restore_transport( TransPtr trans )
{
  CanteraMutex lock;
  CanteraObjects& co = CanteraObjects::self();
  assert( co.hasBeenSetup_ );
  assert( co.hasTransport_ );
  co.soln_.push( co.transMapper_[trans] );
}

//--------------------------------------------------------------------

void
CanteraObjects::restore_thermo( IdealGasPtr gas )
{
  CanteraMutex lock;
  CanteraObjects& co = CanteraObjects::self();
  assert( co.hasBeenSetup_ );
  co.soln_.push( co.gasMapper_[gas] );
}

//--------------------------------------------------------------------

