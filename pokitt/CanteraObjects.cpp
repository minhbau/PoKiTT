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

#include <stdexcept>
#include <cassert>

#include "CanteraObjects.h"

#include <cantera/transport.h>
#include <cantera/IdealGasMix.h>

#include <cantera/base/ct_defs.h>              // value of gas constant
#include <cantera/thermo/speciesThermoTypes.h> // definitions for which polynomial is being used
#include <cantera/kinetics/reaction_defs.h>    // reaction type definitions

#include <boost/foreach.hpp>

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

ThermData::ThermData( const Cantera::SpeciesThermo& spThermo, const int i )
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

RxnData::RxnData( const Cantera::IdealGasMix& gas,
                  const Cantera::Reaction& rxn,
                  const std::vector<double>& MW )
  : reversible( rxn.reversible )
{
  const int rxnType = rxn.reaction_type;

  for( int i=0; i<3; ++i ){
    kFwdCoefs     [i] = 0.0;
    kPressureCoefs[i] = 0.0;
  }
  for( int i=0; i<4; ++i ) troeParams[i]=0;

  // ensure that the type of each reaction is supported by PoKiTT.
  if( rxnType != Cantera::ELEMENTARY_RXN &&
      rxnType != Cantera::THREE_BODY_RXN &&
      rxnType != Cantera::FALLOFF_RXN )
  {
    std::ostringstream msg;
    msg << __FILE__ << " : " << __LINE__
        << "\nUnsupported reaction type encountered for \n"
        << "\t" << rxn.equation() << "\n\t(Rxn type: " << rxnType << ")\n";
    throw std::runtime_error( msg.str() );
  }

  switch( rxnType ){
    case Cantera::ELEMENTARY_RXN:{
      const Cantera::ElementaryReaction& elemRxn = dynamic_cast<const Cantera::ElementaryReaction&>(rxn);
      kFwdCoefs[0] = elemRxn.rate.preExponentialFactor();
      kFwdCoefs[1] = elemRxn.rate.temperatureExponent();
      kFwdCoefs[2] = elemRxn.rate.activationEnergy_R();
      type = ELEMENTARY;
      break;
    }
    case Cantera::FALLOFF_RXN:{
      const Cantera::FalloffReaction& fallRxn = dynamic_cast<const Cantera::FalloffReaction&>(rxn);
      kFwdCoefs     [0] = fallRxn.high_rate.preExponentialFactor();
      kFwdCoefs     [1] = fallRxn.high_rate.temperatureExponent();
      kFwdCoefs     [2] = fallRxn.high_rate.activationEnergy_R();
      kPressureCoefs[0] = fallRxn.low_rate.preExponentialFactor();
      kPressureCoefs[1] = fallRxn.low_rate.temperatureExponent();
      kPressureCoefs[2] = fallRxn.low_rate.activationEnergy_R();
      fallRxn.falloff->getParameters( troeParams );
      thdBdyDefault = fallRxn.third_body.default_efficiency;
      BOOST_FOREACH( const Cantera::Composition::value_type& vt, fallRxn.third_body.efficiencies ){
        const int spIx = gas.speciesIndex(vt.first);
        thdBdySpecies.push_back( SpeciesRxnData( spIx, 0, MW[spIx], vt.second ) );
      }

      switch( fallRxn.falloff->getType() ){
        case Cantera::SIMPLE_FALLOFF: type = LINDEMANN; break;
        case Cantera::TROE_FALLOFF  : type = TROE;      break;
        case Cantera::SRI_FALLOFF   : type = TROE;      break;
        default:{
          std::ostringstream msg;
          msg << __FILE__ << " : " << __LINE__
              << "\nFalloff type not supported,\n Type = " << fallRxn.falloff->getType()
              << "\nSee reaction_defs.h in Cantera\n";
          throw std::runtime_error( msg.str() );
          break;
        }
      }
      break;
    }
    case Cantera::THREE_BODY_RXN:{
      const Cantera::ThreeBodyReaction& tbRxn = dynamic_cast<const Cantera::ThreeBodyReaction&>(rxn);
      kFwdCoefs[0]      = tbRxn.rate.preExponentialFactor();
      kFwdCoefs[1]      = tbRxn.rate.temperatureExponent();
      kFwdCoefs[2]      = tbRxn.rate.activationEnergy_R();
      thdBdyDefault     = tbRxn.third_body.default_efficiency;
      type              = THIRD_BODY;
      BOOST_FOREACH( const Cantera::Composition::value_type& vt, tbRxn.third_body.efficiencies ){
        const int spIx = gas.speciesIndex(vt.first);
        thdBdySpecies.push_back( SpeciesRxnData( spIx, 0, MW[spIx], vt.second ) );
      }
      break;
    }
  }

  /* Here we are reading Cantera's product and reactant stoichiometric coefficients and storing them as ints
   * This helps performance because we only support elementary reaction mechanisms
   * We also check if stoichiometry is consistent with elementary reactions
   */
  BOOST_FOREACH( const Cantera::Composition::value_type& reactant, rxn.reactants ){
    const std::string& spNam = reactant.first;
    const size_t spIx = gas.speciesIndex(spNam);
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

  BOOST_FOREACH( const Cantera::Composition::value_type& prod, rxn.products ){
    const std::string& spNam = prod.first;
    const size_t spIx = gas.speciesIndex(spNam);
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
  BOOST_FOREACH( const SpeciesRxnData& srd, reactants ){
    netOrder += srd.stoich;
    netMap.insert( std::pair<int,SpeciesRxnData>(srd.index, srd) );
  }
  BOOST_FOREACH( const SpeciesRxnData& srd, products ){
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
{
}

//====================================================================

CanteraObjects::CanteraObjects()
  : gasConstant_( Cantera::GasConstant ),
    hasBeenSetup_( false )
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
  while( !available_.empty() ){
    std::pair<IdealGas*,Trans*> p = available_.front();
    delete p.first;
    delete p.second;
    available_.pop();
  }
  gtm_.clear();
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
  for( int i=0; i<ncopies; ++i ) co.build_new();
  co.extract_thermo_data();
  co.extract_kinetics_data();
  co.extract_mix_transport_data();
  co.hasBeenSetup_ = true;
}

//--------------------------------------------------------------------

void
CanteraObjects::build_new()
{
  Cantera::IdealGasMix * gas   = NULL;
  Cantera::Transport       * trans = NULL;
  try{
    gas = new Cantera::IdealGasMix( options_.inputFile, options_.inputGroup ) ;
    trans = Cantera::TransportFactory::factory()->newTransport( options_.transportName, gas );
  }
  catch( Cantera::CanteraError& ){
    Cantera::showErrors();
    throw std::runtime_error("Error initializing cantera.  Check for proper location of input files.\n");
  }
  available_.push( std::make_pair(gas,trans) );
  gtm_.left.insert( GasTransMap::left_value_type( gas, trans ) );
}

//--------------------------------------------------------------------

void
CanteraObjects::extract_thermo_data()
{
  assert( !hasBeenSetup_ );
  IdealGas* gas = available_.front().first;
  phaseName_ = gas->name();
  numSpecies_ = gas->nSpecies();
  molecularWeights_ = gas->molecularWeights();
  const Cantera::SpeciesThermo& spThermo = gas->speciesThermo();
  for( int i = 0; i != numSpecies_; ++i ){
    thermDataMap_.insert( std::pair< int, ThermData >( i, ThermData( spThermo, i ) ) );
    speciesNames_.insert( std::pair< int, std::string>( i, gas->speciesName(i) ) );
    speciesIndices_.insert( std::pair< std::string, int>( gas->speciesName(i), i ) );
  }
}

//--------------------------------------------------------------------

void
CanteraObjects::extract_kinetics_data()
{
  assert( !hasBeenSetup_ );
  IdealGas* gas = available_.front().first;

  const std::vector< Cantera::shared_ptr<Cantera::Reaction> >& rxnVec = gas->getReactionData(); // contains kinetics data for each reaction
  numRxns_ = rxnVec.size();
  for( size_t r=0; r<numRxns_; ++r){
    RxnData rxnDat( *gas, *rxnVec[r], molecularWeights_);
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
  Trans* trans = available_.front().second;
  if( trans->model() != Cantera::cMixtureAveraged){
    std::ostringstream msg;
    msg << __FILE__ << " : " << __LINE__
        << " \n Unsupported model for transport properties detected " << std::endl
        << " Only mixture average transport properties are currently supported" << std::endl
        << " Supplied model = " << trans->model() << std::endl
        << " See Cantera::TransportBase.h for model definitions" << std::endl;
    throw std::runtime_error( msg.str() );
  }
  Cantera::MixTransport* mixTrans = dynamic_cast<Cantera::MixTransport*>( trans ); // cast gas transport object as mix transport
  diffusionCoefs_ = mixTrans->getDiffusionPolyCoefficients(); // diffusion coefficient parameters
  viscosityCoefs_ = mixTrans->getViscosityCoefficients();
  thermalCondCoefs_ = mixTrans->getConductivityCoefficients();
}

//--------------------------------------------------------------------

Cantera::IdealGasMix*
CanteraObjects::get_gasmix()
{
  CanteraMutex lock;
  CanteraObjects& co = CanteraObjects::self();

  assert( co.hasBeenSetup_ );
  if( co.available_.size() == 0 ) co.build_new();
  IdealGas* ig = co.available_.front().first;
  co.available_.pop();
  return ig;
}

//--------------------------------------------------------------------

Cantera::Transport*
CanteraObjects::get_transport()
{
  CanteraMutex lock;
  CanteraObjects& co = CanteraObjects::self();
  assert( co.hasBeenSetup_ );
  if( co.available_.size() == 0 ) co.build_new();
  Trans* tr = co.available_.front().second;
  co.available_.pop();
  return tr;
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
  Cantera::IdealGasMix* thermo = CanteraObjects::get_gasmix();
  const double rp = thermo->refPressure();
  CanteraObjects::restore_gasmix(thermo);
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
  if( co.speciesNames_.find(i) != co.speciesNames_.end() )
    return co.speciesNames_.find(i)->second;
  else{
    std::ostringstream msg;
    msg << "Error in " __FILE__ << " : " << __LINE__ << std::endl
        << " species_name called for a species that doesn't exist " << std::endl
        << " Species # "<< i << std::endl;
    throw std::runtime_error( msg.str() );
  }
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
CanteraObjects::restore_transport( Cantera::Transport* const trans )
{
  CanteraMutex lock;
  CanteraObjects& co = CanteraObjects::self();
  assert( co.hasBeenSetup_ );
  const GasTransMap::right_iterator itg = co.gtm_.right.find(trans);
  assert( itg != co.gtm_.right.end() );
  co.available_.push( std::make_pair(itg->second,itg->first) );
}

//--------------------------------------------------------------------

void
CanteraObjects::restore_gasmix( Cantera::IdealGasMix* const gas )
{
  CanteraMutex lock;
  CanteraObjects& co = CanteraObjects::self();
  assert( co.hasBeenSetup_ );
  const GasTransMap::left_iterator igt = co.gtm_.left.find(gas);
  assert( igt != co.gtm_.left.end() );
  std::pair<IdealGas*,Trans*> pp( igt->first, igt->second );
  co.available_.push( pp );
}

//--------------------------------------------------------------------

