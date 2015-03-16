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

namespace pokitt{
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

ThermData::ThermData( const Cantera::SpeciesThermo& spThermo, const int i ) : index( i ) {
    double refPressure; int cType;
    coefficients.resize( 15 );
    spThermo.reportParams(i, cType, &(coefficients[0]), minTemp, maxTemp, refPressure);
    switch( cType ){
    case SIMPLE:   type = CONST_POLY;   break;
    case NASA2:    type = NASA_POLY;    break;
    case SHOMATE2: type = SHOMATE_POLY; break;
    default:{
      std::ostringstream msg;
      msg << __FILE__ << " : " << __LINE__
          << "\nThermo type not supported,\n Type = " << type
          << ", species # " << i << std::endl
          << "See speciesThermoTypes.h in Cantera\n";
      throw std::runtime_error( msg.str() );
      }
    }
}

RxnData::SpeciesRxnData::SpeciesRxnData( const int index, const int stoich, const double mw,  const double thdBdyEff )
        : index( index ), stoich( stoich ), mw( mw ), invMW( 1/mw ), thdBdyEff( thdBdyEff ) {}

RxnData::RxnData( const Cantera::ReactionData& cDat, const std::vector<double>& MW )  {
  kFwdCoefs      = cDat.rateCoeffParameters;
  kPressureCoefs = cDat.auxRateCoeffParameters;
  thdBdyDefault  = cDat.default_3b_eff;
  troeParams     = cDat.falloffParameters;
  reversible     = cDat.reversible;

  /* Here we are reading Cantera's product and reactant stoichiometric coefficients and storing them as ints
   * This helps performance because we only support elementary reaction mechanisms
   * We also check if stoichiometry is consistent with elementary reactions
   */
  std::vector< double >::const_iterator iSto;
  std::vector< int >::const_iterator iInd = cDat.reactants.begin();
  for( iSto = cDat.rstoich.begin(); iSto!=cDat.rstoich.end(); ++iSto, ++iInd ){
    if(      fabs( *iSto - 1 ) < 1e-2 ) reactants.push_back( SpeciesRxnData( *iInd, 1, MW[*iInd], thdBdyDefault ) );
    else if( fabs( *iSto - 2 ) < 1e-2 ) reactants.push_back( SpeciesRxnData( *iInd, 2, MW[*iInd], thdBdyDefault ) );
    else if( fabs( *iSto - 3 ) < 1e-2 ) reactants.push_back( SpeciesRxnData( *iInd, 3, MW[*iInd], thdBdyDefault ) );
    else{
      std::ostringstream msg;
      msg << "Error in " __FILE__ << " : " << __LINE__ << std::endl
          <<" Non-integer reactant stoichiometric coefficient" << std::endl
          <<" Stoichiometric coefficient = " << *iSto << std::endl;
      throw std::runtime_error( msg.str() );
    }
  }
  if( reactants.size() < 1 || reactants.size() > 3 ){
    std::ostringstream msg;
    msg << "Error in " __FILE__ << " : " << __LINE__ << std::endl
        <<" Number of reactants must be 1 <= n <= 3 " << std::endl
        <<" Number of reactants n = " << reactants.size() << std::endl;
    throw std::runtime_error( msg.str() );
  }

  iInd = cDat.products.begin();
  for( iSto = cDat.pstoich.begin(); iSto!=cDat.pstoich.end(); ++iSto, ++iInd ){
    if(      fabs( *iSto - 1 ) < 1e-2 ) products.push_back( SpeciesRxnData( *iInd, -1, MW[*iInd], thdBdyDefault ) );
    else if( fabs( *iSto - 2 ) < 1e-2 ) products.push_back( SpeciesRxnData( *iInd, -2, MW[*iInd], thdBdyDefault ) );
    else if( fabs( *iSto - 3 ) < 1e-2 ) products.push_back( SpeciesRxnData( *iInd, -3, MW[*iInd], thdBdyDefault ) );
    else{
      std::ostringstream msg;
      msg << "Error in " __FILE__ << " : " << __LINE__ << std::endl
          <<" Non-integer product stoichiometric coefficient" << std::endl
          <<" Stoichiometric coefficient = " << *iSto << std::endl;
      throw std::runtime_error( msg.str() );
    }
  }
  if( products.size() < 1 || products.size() > 3 ){
    std::ostringstream msg;
    msg << "Error in " __FILE__ << " : " << __LINE__ << std::endl
        <<" Number of products must be 1 <= n <= 3 " << std::endl
        <<" Number of products n = " << products.size() << std::endl;
    throw std::runtime_error( msg.str() );
  }

  netOrder = 0; // difference of forward and reverse rxn orders, used for equilibrium K
   std::map< int, SpeciesRxnData > netMap; // maps to find duplicates, we will extract non-0 entries to a vector later
   std::vector< SpeciesRxnData >::const_iterator iR = reactants.begin();
   for(  ; iR != reactants.end(); ++iR ){
     netOrder += iR->stoich;
     netMap.insert( std::pair< int, SpeciesRxnData>(iR->index, *iR) );
   }
   for(iR = products.begin(); iR != products.end(); ++iR ){
     netOrder += iR->stoich;
     if( netMap.find( iR->index ) == netMap.end() ){
       netMap.insert( std::pair< int, SpeciesRxnData>(iR->index, *iR) );
     }
     else{
       SpeciesRxnData oldData = netMap.find( iR->index )->second;
       netMap.erase( iR->index );
       if( iR->stoich - oldData.stoich != 0 )
         netMap.insert( std::pair< int, SpeciesRxnData > ( iR->index, SpeciesRxnData( iR->index, iR->stoich + oldData.stoich, iR->mw, iR->thdBdyEff) ));
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

   std::map< int, double >::const_iterator iMap;
   for( iMap = cDat.thirdBodyEfficiencies.begin(); iMap != cDat.thirdBodyEfficiencies.end(); ++iMap )
     thdBdySpecies.push_back( SpeciesRxnData( iMap->first, 0, MW[iMap->first], iMap->second ) );

   // check if pressure dependant reaction or not
  switch( cDat.reactionType ){
  case Cantera::ELEMENTARY_RXN: type = ELEMENTARY; break;
  case Cantera::THREE_BODY_RXN: type = THIRD_BODY; break;
  case Cantera::FALLOFF_RXN:{
    switch( cDat.falloffType ){
    case Cantera::SIMPLE_FALLOFF: type = LINDEMANN; break;
    case Cantera::TROE3_FALLOFF:  type = TROE; break;
    case Cantera::TROE4_FALLOFF:  type = TROE; break;
    default:{
      std::ostringstream msg;
      msg << __FILE__ << " : " << __LINE__
          << "\nFalloff type not supported,\n Type = " << cDat.falloffType
          << "See reaction_defs.h in Cantera\n";
      throw std::runtime_error( msg.str() );
      break;
    }
    }
    break;
  }
    default:{
      std::ostringstream msg;
      msg << __FILE__ << " : " << __LINE__
          << "\nKinetics type not supported,\n Type = " << cDat.reactionType
          << "See reaction_defs.h in Cantera\n";
      throw std::runtime_error( msg.str() );
    }
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
  : gasConstant_( Cantera::GasConstant )
{
  hasBeenSetup_ = false;
}

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
  co.hasBeenSetup_ = true;
}

//--------------------------------------------------------------------

void
CanteraObjects::build_new()
{
  Cantera_CXX::IdealGasMix * gas   = NULL;
  Cantera::Transport       * trans = NULL;
  try{
    gas = new Cantera_CXX::IdealGasMix( options_.inputFile, options_.inputGroup ) ;
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
  numSpecies_ = gas->nSpecies();
  molecularWeights_ = gas->molecularWeights();
  const Cantera::SpeciesThermo& spThermo = gas->speciesThermo();
  for( int i = 0; i != numSpecies_; ++i ){
    thermDataMap_.insert( std::pair< int, ThermData >( i, ThermData( spThermo, i ) ) );
  }
}

//--------------------------------------------------------------------

void
CanteraObjects::extract_kinetics_data()
{
  assert( !hasBeenSetup_ );
  IdealGas* gas = available_.front().first;

  const std::vector<Cantera::ReactionData>& canteraDataVec = gas->getReactionData(); // contains kinetics data for each reaction
  numRxns_ = canteraDataVec.size();
  for( size_t r=0; r<numRxns_; ++r){
    const Cantera::ReactionData& canteraData = canteraDataVec[r]; // ReactionData for reaction r
    try{
      rxnDataMap_.insert( std::pair< int, RxnData >( r, RxnData( canteraData, molecularWeights_) ) );
    }
    catch( std::runtime_error& err ){
      std::ostringstream msg;
      msg << " \n Error occured when parsing data for reaction " << r << std::endl
          << err.what() << std::endl;
      throw std::runtime_error( msg.str() );
    }
  }
}

//--------------------------------------------------------------------

Cantera_CXX::IdealGasMix*
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
  CanteraMutex lock;
  CanteraObjects& co = CanteraObjects::self();
  assert( co.hasBeenSetup_ );
  return co.gasConstant_;
}

//--------------------------------------------------------------------

int
CanteraObjects::number_species()
{
  CanteraMutex lock;
  CanteraObjects& co = CanteraObjects::self();
  assert( co.hasBeenSetup_ );
  return co.numSpecies_;
}

//--------------------------------------------------------------------

int
CanteraObjects::number_rxns()
{
  CanteraMutex lock;
  CanteraObjects& co = CanteraObjects::self();
  assert( co.hasBeenSetup_ );
  return co.numRxns_;
}

//--------------------------------------------------------------------

const std::vector<double>&
CanteraObjects::molecular_weights()
{
  CanteraMutex lock;
  CanteraObjects& co = CanteraObjects::self();
  assert( co.hasBeenSetup_ );
  return co.molecularWeights_;
}

//--------------------------------------------------------------------

const ThermData&
CanteraObjects::species_thermo( const int i )
{
  CanteraMutex lock;
  CanteraObjects& co = CanteraObjects::self();
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
  CanteraMutex lock;
  CanteraObjects& co = CanteraObjects::self();
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
CanteraObjects::restore_gasmix( Cantera_CXX::IdealGasMix* const gas )
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
} //namespace pokitt
