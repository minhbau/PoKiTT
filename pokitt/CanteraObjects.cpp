#include <stdexcept>
#include <cassert>

#include "CanteraObjects.h"

#include <cantera/Cantera.h>
#include <cantera/transport.h>
#include <cantera/IdealGasMix.h>

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
  assert( !hasBeenSetup_ );
  for( int i=0; i<ncopies; ++i ){
    build_new(options);
  }
  hasBeenSetup_ = true;
  options_ = options;
}

//--------------------------------------------------------------------

void
CanteraObjects::build_new( const Setup& options )
{
  Cantera_CXX::IdealGasMix * gas   = NULL;
  Cantera::Transport       * trans = NULL;
  try{
    gas = new Cantera_CXX::IdealGasMix( options.inputFile, options.inputGroup ) ;
    trans = Cantera::TransportFactory::factory()->newTransport( options.transportName, gas );
  }
  catch( Cantera::CanteraError& ){
    Cantera::showErrors();
    throw std::runtime_error("Error initializing cantera.  Check for proper location of input files.\n");
  }
  available_.push( make_pair(gas,trans) );
  gtm_.left.insert( GasTransMap::left_value_type( gas, trans ) );
}

//--------------------------------------------------------------------

Cantera_CXX::IdealGasMix*
CanteraObjects::get_gasmix()
{
  CanteraMutex lock;
  assert( hasBeenSetup_ );
  if( available_.size() == 0 ) build_new( options_ );
  IdealGas* ig = available_.front().first;
  available_.pop();
  return ig;
}

//--------------------------------------------------------------------

Cantera::Transport*
CanteraObjects::get_transport()
{
  CanteraMutex lock;
  assert( hasBeenSetup_ );
  if( available_.size() == 0 ) build_new( options_ );
  Trans* tr = available_.front().second;
  available_.pop();
  return tr;
}

//--------------------------------------------------------------------

void
CanteraObjects::restore_transport( Cantera::Transport* const trans )
{
  CanteraMutex lock;
  assert( hasBeenSetup_ );
  const GasTransMap::right_iterator itg = gtm_.right.find(trans);
  assert( itg != gtm_.right.end() );
  available_.push( make_pair(itg->second,itg->first) );
}

//--------------------------------------------------------------------

void
CanteraObjects::restore_gasmix( Cantera_CXX::IdealGasMix* const gas )
{
  CanteraMutex lock;
  assert( hasBeenSetup_ );
  const GasTransMap::left_iterator igt = gtm_.left.find(gas);
  assert( igt != gtm_.left.end() );
  std::pair<IdealGas*,Trans*> pp( igt->first, igt->second );
  available_.push( pp );
}

//--------------------------------------------------------------------
