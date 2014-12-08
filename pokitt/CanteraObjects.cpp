#include <stdexcept>
#include <cassert>

#include "CanteraObjects.h"


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
  CanteraObjects& co = CanteraObjects::self();

  assert( !co.hasBeenSetup_ );
  co.options_ = options;
  for( int i=0; i<ncopies; ++i ) co.build_new();
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
  available_.push( make_pair(gas,trans) );
  gtm_.left.insert( GasTransMap::left_value_type( gas, trans ) );
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

void
CanteraObjects::restore_transport( Cantera::Transport* const trans )
{
  CanteraMutex lock;
  CanteraObjects& co = CanteraObjects::self();
  assert( co.hasBeenSetup_ );
  const GasTransMap::right_iterator itg = co.gtm_.right.find(trans);
  assert( itg != co.gtm_.right.end() );
  co.available_.push( make_pair(itg->second,itg->first) );
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
