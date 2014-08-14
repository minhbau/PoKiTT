#ifndef CanteraObjects_h
#define CanteraObjects_h

#include <queue>
#include <string>

#include <boost/bimap.hpp>

namespace Cantera_CXX{ class IdealGasMix; }
namespace Cantera    { class Transport;   }

//====================================================================

class CanteraObjects
{
public:

  struct Setup
  {
    std::string transportName,
      inputFile,
      inputGroup;
    Setup() : transportName(""), inputFile(""), inputGroup("") {}
    Setup( const std::string transName,
	   const std::string inputFileName,
	   const std::string inputGroup = "" );
  };

  Cantera_CXX::IdealGasMix* get_gasmix();
  Cantera::Transport* get_transport();

  void restore_transport( Cantera::Transport* const trans );
  void restore_gasmix( Cantera_CXX::IdealGasMix* const gas );

  static CanteraObjects& self();

  void setup_cantera( const Setup& options,
		      const int ncopies = 1 );

private:

  typedef Cantera_CXX::IdealGasMix    IdealGas;
  typedef Cantera::Transport          Trans;
  typedef boost::bimap< IdealGas*, Trans* >  GasTransMap;

  std::queue< std::pair<IdealGas*,Trans*> >  available_;
  GasTransMap gtm_;

  Setup options_;
  bool hasBeenSetup_;

  CanteraObjects();
  ~CanteraObjects();

  CanteraObjects( const CanteraObjects& );
  CanteraObjects& operator=( const CanteraObjects& );

  void build_new( const Setup& options );

};

//====================================================================

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


#endif // CanteraObjects_h
