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

  static Cantera_CXX::IdealGasMix* get_gasmix();
  static Cantera::Transport* get_transport();

  static void restore_transport( Cantera::Transport* const trans );
  static void restore_gasmix( Cantera_CXX::IdealGasMix* const gas );


  static void setup_cantera( const Setup& options,
                             const int ncopies = 1 );

private:

  static CanteraObjects& self();

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

  void build_new();

};

//====================================================================

#endif // CanteraObjects_h
