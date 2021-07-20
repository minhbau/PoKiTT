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

#ifndef CanteraObjects_h
#define CanteraObjects_h

#include <queue>
#include <vector>
#include <map>
#include <string>

#include <boost/bimap.hpp>

#include <cantera/kinetics/Reaction.h>

//====================================================================

namespace Cantera{
  class Transport;
  class ThermoPhase;
  class Solution;
  class MultiSpeciesThermo;
}

typedef std::shared_ptr<Cantera::Solution   >  SolnPtr;
typedef std::shared_ptr<Cantera::ThermoPhase>  IdealGasPtr;
typedef std::shared_ptr<Cantera::Transport  >  TransPtr;

enum ThermoPoly
{
  UNKNOWN_POLY = 0,
  CONST_POLY,
  NASA_POLY
};

struct ThermData
{
  ThermData( const Cantera::MultiSpeciesThermo& spThermo, const int i );
  int index;
  double minTemp;
  double maxTemp;
  ThermoPoly type;
  std::vector< double > coefficients;
};

enum ReactionType{
  UNKNOWN_RXN = 0,
  ELEMENTARY,
  THIRD_BODY,
  LINDEMANN,
  TROE
};

struct RxnData
{
  RxnData( const IdealGasPtr gas,
           const Cantera::Reaction& rxn,
           const std::vector<double>& MW );

  struct SpeciesRxnData{
    SpeciesRxnData( int index,
                    int stoich,
                    double mw,
                    double thdBdyEff );
    int index;
    int stoich;
    double mw;
    double invMW;
    double thdBdyEff;
  };

  const bool reversible;
  std::vector< SpeciesRxnData > reactants;
  std::vector< SpeciesRxnData > products;
  std::vector< SpeciesRxnData > netSpecies;
  std::vector< SpeciesRxnData > thdBdySpecies;
  double kFwdCoefs[3];
  double kPressureCoefs[3];
  double troeParams[4];
  double thdBdyDefault;
  std::map<std::string,double> tbeff;
  int netOrder;
  ReactionType type;
};


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

  static IdealGasPtr get_thermo();
  static TransPtr get_transport();
  static SolnPtr get_cantera_solution();

  static void restore_soln( SolnPtr soln );
  static void restore_transport( TransPtr trans );
  static void restore_thermo(IdealGasPtr gas );

  static void setup_cantera( const Setup& options,
                             const int ncopies = 1 );

  static const std::string& phase_name();
  static double reference_pressure();
  static double gas_constant();
  static int number_species();
  static int number_rxns();
  static const std::vector< double >& molecular_weights();

  static std::string species_name( const int i );
  static const std::vector<std::string>& species_names();
  static int species_index( const std::string& name );
  static const ThermData& species_thermo( const int i );
  static const RxnData& rxn_data( const int r );

  typedef std::vector< std::vector< double > > TransportCoefsT;
  static const TransportCoefsT& diffusion_coefs();
  static const TransportCoefsT& viscosity_coefs();
  static const TransportCoefsT& thermal_cond_coefs();

private:

  static CanteraObjects& self();

  bool hasTransport_;

  std::queue<SolnPtr> soln_;
  std::map<IdealGasPtr,SolnPtr> gasMapper_;
  std::map<TransPtr,SolnPtr> transMapper_;

  std::string phaseName_;
  const double gasConstant_;
  int numSpecies_;
  int numRxns_;
  std::vector< double > molecularWeights_;

  std::vector< std::string > speciesNames_;
  std::map< std::string, int > speciesIndices_;
  std::map< int, ThermData > thermDataMap_;
  std::map< int, RxnData > rxnDataMap_;
  TransportCoefsT diffusionCoefs_;
  TransportCoefsT viscosityCoefs_;
  TransportCoefsT thermalCondCoefs_;

  Setup options_;
  bool hasBeenSetup_;

  std::ostringstream transportErrorMessage_;

  CanteraObjects();
  ~CanteraObjects();

  CanteraObjects( const CanteraObjects& );
  CanteraObjects& operator=( const CanteraObjects& );

  void build_new_soln();
  void extract_thermo_data();
  void extract_kinetics_data();
  void extract_mix_transport_data();
};

//====================================================================

#endif // CanteraObjects_h
