/*
 * Thermo_test.cpp
 *
 *  Created on: August 28, 2014
 *      Author: Nathan Yonkee
 */

#include <iostream>
#include <stdio.h>
#include <fstream>
#include "TestHelper.h"

#include <test/TemperaturePowers.h>
#include <pokitt/thermo/HeatCapacity_Cp.h>
#include <pokitt/thermo/HeatCapacity_Cv.h>
#include <pokitt/thermo/Enthalpy.h>

#include <expression/ExprLib.h>

#include <spatialops/structured/Grid.h>
#include <spatialops/structured/FieldComparisons.h>
#include <spatialops/util/TimeLogger.h>

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>

#include <cantera/IdealGasMix.h>

namespace So = SpatialOps;
typedef So::SVolField   CellField;
typedef So::SpatFldPtr<CellField> CellFieldPtrT;

namespace Cantera_CXX{ class IdealGasMix; } //location of polynomial

namespace po = boost::program_options;

using namespace pokitt;

//==============================================================================

enum ThermoQuantity{
  CP,
  CV,
  ENTH
};

std::string thermo_name( const ThermoQuantity q )
{
  switch(q){
    case CP  : return "cp";
    case CV  : return "cv";
    case ENTH: return "h";
  }
}

//==============================================================================

std::set< Expr::ExpressionID >
register_thermo_id( const bool mix,
                    const ThermoQuantity thermoQuantity,
                    Expr::ExpressionFactory& exprFactory,
                    const Expr::TagList& thermoTags,
                    const Expr::Tag& tTag,
                    const Expr::Tag& yiTag,
                    const int nSpec)
{
  std::set< Expr::ExpressionID > thermoID;

  if( mix ){
    typedef HeatCapacity_Cp<CellField>::Builder Cp;
    typedef HeatCapacity_Cv<CellField>::Builder Cv;
    typedef Enthalpy       <CellField>::Builder Enthalpy;
    Expr::ExpressionBuilder* builder = NULL;
    switch( thermoQuantity ){
      case CP  : builder = new       Cp( thermoTags[0], tTag, yiTag ); break;
      case CV  : builder = new       Cv( thermoTags[0], tTag, yiTag ); break;
      case ENTH: builder = new Enthalpy( thermoTags[0], tTag, yiTag ); break;
    } // switch( thermoQuantity )
    thermoID.insert( exprFactory.register_expression(builder) );
  }
  else{ // species
    typedef SpeciesHeatCapacity_Cp<CellField>::Builder SpecCp;
    typedef SpeciesHeatCapacity_Cv<CellField>::Builder SpecCv;
    typedef SpeciesEnthalpy       <CellField>::Builder SpecEnth;
    for( size_t n=0; n<nSpec; ++n ){
      switch( thermoQuantity ){
        case CP  : thermoID.insert(exprFactory.register_expression( new SpecCp  (thermoTags[n], tTag, n) )); break;
        case CV  : thermoID.insert(exprFactory.register_expression( new SpecCv  (thermoTags[n], tTag, n) )); break;
        case ENTH: thermoID.insert(exprFactory.register_expression( new SpecEnth(thermoTags[n], tTag, n) )); break;
      } // switch(thermoQuantity)
    } // species loop
  }
  return thermoID;
}

//==============================================================================

void write_tree( const bool mix, const ThermoQuantity thermoQuantity, const Expr::ExpressionTree& tree)
{
  if( mix ){
    std::ofstream out( (thermo_name(thermoQuantity) + "Mixture.dot").c_str() );
    tree.write_tree(out);
  }
  else{
    std::ofstream out( (thermo_name(thermoQuantity) + "Species.dot").c_str() );
    tree.write_tree(out);
  }
}

//==============================================================================

const std::vector< std::vector<double> >
calculate_mass_fracs(const int nSpec, CellField& xcoord, const Expr::TagList yiTags, Expr::FieldManagerList& fml){
  Expr::FieldMgrSelector<CellField>::type& cellFM = fml.field_manager<CellField>();
  CellField& yi = cellFM.field_ref(yiTags[0]);
  CellFieldPtrT sum = So::SpatialFieldStore::get<CellField>(yi);
  *sum<<=0.0;
  for( size_t n=0; n<nSpec; ++n ){
    CellField& yi = cellFM.field_ref(yiTags[n]);
    yi <<= n + 1 + xcoord;
    *sum <<= *sum + yi;
  }
  BOOST_FOREACH( Expr::Tag yiTag, yiTags){
    CellField& yi = cellFM.field_ref(yiTag);
    yi <<= yi / *sum;
  }
  std::vector< std::vector<double> > massFracs;
  size_t nPts = xcoord.window_with_ghost().glob_npts();
  for( size_t i=0; i<nPts; ++i ){
    std::vector<double> massFrac( nSpec, 0.0 );
    massFracs.push_back(massFrac);
  }
  for( size_t n=0; n<nSpec; ++n ){
    CellField& yi = cellFM.field_ref(yiTags[n]);
#   ifdef ENABLE_CUDA
    yi.set_device_as_active( CPU_INDEX );
#   endif
    size_t i=0;
    for( CellField::iterator iY = yi.begin(); iY != yi.end(); ++iY, ++i ){
      massFracs[i][n] = *iY;
    }
#   ifdef ENABLE_CUDA
    yi.set_device_as_active( GPU_INDEX );
#   endif
  }
  return massFracs;
}

//==============================================================================

std::vector< CellFieldPtrT >
get_cantera_results( const bool mix,
                     const bool timings,
                     const bool canteraReps,
                     const ThermoQuantity thermoQuantity,
                     Cantera_CXX::IdealGasMix& gasMix,
                     const std::vector< std::vector<double> >& massFracs,
                     const CellField& temp)
{
  using namespace SpatialOps;

  const double refPressure = gasMix.pressure();
  const std::vector<double>& molecularWeights = gasMix.molecularWeights();
  const int nSpec = gasMix.nSpecies();

  std::vector< CellFieldPtrT > canteraResults;
  CellField::const_iterator iTemp = temp.begin();
  std::vector< std::vector<double> >::const_iterator iMass = massFracs.begin();

  Timer thermoTimer;

  if( mix ){
    CellFieldPtrT canteraResult = SpatialFieldStore::get<CellField>(temp);
    CellField::iterator iCantEnd = canteraResult->end();

    thermoTimer.start();
    for( size_t rep=0; rep < canteraReps; ++rep ){
      iTemp = temp.begin();
      iMass = massFracs.begin();
      for( CellField::iterator iCant = canteraResult->begin(); iCant!=iCantEnd; ++iTemp, ++iMass, ++iCant ){
        gasMix.setState_TPY( *iTemp, refPressure, &(*iMass)[0]);
        switch(thermoQuantity){
        case CP  : *iCant=gasMix.cp_mass();       break;
        case CV  : *iCant=gasMix.cv_mass();       break;
        case ENTH: *iCant=gasMix.enthalpy_mass(); break;
        } // switch(thermoQuantity)
      }
    }
    thermoTimer.stop();
    canteraResults.push_back( canteraResult );
  }

  else{ // species
    for( size_t n=0; n < nSpec; ++n){
      canteraResults.push_back( SpatialFieldStore::get<CellField>(temp) );
    }
    std::vector<double> thermoResult(nSpec,0.0);
    const CellField::const_iterator iTempEnd = temp.end();
    thermoTimer.start();
    for( size_t rep=0; rep < canteraReps; ++rep ){
      size_t i = 0;
      iMass = massFracs.begin();
      for( iTemp = temp.begin(); iTemp != iTempEnd; ++iTemp, ++iMass, ++i){
        gasMix.setState_TPY( *iTemp, refPressure, &(*iMass)[0]);
        switch( thermoQuantity ){
        case CP  : gasMix.getPartialMolarCp(&thermoResult[0]);         break;
        case CV  : gasMix.getPartialMolarCp(&thermoResult[0]);         break;
        case ENTH: gasMix.getPartialMolarEnthalpies(&thermoResult[0]); break;
        }
        for( size_t n=0; n<nSpec; ++n ){
          (*canteraResults[n])[i] = thermoResult[n];
        }
      }
      for( size_t n=0; n<nSpec; ++n){
        switch( thermoQuantity ){
        case CP  : *canteraResults[n] <<=   *canteraResults[n] / molecularWeights[n];                           break; // convert to mass basis for field comparison
        case CV  : *canteraResults[n] <<= ( *canteraResults[n] - Cantera::GasConstant ) / molecularWeights[n];  break; // convert from molar cp to mass cv for field comparison
        case ENTH: *canteraResults[n] <<=   *canteraResults[n] / molecularWeights[n];                           break; // convert to mass basis for field comparison
        }
      }
    }
    thermoTimer.stop();
  }

  if( timings ) std::cout << "Cantera " + thermo_name(thermoQuantity) + " time " << thermoTimer.elapsed_time()/canteraReps << std::endl;
  return canteraResults;
}

//==============================================================================

bool driver( const bool timings,
             const size_t pokittReps,
             const size_t canteraReps,
             const bool mix,
             const ThermoQuantity thermoQuantity )
{
  TestHelper status( !timings );

  Cantera_CXX::IdealGasMix* const gasMix = CanteraObjects::get_gasmix();

  const int nSpec=gasMix->nSpecies();

  typedef Expr::PlaceHolder < CellField > Temp;
  typedef TemperaturePowers < CellField > TemperaturePowers;
  typedef Expr::PlaceHolder < CellField > MassFracs;

  const Expr::Tag tTag  ( "Temperature", Expr::STATE_NONE );
  const Expr::Tag yiTag ( "yi", Expr::STATE_NONE );
  Expr::TagList yiTags;
  for( size_t n=0; n<nSpec; ++n ){
    std::ostringstream name;
    name << "yi_" << n;
    yiTags.push_back( Expr::Tag( name.str(), Expr::STATE_NONE ) );
  }
  Expr::TagList thermoTags;
  if( mix )
    thermoTags.push_back(Expr::Tag( thermo_name(thermoQuantity) + " mix", Expr::STATE_NONE));
  else{
    for( size_t n=0; n<nSpec; ++n ){
      thermoTags.push_back( Expr::Tag( thermo_name(thermoQuantity) + boost::lexical_cast<std::string>(n), Expr::STATE_NONE ) );
    }
  }

  Expr::ExpressionFactory exprFactory;

  exprFactory.register_expression( new Temp             ::Builder(tTag) );
  exprFactory.register_expression( new TemperaturePowers::Builder(tTag) );
  BOOST_FOREACH( Expr::Tag yiTag, yiTags){
    exprFactory.register_expression( new MassFracs::Builder (yiTag) );
  }
  std::set< Expr::ExpressionID > thermoID = register_thermo_id( mix,
                                                                thermoQuantity,
                                                                exprFactory,
                                                                thermoTags, tTag, yiTag,
                                                                nSpec);

  std::vector<So::IntVec> ptvec;
  if( timings ){
    ptvec.push_back( So::IntVec(  6,  6,  6) );
    ptvec.push_back( So::IntVec( 14, 14, 14) );
    ptvec.push_back( So::IntVec( 30, 30, 30) );
    ptvec.push_back( So::IntVec( 62, 62, 62) );
    ptvec.push_back( So::IntVec(126,126,126) );
  }
  else{
    ptvec.push_back( So::IntVec(20,1,1) );
  }

  for( std::vector<So::IntVec>::iterator iPts = ptvec.begin(); iPts!= ptvec.end(); ++iPts){

    Expr::ExpressionTree thermoTree( thermoID, exprFactory, 0 );
    write_tree( mix, thermoQuantity, thermoTree);

    So::IntVec nPts = *iPts;
    const So::BoundaryCellInfo cellBCInfo = So::BoundaryCellInfo::build<CellField>(false,false,false);
    const So::GhostData cellGhosts(1);
    const So::MemoryWindow vwindow( So::get_window_with_ghost(nPts,cellGhosts,cellBCInfo) );
    CellField xcoord( vwindow, cellBCInfo, cellGhosts, NULL );

    std::vector<double> length(3,1.0);
    So::Grid grid( nPts, length );
    grid.set_coord<SpatialOps::XDIR>( xcoord );
#   ifdef ENABLE_CUDA
    xcoord.add_device( GPU_INDEX );
#   endif

    Expr::FieldManagerList fml;

    thermoTree.register_fields( fml );
    fml.allocate_fields( Expr::FieldAllocInfo( nPts, 0, 0, false, false, false ) );
    thermoTree.bind_fields( fml );

    using namespace SpatialOps;
    Expr::FieldMgrSelector<CellField>::type& cellFM = fml.field_manager< CellField>();

    CellField& temp = cellFM.field_ref(tTag);
    temp <<= 500.0 + 1000.0 * xcoord;

    std::vector< std::vector<double> > massFracs;
    if(mix){
      massFracs = calculate_mass_fracs( nSpec, xcoord, yiTags, fml );
    }
    else{
      for( size_t i=0; i < vwindow.glob_npts(); ++i ){
        std::vector<double> massFrac( nSpec, 1.0 );
        massFracs.push_back(massFrac);
      }
    }

    thermoTree.lock_fields( fml );  // prevent fields from being deallocated so that we can get them after graph execution.

    if( timings ) std::cout << std::endl << thermo_name(thermoQuantity) << " test - " << vwindow.glob_npts() << std::endl;

    Timer thermoTimer;  thermoTimer.start();
    for( size_t rep = 0; rep < pokittReps; ++rep ){
      thermoTree.execute_tree();
    }
    thermoTimer.stop();

    if( timings ) std::cout << "PoKiTT  " + thermo_name(thermoQuantity) + " time " << thermoTimer.elapsed_time()/pokittReps << std::endl;

#   ifdef ENABLE_CUDA
    BOOST_FOREACH( Expr::Tag thermoTag, thermoTags){
      CellField& thermo = fml.field_manager<CellField>().field_ref(thermoTag);
      thermo.set_device_as_active(CPU_INDEX);
    }
    temp.set_device_as_active( CPU_INDEX );
#   endif

    const std::vector< CellFieldPtrT > canteraResults = get_cantera_results( mix,
                                                                             timings,
                                                                             canteraReps,
                                                                             thermoQuantity,
                                                                             *gasMix,
                                                                             massFracs,
                                                                             temp );

    std::vector< CellFieldPtrT >::const_iterator iCantera = canteraResults.begin();
    BOOST_FOREACH( const Expr::Tag& thermoTag, thermoTags ){
      CellField& thermo = cellFM.field_ref(thermoTag);
      switch( thermoQuantity ){
        case CP  :
        case CV  : {
          status( field_equal( thermo, **iCantera, 1e-14 ) || field_equal_abs( thermo, **iCantera, 1e-11 ), thermoTag.name() );
          break;
        }
        case ENTH: {
          status( field_equal( thermo, **iCantera, 1e-11 ) || field_equal_abs( thermo, **iCantera, 1e-8 ), thermoTag.name() );
          break;
        }
      }
      ++iCantera;
    }

  } // number of points

  return status.ok();
}

//==============================================================================

int main( int iarg, char* carg[] )
{
  std::string inputFileName;
  std::string inpGroup;
  bool mix     = false;
  bool timings = false;
  size_t pokittReps = 1;
  size_t canteraReps = 1;

  // parse the command line options input describing the problem
  try {
    po::options_description desc("Supported Options");
    desc.add_options()
           ( "help", "print help message" )
           ( "xml-input-file", po::value<std::string>(&inputFileName), "Cantera xml input file name" )
           ( "phase", po::value<std::string>(&inpGroup), "name of phase in Cantera xml input file" )
           ( "mix", "Triggers mixture heat capacity test.  Otherwise, species heat capacities are tested." )
           ( "timings", "Generate comparison timings between Cantera and PoKiTT across several problem sizes" )
           ( "pokitt-reps", po::value<size_t>(&pokittReps), "Repeat the PoKiTT tests and report the average execution time")
           ( "cantera-reps", po::value<size_t>(&canteraReps), "Repeat the Cantera tests and report the average execution time");

    po::variables_map args;
    po::store( po::parse_command_line(iarg,carg,desc), args );
    po::notify(args);

    mix = args.count("mix") > 0;
    timings = args.count("timings") > 0;

    if (!args.count("xml-input-file")){
      std::cout << "You must enter an xml input file for Cantera" << std::endl;
      std::cout << desc << std::endl;
      return 1;
    }

    if (args.count("help")) {
      std::cout << desc << "\n";
      return -1;
    }
  }

  catch( std::exception& err ){
    std::cout << "Error parsing input arguments\n" << err.what() << std::endl;
    return -2;
  }

  try{
    const CanteraObjects::Setup setup( "Mix", inputFileName, inpGroup );
    CanteraObjects::setup_cantera( setup );

    TestHelper status( !timings );
    status( driver(  timings, pokittReps, canteraReps, mix, CP   ), thermo_name(CP  ) );
    status( driver(  timings, pokittReps, canteraReps, mix, CV   ), thermo_name(CV  ) );
    status( driver(  timings, pokittReps, canteraReps, mix, ENTH ), thermo_name(ENTH) );

    if( status.ok() ){
      std::cout << "\nPASS\n";
      return 0;
    }
  }
  catch( Cantera::CanteraError& ){
    Cantera::showErrors();
  }
  catch( std::exception& err ){
    std::cout << err.what() << std::endl;
  }

  std::cout << "\nFAIL\n";
  return -1;
}

//==============================================================================
