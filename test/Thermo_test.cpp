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

#include <boost/lexical_cast.hpp>
#include <boost/timer.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>

#include <cantera/IdealGasMix.h>

namespace So = SpatialOps;
typedef So::SVolField   CellField;

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
    typedef Enthalpy       <CellField>::Builder Enthalpy;
    typedef HeatCapacity_Cv<CellField>::Builder Cv;
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
    typedef SpeciesEnthalpy       <CellField>::Builder SpecEnth;
    typedef SpeciesHeatCapacity_Cv<CellField>::Builder SpecCv;
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

std::vector< So::SpatFldPtr<CellField> >
get_cantera_results( const bool mix,
                     const bool timings,
                     const ThermoQuantity thermoQuantity,
                     Cantera_CXX::IdealGasMix& gasMix,
                     const int npts,
                     const CellField& temp)
{
  using namespace SpatialOps;

  const double refPressure=gasMix.pressure();
  const std::vector<double>& molecularWeights = gasMix.molecularWeights();
  const int nSpec = gasMix.nSpecies();

  // we need to calculate massfracs manually because they don't get registered for the pure species case
  std::vector< std::vector<double> > massfracs;
  for( size_t i=0; i<npts+2; ++i){
    std::vector<double> massfrac;
    double sum = 0.0;
    for( size_t n=0; n<nSpec; ++n){
      massfrac.push_back(1 + n + (i-0.5)/ npts);
      sum+=massfrac[n];
    }
    for( size_t n=0; n<nSpec; ++n)
      massfrac[n] = massfrac[n]/sum;
    massfracs.push_back(massfrac);
  }

  std::vector< SpatFldPtr<CellField> > canteraResults;
  CellField::const_iterator iTemp = temp.begin();
  std::vector< std::vector<double> >::iterator iMass = massfracs.begin();
  double evalTime;

  if( mix ){
    SpatFldPtr<CellField> canteraResult = SpatialFieldStore::get<CellField>(temp);
    CellField::iterator iCantEnd = canteraResult->end();
    boost::timer thermoTimer;

    for( CellField::iterator iCant = canteraResult->begin(); iCant!=iCantEnd; ++iTemp, ++iMass, ++iCant ){
      gasMix.setState_TPY( *iTemp, refPressure, &(*iMass)[0]);
      switch(thermoQuantity){
        case CP  : *iCant=gasMix.cp_mass();       break;
        case CV  : *iCant=gasMix.cv_mass();       break;
        case ENTH: *iCant=gasMix.enthalpy_mass(); break;
      } // switch(thermoQuantity)
    }
    evalTime = thermoTimer.elapsed();
    canteraResults.push_back( canteraResult );
  }

  else{ // species
    for( size_t n=0; n < nSpec; ++n){
      canteraResults.push_back( SpatialFieldStore::get<CellField>(temp) );
    }
    std::vector<double> thermoResult(nSpec,0.0);
    boost::timer thermoTimer;

    for( size_t i=0; i<npts+2; ++iTemp, ++iMass, ++i){
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
    evalTime = thermoTimer.elapsed();
    for( size_t n=0; n<nSpec; ++n){
      switch( thermoQuantity ){
        case CP  : *canteraResults[n] <<=   *canteraResults[n] / molecularWeights[n];                           break; // convert to mass basis for field comparison
        case CV  : *canteraResults[n] <<= ( *canteraResults[n] - Cantera::GasConstant ) / molecularWeights[n];  break; // convert from molar cp to mass cv for field comparison
        case ENTH: *canteraResults[n] <<=   *canteraResults[n] / molecularWeights[n];                           break; // convert to mass basis for field comparison
      }
    }
  }

  if( timings ) std::cout << "Cantera " + thermo_name(thermoQuantity) + " time " << evalTime << std::endl;
  return canteraResults;
}

//==============================================================================

bool driver( const bool timings,
             const bool mix,
             const ThermoQuantity thermoQuantity )
{
  TestHelper status( !timings );

  Cantera_CXX::IdealGasMix* const gasMix = CanteraObjects::get_gasmix();

  const int nSpec=gasMix->nSpecies();
  size_t n;

  typedef Expr::PlaceHolder < CellField > Temp;
  typedef TemperaturePowers < CellField > TemperaturePowers;
  typedef Expr::PlaceHolder < CellField > MassFracs;

  const Expr::Tag tTag ( "Temperature", Expr::STATE_NONE );
  const Expr::Tag yiTag ( "yi", Expr::STATE_NONE );
  Expr::TagList yiTags;
  for( n=0; n<nSpec; ++n ){
    std::ostringstream name;
    name << "yi_" << n;
    yiTags.push_back( Expr::Tag( name.str(), Expr::STATE_NONE ) );
  }
  Expr::TagList thermoTags;
  if( mix )
    thermoTags.push_back(Expr::Tag( thermo_name(thermoQuantity) + " mix", Expr::STATE_NONE));
  else{
    for( n=0; n<nSpec; ++n ){
      thermoTags.push_back( Expr::Tag( thermo_name(thermoQuantity) + boost::lexical_cast<std::string>(n), Expr::STATE_NONE ) );
    }
  }

  Expr::ExpressionFactory exprFactory;

  exprFactory.register_expression( new Temp::Builder(tTag) );
  exprFactory.register_expression( new TemperaturePowers::Builder(tTag) );
  BOOST_FOREACH( Expr::Tag yiTag, yiTags){
    exprFactory.register_expression( new MassFracs::Builder (yiTag) );
  }
  std::set< Expr::ExpressionID > thermoID = register_thermo_id( mix,
                                                                thermoQuantity,
                                                                exprFactory,
                                                                thermoTags, tTag, yiTag,
                                                                nSpec);

  std::vector<int> ptvec;
  if( timings ){
    ptvec.push_back(8*8*8);
    ptvec.push_back(16*16*16);
    ptvec.push_back(32*32*32);
    ptvec.push_back(64*64*64);
#   ifdef ENABLE_CUDA
    ptvec.push_back(128*128*128);
#   endif
  }
  else{
    ptvec.push_back(10);
  }

  for( std::vector<int>::iterator iPts = ptvec.begin(); iPts!= ptvec.end(); ++iPts){

    Expr::ExpressionTree thermoTree( thermoID, exprFactory, 0 );
    write_tree( mix, thermoQuantity, thermoTree);

    So::IntVec npts(*iPts,1,1);
    const So::BoundaryCellInfo cellBCInfo = So::BoundaryCellInfo::build<CellField>(false,false,false);
    const So::GhostData cellGhosts(1);
    const So::MemoryWindow vwindow( So::get_window_with_ghost(npts,cellGhosts,cellBCInfo) );
    CellField xcoord( vwindow, cellBCInfo, cellGhosts, NULL );

    std::vector<double> length(3,1.0);
    So::Grid grid( npts, length );
    grid.set_coord<SpatialOps::XDIR>( xcoord );
#   ifdef ENABLE_CUDA
    xcoord.add_device( GPU_INDEX );
#   endif

    Expr::FieldManagerList fml;

    thermoTree.register_fields( fml );
    fml.allocate_fields( Expr::FieldAllocInfo( npts, 0, 0, false, false, false ) );
    thermoTree.bind_fields( fml );

    using namespace SpatialOps;
    Expr::FieldMgrSelector<CellField>::type& cellFM = fml.field_manager< CellField>();

    CellField& temp = cellFM.field_ref(tTag);
    temp <<= 500.0 + 1000.0 * xcoord;

    if(mix){
      SpatFldPtr<CellField> sum  = SpatialFieldStore::get<CellField>(temp);
      *sum<<=0.0;
      for( n=0; n<nSpec; ++n ){
        CellField& yi = cellFM.field_ref(yiTags[n]);
        yi <<= n + 1 + xcoord;
        *sum <<= *sum + yi;
      }
      BOOST_FOREACH( Expr::Tag yiTag, yiTags){
        CellField& yi = cellFM.field_ref(yiTag);
        yi <<= yi / *sum;
      }
    }

    thermoTree.lock_fields( fml );  // prevent fields from being deallocated so that we can get them after graph execution.

    if( timings ) std::cout << std::endl << thermo_name(thermoQuantity) << " test - " << *iPts << std::endl;

    boost::timer thermoTimer;
    thermoTree.execute_tree();

    if( timings ) std::cout << "PoKiTT  " + thermo_name(thermoQuantity) + " time " << thermoTimer.elapsed() << std::endl;

#   ifdef ENABLE_CUDA
    BOOST_FOREACH( Expr::Tag thermoTag, thermoTags){
      CellField& thermo = fml.field_manager<CellField>().field_ref(thermoTag);
      thermo.set_device_as_active(CPU_INDEX);
    }
    temp.set_device_as_active( CPU_INDEX );
#   endif

    const std::vector< SpatFldPtr<CellField> > canteraResults = get_cantera_results( mix,
                                                                                     timings,
                                                                                     thermoQuantity,
                                                                                     *gasMix,
                                                                                     *iPts,
                                                                                     temp );

    std::vector< SpatFldPtr<CellField> >::const_iterator iCantera = canteraResults.begin();
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

  // parse the command line options input describing the problem
  try {
    po::options_description desc("Supported Options");
    desc.add_options()
           ( "help", "print help message" )
           ( "xml-input-file", po::value<std::string>(&inputFileName), "Cantera xml input file name" )
           ( "phase", po::value<std::string>(&inpGroup), "name of phase in Cantera xml input file" )
           ( "mix", "Triggers mixture heat capacity test.  Otherwise, species heat capacities are tested." )
           ( "timings", "Generate comparison timings between Cantera and PoKiTT across several problem sizes" );

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
    status( driver(  timings, mix, CP   ), thermo_name(CP  ) );
    status( driver(  timings, mix, CV   ), thermo_name(CV  ) );
    status( driver(  timings, mix, ENTH ), thermo_name(ENTH) );

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
