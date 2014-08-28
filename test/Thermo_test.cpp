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

#include <pokitt/thermo/TemperaturePowers.h>
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

//==============================================================================

std::set< Expr::ExpressionID > register_thermo_id( const bool mix,
                                                   const std::string& thermo_quantity,
                                                   Expr::ExpressionFactory& exprFactory,
                                                   const Expr::TagList& thermoTags, const Expr::Tag& tTag, const Expr::Tag& yiTag,
                                                   const int nSpec)
{
  std::set< Expr::ExpressionID > thermo_id;
  if( mix )
  {
    if(      thermo_quantity == "cp" ){
      typedef HeatCapacity_Cp < CellField > HeatCapacity;
      thermo_id.insert(exprFactory.register_expression( new HeatCapacity::Builder( thermoTags[0], tTag, yiTag )));
    }
    else if( thermo_quantity == "cv" ){
      typedef HeatCapacity_Cv < CellField > HeatCapacity;
      thermo_id.insert(exprFactory.register_expression( new HeatCapacity::Builder( thermoTags[0], tTag, yiTag )));
    }
    else if( thermo_quantity == "h" ){
      typedef Enthalpy < CellField > Enthalpy;
      thermo_id.insert(exprFactory.register_expression( new Enthalpy    ::Builder( thermoTags[0], tTag, yiTag )));
    }
  }
  else // species
  {
    if( thermo_quantity == "cp" ){
      typedef SpeciesHeatCapacity_Cp < CellField > SpeciesHeatCapacity;
      for( size_t n=0; n<nSpec; ++n ){
        thermo_id.insert(exprFactory.register_expression( new SpeciesHeatCapacity::Builder(thermoTags[n], tTag, n) ));
      }
    }
    else if( thermo_quantity == "h" ){
      typedef SpeciesEnthalpy   < CellField > SpeciesEnthalpy;
      for( size_t n=0; n<nSpec; ++n ){
        thermo_id.insert(exprFactory.register_expression( new SpeciesEnthalpy::Builder(thermoTags[n], tTag, n) ));
      }
    }
    else if( thermo_quantity == "cv" ){
      typedef SpeciesHeatCapacity_Cv < CellField > SpeciesHeatCapacity;
      for( size_t n=0; n<nSpec; ++n ){
        thermo_id.insert(exprFactory.register_expression( new SpeciesHeatCapacity::Builder(thermoTags[n], tTag, n) ));
      }
    }
  }
  return thermo_id;
}

//==============================================================================

void write_tree( const bool mix, const std::string& thermo_quantity, const Expr::ExpressionTree& tree)
{
  if( mix ){
    std::ofstream out( thermo_quantity + "Mixture.dot" );
    tree.write_tree(out);
  }
  else{
    std::ofstream out( thermo_quantity + "Species.dot" );
    tree.write_tree(out);
  }
}

//==============================================================================

std::vector< So::SpatFldPtr<CellField> >
get_cantera_results( const bool mix,
                     const std::string& thermo_quantity,
                     const bool timings,
                     Cantera_CXX::IdealGasMix& gasMix,
                     const int npts,
                     const int nSpec,
                     const CellField& prototype )
{
  using namespace SpatialOps;

  const double refPressure=gasMix.pressure();
  const std::vector<double>& molecularWeights = gasMix.molecularWeights();

  std::vector<double> tVec;
  for( size_t i=0; i<npts+2; ++i)
    tVec.push_back( 500.0 + 1000.0 * (i-0.5)/ npts);

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
  std::vector<double>::const_iterator itemp = tVec.begin();
  std::vector< std::vector<double> >::iterator imass = massfracs.begin();
  double eval_time;

  if( mix )
  {
    SpatFldPtr<CellField> canteraResult = SpatialFieldStore::get<CellField>(prototype);
    CellField::iterator icantend = canteraResult->end();
    boost::timer thermoTimer;

    if( thermo_quantity == "cp" ){
      for( CellField::iterator icant = canteraResult->begin(); icant!=icantend; ++itemp, ++imass, ++icant ){
        gasMix.setState_TPY( *itemp, refPressure, &(*imass)[0]);
        *icant=gasMix.cp_mass();
      }
      eval_time = thermoTimer.elapsed();
    }

    else if( thermo_quantity == "cv" ){
      for( CellField::iterator icant = canteraResult->begin(); icant!=icantend; ++itemp, ++imass, ++icant ){
        gasMix.setState_TPY( *itemp, refPressure, &(*imass)[0]);
        *icant=gasMix.cv_mass();
      }
      eval_time = thermoTimer.elapsed();
    }

    else if( thermo_quantity == "h" ){
      for( CellField::iterator icant = canteraResult->begin(); icant!=icantend; ++itemp, ++imass, ++icant ){
        gasMix.setState_TPY( *itemp, refPressure, &(*imass)[0]);
        *icant=gasMix.enthalpy_mass();
      }
      eval_time = thermoTimer.elapsed();
    }

    canteraResults.push_back( canteraResult );
  }

  else // species
  {
    for( size_t n=0; n < nSpec; ++n){
      canteraResults.push_back( SpatialFieldStore::get<CellField>(prototype) );
    }
    std::vector<double> thermo_result(nSpec,0.0);
    boost::timer thermoTimer;

    if( thermo_quantity == "cp" ){
      for( size_t i=0; i<npts+2; ++itemp, ++imass, ++i){
        gasMix.setState_TPY( *itemp, refPressure, &(*imass)[0]);
        gasMix.getPartialMolarCp(&thermo_result[0]);
        for( size_t n=0; n<nSpec; ++n){
          (*canteraResults[n])[i] = thermo_result[n];
        }
      }
      eval_time = thermoTimer.elapsed();
      for( size_t n=0; n<nSpec; ++n){
        *canteraResults[n] <<= *canteraResults[n] / molecularWeights[n]; // convert to mass basis for field comparison
      }
    }

    else if( thermo_quantity == "cv" ){
      for( size_t i=0; i<npts+2; ++itemp, ++imass, ++i){
        gasMix.setState_TPY( *itemp, refPressure, &(*imass)[0]);
        gasMix.getPartialMolarCp(&thermo_result[0]);
        for( size_t n=0; n<nSpec; ++n){
          (*canteraResults[n])[i] = thermo_result[n];
        }
      }
      eval_time = thermoTimer.elapsed();
      for( size_t n=0; n<nSpec; ++n){
        *canteraResults[n] <<= ( *canteraResults[n] - Cantera::GasConstant ) / molecularWeights[n]; // convert from molar cp to mass cv for field comparison
      }
    }

    else if( thermo_quantity == "h" ){
      for( size_t i=0; i<npts+2; ++itemp, ++imass, ++i){
        gasMix.setState_TPY( *itemp, refPressure, &(*imass)[0]);
        gasMix.getPartialMolarEnthalpies(&thermo_result[0]);
        for( size_t n=0; n<nSpec; ++n){
          (*canteraResults[n])[i] = thermo_result[n];
        }
      }
      eval_time = thermoTimer.elapsed();
      for( size_t n=0; n<nSpec; ++n){
        *canteraResults[n] <<= *canteraResults[n] / molecularWeights[n]; // convert to mass basis for field comparison
      }
    }
  }

  if( timings ) std::cout << "Cantera " + thermo_quantity + " time " << eval_time << std::endl;
  return canteraResults;
}

//==============================================================================

int main( int iarg, char* carg[] )
{
  try {
    std::string inputFileName;
    std::string inpGroup;
    bool mix     = false;
    bool timings = false;
    std::string thermo_quantity;

    // parse the command line options input describing the problem
    {
      po::options_description desc("Supported Options");
      desc.add_options()
           ( "help", "print help message" )
           ( "xml-input-file", po::value<std::string>(&inputFileName), "Cantera xml input file name" )
           ( "phase", po::value<std::string>(&inpGroup), "name of phase in Cantera xml input file" )
           ( "cv", "test isometric heat capacity" )
           ( "cp", "test isobaric heat capcity" )
           ( "h", "test enthalpy")
           ( "mix", "Triggers mixture heat capacity test.  Otherwise, species heat capacities are tested." )
           ( "timings", "Generate comparison timings between Cantera and PoKiTT across several problem sizes" );

      po::variables_map args;
      po::store( po::parse_command_line(iarg,carg,desc), args );
      po::notify(args);

      if( args.count("cv") ) thermo_quantity = "cv";
      else if( args.count("cp") ) thermo_quantity = "cp";
      else if( args.count("h") ) thermo_quantity = "h";
      else{
        std::cout << "You must enter the thermodynamic quantity you wish to calculate" << std::endl;
        std::cout << desc << std::endl;
        return 1;
      }

      mix = args.count("mix") > 0;
      timings = args.count("timings") > 0;

      if (!args.count("xml-input-file")){
        std::cout << "You must enter an xml input file for Cantera" << std::endl;
        std::cout << desc << std::endl;
        return 1;
      }

      if (args.count("help")) {
        std::cout << desc << "\n";
        return 1;
      }
    }

    TestHelper status( !timings ); // we don't need test helper output if we're running timings

    const CanteraObjects::Setup setup( "Mix", inputFileName, inpGroup );
    CanteraObjects::setup_cantera(setup);
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
      thermoTags.push_back(Expr::Tag( thermo_quantity + " mix", Expr::STATE_NONE));
    else{
      for( n=0; n<nSpec; ++n ){
        thermoTags.push_back( Expr::Tag( thermo_quantity + boost::lexical_cast<std::string>(n), Expr::STATE_NONE ) );
      }
    }

    Expr::ExpressionFactory exprFactory;

    exprFactory.register_expression( new Temp::Builder(tTag) );
    exprFactory.register_expression( new TemperaturePowers::Builder(tTag) );
    BOOST_FOREACH( Expr::Tag yiTag, yiTags){
      exprFactory.register_expression( new MassFracs::Builder (yiTag) );
    }
    std::set< Expr::ExpressionID > thermo_id = register_thermo_id( mix,
                                                                   thermo_quantity,
                                                                   exprFactory,
                                                                   thermoTags, tTag, yiTag,
                                                                   nSpec);

    std::vector<int> ptvec;
    if( timings ){
      ptvec.push_back(8*8*8);
      ptvec.push_back(16*16*16);
      ptvec.push_back(32*32*32);
      ptvec.push_back(64*64*64);
#     ifdef ENABLE_CUDA
      ptvec.push_back(128*128*128);
#     endif
    }
    else{
      ptvec.push_back(10);
    }

    for( std::vector<int>::iterator ptit = ptvec.begin(); ptit!= ptvec.end(); ++ptit){

      Expr::ExpressionTree thermoTree( thermo_id, exprFactory, 0 );
      write_tree( mix, thermo_quantity, thermoTree);

      So::IntVec npts(*ptit,1,1);
      const So::BoundaryCellInfo cellBCInfo = So::BoundaryCellInfo::build<CellField>(false,false,false);
      const So::GhostData cellGhosts(1);
      const So::MemoryWindow vwindow( So::get_window_with_ghost(npts,cellGhosts,cellBCInfo) );
      CellField xcoord( vwindow, cellBCInfo, cellGhosts, NULL );

      std::vector<double> length(3,1.0);
      So::Grid grid( npts, length );
      grid.set_coord<SpatialOps::XDIR>( xcoord );
#     ifdef ENABLE_CUDA
      xcoord.add_device( GPU_INDEX );
#     endif

      Expr::FieldManagerList fml;

      thermoTree.register_fields( fml );
      fml.allocate_fields( Expr::FieldAllocInfo( npts, 0, 0, false, false, false ) );
      thermoTree.bind_fields( fml );

      using namespace SpatialOps;
      Expr::FieldMgrSelector<CellField>::type& cellFM = fml.field_manager< CellField>();

      CellField& temp = cellFM.field_ref(tTag);
      temp <<= 500.0 + 1000.0*xcoord;

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

      if( timings ) std::cout << std::endl << setup.inputFile << " - " << *ptit << std::endl;

      boost::timer thermoTimer;
      thermoTree.execute_tree();

      if( timings ) std::cout << "PoKiTT " + thermo_quantity + " time  " << thermoTimer.elapsed() << std::endl;

#     ifdef ENABLE_CUDA
      BOOST_FOREACH( Expr::Tag thermoTag, thermoTags){
        CellField& thermo = fml.field_manager<CellField>().field_ref(thermoTag);
        thermo.add_device(CPU_INDEX);
      }
#     endif

      const std::vector< SpatFldPtr<CellField> > canteraResults = get_cantera_results( mix, thermo_quantity, timings, *gasMix, *ptit, nSpec, xcoord );

      std::vector< SpatFldPtr<CellField> >::const_iterator icantera = canteraResults.begin();
      BOOST_FOREACH( const Expr::Tag& thermoTag, thermoTags ){
        CellField& thermo = cellFM.field_ref(thermoTag);
        status( field_equal( thermo, **icantera, 1e-14 ), thermoTag.name() );
        ++icantera;
      }

    } // number of points

    if( status.ok() ){
      std::cout << "PASS\n";
      return 0;
    }

  }
  catch( Cantera::CanteraError& ){
    Cantera::showErrors();
  }
  catch( std::exception& err ){
    std::cout << err.what() << std::endl;
  }

  std::cout << "FAIL\n";
  return -1;
}

//==============================================================================
