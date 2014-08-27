/*
 * HeatCapacity_Cp_test.cpp
 *
 *  Created on: Mar 31, 2014
 *      Author: Nathan Yonkee
 */

#include <iostream>
#include <stdio.h>
#include <fstream>
#include "TestHelper.h"

#include <pokitt/thermo/TemperaturePowers.h>
#include <pokitt/thermo/HeatCapacity_Cp.h>

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

// helper functions which differentiate between mixture and species test cases
std::set< Expr::ExpressionID > register_cp( bool mix, Expr::ExpressionFactory& exprFactory, Expr::TagList cpTags, Expr::Tag tTag, Expr::Tag yiTag, int nSpec);
void write_tree( bool mix, Expr::ExpressionTree& tree);
std::vector< So::SpatFldPtr<CellField> > get_Cantera_results( bool mix, bool timings, Cantera_CXX::IdealGasMix* gasMix, int npts, int nSpec, CellField prototype);

int main( int iarg, char* carg[] )
{
  try {
    std::string inputFileName;
    std::string inpGroup;
    bool mix;
    bool timings;
    // parse the command line options input describing the problem
    {
      po::options_description desc("Supported Options");
      desc.add_options()
                ( "help", "print help message" )
                ( "input", po::value<std::string>(&inputFileName)->default_value("thermo_tester.xml"), "Cantera input file name" )
                ( "phase", po::value<std::string>(&inpGroup)->default_value("const_cp"), "name of phase in Cantera input" )
                ( "mix", po::value<bool>(&mix)->default_value(true), "false for species heat capacity" )
                ( "time", po::value<bool>(&timings)->default_value(false), "true to run with timings" );

      po::variables_map args;
      po::store( po::parse_command_line(iarg,carg,desc), args );
      po::notify(args);

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
    Expr::TagList cpTags;
    if( mix )
      cpTags.push_back(Expr::Tag("cp mix", Expr::STATE_NONE));
    else{
      for( n=0; n<nSpec; ++n ){
        cpTags.push_back( Expr::Tag( "cp" + boost::lexical_cast<std::string>(n), Expr::STATE_NONE ) );
      }
    }

    Expr::ExpressionFactory exprFactory;

    exprFactory.register_expression( new Temp::Builder(tTag) );
    exprFactory.register_expression( new TemperaturePowers::Builder(tTag) );
    BOOST_FOREACH( Expr::Tag yiTag, yiTags){
      exprFactory.register_expression( new MassFracs::Builder (yiTag) );
    }
    std::set< Expr::ExpressionID > cp_id = register_cp( mix, exprFactory, cpTags, tTag, yiTag, nSpec);

    std::vector<int> ptvec;
    if(timings){
      ptvec.push_back(8*8*8);
      ptvec.push_back(16*16*16);
      ptvec.push_back(32*32*32);
      ptvec.push_back(64*64*64);
#     ifdef ENABLE_CUDA
      ptvec.push_back(128*128*128);
#     endif
    }
    else
      ptvec.push_back(10);

    for( std::vector<int>::iterator ptit = ptvec.begin(); ptit!= ptvec.end(); ++ptit){

      Expr::ExpressionTree cpTree( cp_id, exprFactory, 0 );
      write_tree( mix, cpTree);

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

      cpTree.register_fields( fml );
      fml.allocate_fields( Expr::FieldAllocInfo( npts, 0, 0, false, false, false ) );
      cpTree.bind_fields( fml );

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

      cpTree.lock_fields( fml );  // prevent fields from being deallocated so that we can get them after graph execution.

      if(timings)
        std::cout << std::endl << setup.inputFile << " - " << *ptit << std::endl;

      boost::timer cpTimer;
      cpTree.execute_tree();

      if(timings)
        std::cout << "PoKiTT cp time  " << cpTimer.elapsed() << std::endl;

#     ifdef ENABLE_CUDA
      BOOST_FOREACH( Expr::Tag cpTag, cpTags){
        CellField& cp = fml.field_manager<CellField>().field_ref(cpTag);
        cp.add_device(CPU_INDEX);
      }
#     endif

      std::vector< SpatFldPtr<CellField> > canteraResults = get_Cantera_results( mix, timings, gasMix, *ptit, nSpec, xcoord);

      std::vector< SpatFldPtr<CellField> >::iterator icantera = canteraResults.begin();
      BOOST_FOREACH( Expr::Tag cpTag, cpTags){
        CellField& cp = cellFM.field_ref(cpTag);
        status(field_equal(cp, **icantera, 1e-14), cpTag.name());
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

std::set< Expr::ExpressionID > register_cp( bool mix, Expr::ExpressionFactory& exprFactory, Expr::TagList cpTags, Expr::Tag tTag, Expr::Tag yiTag, int nSpec)
{
  typedef HeatCapacity_Cp        < CellField > HeatCapacity;
  typedef SpeciesHeatCapacity_Cp < CellField > SpeciesHeatCapacity;

  std::set< Expr::ExpressionID > cp_id;
  if( mix )
    cp_id.insert(exprFactory.register_expression( new HeatCapacity::Builder( cpTags[0], tTag, yiTag )));
  else{
    for( size_t n=0; n<nSpec; ++n ){
      cp_id.insert(exprFactory.register_expression( new SpeciesHeatCapacity::Builder(cpTags[n], tTag, n) ));
    }
  }
  return cp_id;
}

void write_tree( bool mix, Expr::ExpressionTree& tree)
{
  if(mix){
    std::ofstream out( "CpMixture.dot" );
    tree.write_tree(out);
  }
  else{
    std::ofstream out( "CpSpecies.dot" );
    tree.write_tree(out);
  }
}

std::vector< So::SpatFldPtr<CellField> > get_Cantera_results( bool mix, bool timings, Cantera_CXX::IdealGasMix* gasMix, int npts, int nSpec, CellField prototype)
{
  using namespace SpatialOps;
  size_t i;
  size_t n;

  const double refPressure=gasMix->pressure();
  const std::vector<double>& molecularWeights = gasMix->molecularWeights();

  std::vector<double> tVec;
  for( i=0; i<npts+2; ++i)
    tVec.push_back( 500.0 + 1000.0 * (i-0.5)/ npts);

  std::vector< std::vector<double> > massfracs;
  for( i=0; i<npts+2; ++i){
    std::vector<double> massfrac;
    double sum = 0.0;
    for( n=0; n<nSpec; ++n){
      massfrac.push_back(1 + n + (i-0.5)/ npts);
      sum+=massfrac[n];
    }
    for( n=0; n<nSpec; ++n)
      massfrac[n] = massfrac[n]/sum;
    massfracs.push_back(massfrac);
  }

  std::vector< SpatFldPtr<CellField> > canteraResults;
  std::vector<double>::const_iterator itemp = tVec.begin();
  std::vector< std::vector<double> >::iterator imass = massfracs.begin();

  if( mix ){
    SpatFldPtr<CellField> canteraResult = SpatialFieldStore::get<CellField>(prototype);
    boost::timer cpTimer;
    for(CellField::iterator icant = canteraResult->begin(); icant!=canteraResult->end(); ++itemp, ++imass, ++icant){
      gasMix->setState_TPY( *itemp, refPressure, &(*imass)[0]);
      *icant=gasMix->cp_mass();
    }
    if(timings)
      std::cout << "Cantera cp time " << cpTimer.elapsed() << std::endl;
    canteraResults.push_back(canteraResult);
  }
  else{
    for( n=0; n < nSpec; ++n){
      canteraResults.push_back(SpatialFieldStore::get<CellField>(prototype));
    }
    std::vector<double> cp_result(nSpec,0.0);
    boost::timer cpTimer;
    for( i=0; i<npts+2; ++itemp, ++imass, ++i){
      gasMix->setState_TPY( *itemp, refPressure, &(*imass)[0]);
      gasMix->getPartialMolarCp(&cp_result[0]);
      for( n=0; n<nSpec; ++n){
        (*canteraResults[n])[i] = cp_result[n];
      }
    }
    if( timings)
      std::cout << "Cantera cp time " << cpTimer.elapsed() << std::endl;
    for( n=0; n<nSpec; ++n){
      *canteraResults[n] <<= *canteraResults[n] / molecularWeights[n]; // convert to mass basis for field comparison
    }
  }
  return canteraResults;
}
