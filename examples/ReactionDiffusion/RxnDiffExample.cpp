/*
 * RxnDiffExample.cpp
 *
 *  Created on: Feb 4, 2014
 *      Author: Nathan Yonkee
 */

#include <iostream>

#include <expression/ExprLib.h>

#include <spatialops/structured/Grid.h>
#include <spatialops/OperatorDatabase.h>
#include <spatialops/structured/FieldComparisons.h>
#include <spatialops/util/TimeLogger.h>

#include <pokitt/CanteraObjects.h>
#include "test/TestHelper.h"
#include "SpeciesRHS.h"
#include "MassFractions.h"
#include <pokitt/thermo/Temperature.h>
#include <pokitt/thermo/Pressure.h>
#include <pokitt/thermo/Density.h>
#include <pokitt/MixtureMolWeight.h>
#include <pokitt/kinetics/ReactionRates.h>
#include <pokitt/thermo/Enthalpy.h>

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>

namespace Cantera_CXX{ class IdealGasMix; }

namespace SO = SpatialOps;
typedef SO::SVolField  CellField;
typedef SO::SpatFldPtr<CellField> CellFieldPtrT;
typedef SO::BasicOpTypes<CellField>::GradX GradOp;
typedef SO::BasicOpTypes<CellField>::DivX  DivOp;

namespace po = boost::program_options;

using namespace pokitt;

bool driver( const bool timings,
             const size_t nSteps,
             const double dt )
{
  TestHelper status( !timings );
  Cantera_CXX::IdealGasMix* const gasMix = CanteraObjects::get_gasmix();
  const int nSpec=gasMix->nSpecies();

  const Expr::Tag fluxTag ( "for testing",Expr::STATE_NONE );
  const Expr::Tag xTag  ( "XCoord",      Expr::STATE_NONE );
  const Expr::Tag tTag  ( "Temp",      Expr::STATE_NONE );
  const Expr::Tag pTag  ( "Pressure",    Expr::STATE_NONE);
  const Expr::Tag rhoTag  ( "Density",    Expr::STATE_NONE);
  const Expr::Tag mmwTag( "mmw",         Expr::STATE_NONE);
  const Expr::Tag hTag    ( "enthalpy",           Expr::STATE_NONE);
  Expr::TagList yiTags;
  Expr::TagList rhoYiTags;
  Expr::TagList specRHSTags;
  Expr::TagList rTags;
  for( size_t n=0; n<nSpec; ++n ){
    yiTags.push_back( Expr::Tag( "yi_" + boost::lexical_cast<std::string>(n), Expr::STATE_NONE ) );

    rTags.push_back(  Expr::Tag( "ri_" + boost::lexical_cast<std::string>(n), Expr::STATE_NONE ) );
  }
  for( size_t n=0; n<nSpec-1; ++n ){
    specRHSTags.push_back( Expr::Tag( "specRHS_" + boost::lexical_cast<std::string>(n), Expr::STATE_N ) );
    rhoYiTags.push_back( Expr::Tag( "rhoYi_" + boost::lexical_cast<std::string>(n), Expr::STATE_N ) );
  }

  const double rho0 = 1;
  const double t0 = 1200;

  Expr::ExpressionFactory initFactory;
  std::set<Expr::ExpressionID> initIDs;
  Expr::ExpressionID hID; // temperature
  Expr::ExpressionID rhoID; //density
  {
    typedef Expr::ConstantExpr  <CellField>::Builder Density;
    typedef Expr::LinearFunction<CellField>::Builder rhoYi;
    typedef Expr::ConstantExpr  <CellField>::Builder Temperature;
    typedef Enthalpy            <CellField>::Builder Enthalpy;
    typedef MassFractions<CellField>::Builder MassFracs;
    initFactory.register_expression( new Temperature( tTag, t0 )    );
    initFactory.register_expression( new Density( rhoTag, rho0 )    );
    initFactory.register_expression( new MassFracs( yiTags, rhoTag, rhoYiTags )    );
    for( size_t n = 0; n < (nSpec-1); ++n ){
      if( n != 0 && n != 3 )
        initFactory.register_expression( new rhoYi ( rhoYiTags[n], rhoTag, 1e-12, 0.0  ) );
    }
    initFactory.register_expression( new rhoYi ( rhoYiTags[0], rhoTag, 0.01, 0.0  ) );
    initFactory.register_expression( new rhoYi ( rhoYiTags[3], rhoTag, 0.20, 0.0  ) );
    hID = initFactory.register_expression( new Enthalpy  ( hTag , tTag, yiTags )    );
    initIDs.insert( hID );
  }
  Expr::ExpressionTree initTree( initIDs, initFactory, 0 );

  Expr::ExpressionFactory execFactory;
  {
    typedef Expr::PlaceHolder     <CellField>::Builder XCoord;
    typedef Temperature <CellField>::Builder Temperature;
//    typedef Expr::ConstantExpr<CellField>::Builder Temperature;
    typedef MassFractions<CellField>::Builder MassFracs;
    typedef MixtureMolWeight<CellField>::Builder MixMolWeight;
    typedef       ReactionRates <CellField>::Builder ReactionRates;
    typedef Expr::PlaceHolder<CellField>::Builder Density;
    typedef Pressure<CellField>::Builder Pressure;
    typedef Expr::PlaceHolder<CellField>::Builder Enthalpy;
    typedef Expr::ConstantExpr<DivOp::SrcFieldType>::Builder Flux;
    typedef SpeciesRHS<DivOp>::Builder SpeciesRHS;
    execFactory.register_expression( new Density( rhoTag ) );
    execFactory.register_expression( new Enthalpy( hTag ) );
    execFactory.register_expression( new MassFracs( yiTags, rhoTag, rhoYiTags));
    execFactory.register_expression( new MixMolWeight ( mmwTag, yiTags )           );
    execFactory.register_expression( new ReactionRates(rTags, tTag, pTag, yiTags, mmwTag) );
    execFactory.register_expression( new Pressure ( pTag, tTag, rhoTag, mmwTag) );
    execFactory.register_expression( new Temperature ( tTag, yiTags, hTag ) );;
    execFactory.register_expression( new Flux( fluxTag, 0 ) );
    for( size_t n = 0; n < (nSpec-1); ++n ){
      execFactory.register_expression( new SpeciesRHS ( specRHSTags[n], fluxTag, rTags[n]  ) );
    }
  }

  std::vector<SO::IntVec> sizeVec;
  if( timings ){
    sizeVec.push_back( SO::IntVec(126,126,126) );
    sizeVec.push_back( SO::IntVec( 62, 62, 62) );
    sizeVec.push_back( SO::IntVec( 30, 30, 30) );
    sizeVec.push_back( SO::IntVec( 14, 14, 14) );
    sizeVec.push_back( SO::IntVec(  6,  6,  6) );
  }
  else{
    sizeVec.push_back( SO::IntVec( 20,  1,  1) );
  }

  for( std::vector<SO::IntVec>::iterator iSize = sizeVec.begin(); iSize!= sizeVec.end(); ++iSize){

    SO::IntVec gridSize = *iSize;
    Expr::ExprPatch patch( gridSize[0], gridSize[1], gridSize[2] );

    Expr::FieldManagerList& fml = patch.field_manager_list();
    initTree.register_fields( fml );
    initTree.bind_fields( fml );
    initTree.lock_fields( fml );

    fml.allocate_fields( patch.field_info() );
    Expr::TimeStepper timeIntegrator( execFactory, Expr::FORWARD_EULER, patch.id() );

    for( size_t n = 0; n < (nSpec-1); ++n ){
      timeIntegrator.add_equation<CellField>( (rhoYiTags[n].name()), specRHSTags[n], SpatialOps::GhostData(0) );
    }
    SO::Grid grid( gridSize, SO::DoubleVec(1,1,1) );
    SpatialOps::OperatorDatabase& opDB = patch.operator_database();
    SpatialOps::build_stencils( grid, opDB );

    timeIntegrator.finalize( fml, patch.operator_database(), patch.field_info() );
    {
      std::ofstream out( "rxn_example.dot" );
      timeIntegrator.get_tree()->write_tree( out );
    }
    timeIntegrator.get_tree()->lock_fields(fml);

    initTree.execute_tree();

    Timer timer;
    timer.start();
    for( size_t s = 0; s < nSteps; ++s ){
      timeIntegrator.step( dt );
    }
    timer.stop();

    CellField& t = fml.field_ref< CellField >( tTag );
    const double tMean = SO::nebo_sum_interior( t ) / ( gridSize[0] * gridSize[1] * gridSize[2] );
    status( tMean >= (t0 - 100) && tMean < 5000 );

    fml.deallocate_fields();
  } // number of points

  return status.ok();
}

//==============================================================================

int main( int iarg, char* carg[] )
{
  std::string inputFileName;
  std::string inpGroup;
  bool timings = false;
  size_t nSteps = 1000;
  double dt = 1e-12;

  // parse the command line options input describing the problem
  try {
    po::options_description desc("Supported Options");
    desc.add_options()
           ( "help", "print help message" )
           ( "xml-input-file", po::value<std::string>(&inputFileName), "Cantera xml input file name" )
           ( "phase", po::value<std::string>(&inpGroup), "name of phase in Cantera xml input file" )
           ( "timings", "Generate comparison timings between Cantera and PoKiTT across several problem sizes" )
           ( "nsteps", po::value<size_t>(&nSteps), "How many time steps to take" )
           ( "dt",     po::value<double>(&dt),     "Size of time steps to take"  );

    po::variables_map args;
    po::store( po::parse_command_line(iarg,carg,desc), args );
    po::notify(args);

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
    status( driver( timings, nSteps, dt ), "Reaction Diffusion" );

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
