/*
 * RxnDiffExample.cpp
 *
 *  Created on: Feb 4, 2014
 *      Author: Nathan Yonkee
 */

#include <iostream>

#include <expression/ExprLib.h>
#include <expression/BoundaryConditionExpression.h>

#include <spatialops/structured/Grid.h>
#include <spatialops/OperatorDatabase.h>
#include <spatialops/structured/FieldComparisons.h>
#include <spatialops/util/TimeLogger.h>

#include <pokitt/CanteraObjects.h>
#include "test/TestHelper.h"
#include "EnthalpyTransport.h"
#include "SpeciesTransport.h"

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>

namespace Cantera_CXX{ class IdealGasMix; }

namespace SO = SpatialOps;
typedef SO::SVolField  CellField;
typedef SO::FaceTypes<CellField> FaceTypes;
typedef FaceTypes::XFace XFluxT;
typedef FaceTypes::YFace YFluxT;
using SO::IntVec;
using Expr::Tag;
using Expr::TagList;

namespace po = boost::program_options;

using namespace pokitt;

void
initialize_mask_points( const SO::Grid& grid,
                        std::vector<IntVec>& leftSet,
                        std::vector<IntVec>& rightSet,
                        std::vector<IntVec>& topSet,
                        std::vector<IntVec>& bottomSet )
{
  for( int i=0; i < grid.extent(1); ++i ){
    leftSet .push_back( IntVec(0, i, 0) );
    rightSet.push_back( IntVec(grid.extent(0), i, 0) );
  }

  for(int i = 0; i < grid.extent(0); i++){
    topSet   .push_back( IntVec(i, 0, 0) );
    bottomSet.push_back( IntVec(i, grid.extent(1), 0) );
  }
}

//==============================================================================

void setup_bcs( Expr::ExpressionFactory& execFactory,
                const SO::Grid& grid,
                const SO::GhostData ghosts,
                const Tag& fieldTag,
                const double value,
                const SO::BCType type)
{
  std::vector<IntVec> xminusPts, xplusPts, yminusPts, yplusPts;
  initialize_mask_points( grid, xminusPts, xplusPts, yminusPts, yplusPts );

  const SO::BoundaryCellInfo xInfo = SO::BoundaryCellInfo::build<XFluxT>(true,true,true);
  const SO::BoundaryCellInfo yInfo = SO::BoundaryCellInfo::build<YFluxT>(true,true,true);
  const SO::MemoryWindow xwindow( SO::get_window_with_ghost( grid.extent(), ghosts, xInfo ) );
  const SO::MemoryWindow ywindow( SO::get_window_with_ghost( grid.extent(), ghosts, yInfo ) );

  typedef Expr::ConstantBCOpExpression<CellField,SO::XDIR>::Builder XBC;
  typedef Expr::ConstantBCOpExpression<CellField,SO::YDIR>::Builder YBC;

  XBC::MaskPtr xminus( new XBC::MaskType( xwindow, xInfo, ghosts, xminusPts ) );
  XBC::MaskPtr xplus ( new XBC::MaskType( xwindow, xInfo, ghosts, xplusPts  ) );
  YBC::MaskPtr yminus( new YBC::MaskType( ywindow, yInfo, ghosts, yminusPts ) );
  YBC::MaskPtr yplus ( new YBC::MaskType( ywindow, yInfo, ghosts, yplusPts  ) );

# ifdef ENABLE_CUDA
  // Masks are created on CPU so we need to explicitly transfer them to GPU
  xminus->add_consumer( GPU_INDEX );
  xplus ->add_consumer( GPU_INDEX );
  yminus->add_consumer( GPU_INDEX );
  yplus ->add_consumer( GPU_INDEX );
# endif

  Tag xmbcTag( fieldTag.name() + "xmbc", Expr::STATE_NONE );
  Tag xpbcTag( fieldTag.name() + "xpbc", Expr::STATE_NONE );
  Tag ymbcTag( fieldTag.name() + "ymbc", Expr::STATE_NONE );
  Tag ypbcTag( fieldTag.name() + "ypbc", Expr::STATE_NONE );

  execFactory.register_expression( new XBC( xmbcTag, xminus, type, SO::MINUS_SIDE,  value ) );
  execFactory.register_expression( new XBC( xpbcTag, xplus , type, SO::PLUS_SIDE,  value ) );
  execFactory.register_expression( new YBC( ymbcTag, yminus, type, SO::MINUS_SIDE,  value ) );
  execFactory.register_expression( new YBC( ypbcTag, yplus , type, SO::PLUS_SIDE,  value ) );

  execFactory.attach_modifier_expression( xmbcTag, fieldTag );
  execFactory.attach_modifier_expression( xpbcTag, fieldTag );
  execFactory.attach_modifier_expression( ymbcTag, fieldTag );
  execFactory.attach_modifier_expression( ypbcTag, fieldTag );
}

bool driver( const bool timings,
             const size_t nSteps,
             const double dt )
{
  TestHelper status( !timings );
  Cantera_CXX::IdealGasMix* const gasMix = CanteraObjects::get_gasmix();
  const int nSpec=gasMix->nSpecies();

  const double rho0 = 0.3;
  const double length = 1e-3;
  const double tMax = 1700;
  const double tDev = 0.2*length;
  const double tMean = length/2;
  const double tBase = 400;
  std::vector<double> yi(nSpec, 0.0);
  yi[0]=0.01;
  yi[3]=0.20;
  std::cout<<"here00\n";
  const Tag xTag    ( "XCoord",         Expr::STATE_NONE );
  const Tag yTag    ( "YCoord",         Expr::STATE_NONE );
  const Tag rhoTag  ( "Density",        Expr::STATE_NONE );
  const Tag mmwTag  ( "mix MW ",        Expr::STATE_NONE );

  const Tag tTag    ( "Temp",           Expr::STATE_NONE );
  const Tag hTag    ( "enthalpy",       Expr::STATE_N    );

  const Tag jTag    ( "J",              Expr::STATE_NONE );
  const Tag rhoYiTag( "rhoYi",          Expr::STATE_N    );
  TagList yiTags;
  for( size_t n=0; n<nSpec; ++n ){
    std::string spec = boost::lexical_cast<std::string>(n);
    yiTags.push_back( Tag("yi"+spec,    Expr::STATE_NONE ) );
  }
  std::cout<<"here01\n";
  Expr::ExpressionFactory initFactory;
  Expr::ExpressionFactory execFactory;

  typedef EnthalpyTransport<CellField> EnthalpyTransport;
  typedef SpeciesTransport <CellField> SpeciesTransport;
  std::cout<<"here02\n";
  std::set<Expr::ExpressionID> initIDs;
  std::cout<<"here03\n";
  std::list<SpeciesTransport*> specEqns;
  for( size_t n=0; n<(nSpec-1); ++n ){
    std::string spec = boost::lexical_cast<std::string>(n);
    SpeciesTransport* specTrans = new SpeciesTransport( execFactory, rhoYiTag, jTag, spec );
    initIDs.insert( specTrans->initial_condition( initFactory, yiTags[n], yi[n], rho0 ) );
    specEqns.push_back( specTrans );
  }
  specEqns.front()->register_one_time_expressions(execFactory, rhoTag, yiTags, mmwTag, tTag );
  std::cout<<"here04\n";
  EnthalpyTransport* hTrans = new EnthalpyTransport( execFactory, hTag, tTag, rhoTag, yiTags, jTag, mmwTag );
  std::cout<<"here041\n";
  initIDs.insert( hTrans->initial_condition( initFactory, xTag, yTag, tMax, tDev, tMean, tBase ) );
std::cout<<"here05\n";
  Expr::ExpressionTree initTree( initIDs, initFactory, 0 );

  std::vector< IntVec > sizeVec;
  if( timings ){
    sizeVec.push_back( IntVec(126,126, 1 ) );
    sizeVec.push_back( IntVec( 62, 62, 1 ) );
    sizeVec.push_back( IntVec( 30, 30, 1 ) );
    sizeVec.push_back( IntVec( 14, 14, 1 ) );
    sizeVec.push_back( IntVec(  6,  6, 1 ) );
  }
  else{
    sizeVec.push_back( IntVec( 6,  6,  1) );
  }
  std::cout<<"here1\n";
  for( std::vector< IntVec >::iterator iSize = sizeVec.begin(); iSize!= sizeVec.end(); ++iSize){

    IntVec gridSize = *iSize;
    Expr::ExprPatch patch( gridSize[0], gridSize[1], gridSize[2] );
    std::cout<<"here2\n";
    Expr::FieldManagerList& fml = patch.field_manager_list();
    initTree.register_fields( fml );
    initTree.bind_fields( fml );
    initTree.lock_fields( fml );
    std::cout<<"here3\n";
    fml.allocate_fields( patch.field_info() );
    Expr::TimeStepper timeIntegrator( execFactory, Expr::FORWARD_EULER, patch.id() );
    SO::GhostData ghosts( IntVec( 1, 1, 0), IntVec( 1, 1, 0 ) ); // 1 on +-x and +- y and 0 on z

    SO::Grid grid( gridSize, SO::DoubleVec( length, length, 1 ) );
    SO::OperatorDatabase& opDB = patch.operator_database();
    SO::build_stencils( grid, opDB );
    std::cout<<"here4\n";
    for( std::list<SpeciesTransport*>::iterator iEqn=specEqns.begin(); iEqn!=specEqns.end(); ++iEqn ){
      timeIntegrator.add_equation<CellField>( (*iEqn)->solution_variable_name(), (*iEqn)->get_rhs_tag(), ghosts );
      (*iEqn)->setup_boundary_conditions( grid, execFactory );
    }
    timeIntegrator.add_equation<CellField>( hTrans->solution_variable_name(), hTrans->get_rhs_tag(), ghosts );
    hTrans->setup_boundary_conditions( grid, execFactory );
    std::cout<<"here5\n";
    timeIntegrator.finalize( fml, patch.operator_database(), patch.field_info() );
    {
      std::ofstream out( "rxn_example.dot" );
      timeIntegrator.get_tree()->write_tree( out );
    }
    timeIntegrator.get_tree()->lock_fields(fml);

    CellField& xcoord = fml.field_ref< CellField >( xTag );
    CellField& ycoord = fml.field_ref< CellField >( yTag );
    grid.set_coord<SO::XDIR>( xcoord );
    grid.set_coord<SO::YDIR>( ycoord );
#   ifdef ENABLE_CUDA
    xcoord.set_device_as_active( GPU_INDEX );
    ycoord.set_device_as_active( GPU_INDEX );
#   endif

    initTree.execute_tree();

    Timer timer;
    timer.start();
    for( size_t s = 0; s < nSteps; ++s ){
      timeIntegrator.step( dt );
    }
    timer.stop();

    CellField& t = fml.field_ref< CellField >( tTag );
    const double tMean = SO::nebo_sum_interior( t ) / ( gridSize[0] * gridSize[1] * gridSize[2] );
    status( tMean >= 200 && tMean <= 3500 ); // bounds on NASA polynomials

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
  double dt = 1e-11;

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
