/*
 * The MIT License
 *
 * Copyright (c) 2015 The University of Utah
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
#include "EnthalpyTransport.h"
#include "SpeciesTransport.h"
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>
#include "TagManager.h"

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

void print_fields( Expr::FieldManagerList& fml, const TagList& fieldTags ){
  BOOST_FOREACH( const Tag iTag, fieldTags){
    CellField& field = fml.field_ref< CellField >( iTag );
    std::cout << iTag.name() << std::endl;
    SO::print_field( field, std::cout );
  }
}

bool driver( const bool timings,
             const bool print,
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
  const double tBase = 800;
  std::vector<double> yi(nSpec, 0.0);
  yi[0]=0.01;
  yi[3]=0.15;

  const TagManager tagMgr( nSpec,           // enum values:
      Tag( "rhoYi",     Expr::STATE_N    ), // RHOYI
      Tag( "enthalpy",  Expr::STATE_N    ), // H
      Tag( "yi",        Expr::STATE_NONE ), // YI
      Tag( "temp",      Expr::STATE_NONE ), // T
      Tag( "pressure",  Expr::STATE_NONE ), // P
      Tag( "density",   Expr::STATE_NONE ), // RHO
      Tag( "mix MW ",   Expr::STATE_NONE ), // MMW
      Tag( "r",         Expr::STATE_NONE ), // R
      Tag( "lambda",    Expr::STATE_NONE ), // LAM
      Tag( "D",         Expr::STATE_NONE ), // D
      Tag( "J",         Expr::STATE_NONE ), // J
      Tag( "heat flux", Expr::STATE_NONE ), // Q
      Tag( "XCoord",    Expr::STATE_NONE ), // XCOORD
      Tag( "YCoord",    Expr::STATE_NONE ) // YCOORD
  );

  std::vector< IntVec > sizeVec;
  if( timings ){
    sizeVec.push_back( IntVec(1022, 1022, 1 ) );
    sizeVec.push_back( IntVec(510,  510,  1 ) );
    sizeVec.push_back( IntVec(254,  254,  1 ) );
    sizeVec.push_back( IntVec(126,  126,  1 ) );
    sizeVec.push_back( IntVec( 62,  62,   1 ) );
    sizeVec.push_back( IntVec( 30,  30,   1 ) );
  }
  sizeVec.push_back( IntVec( 6,  6,  1 ) );
  for( std::vector< IntVec >::iterator iSize = sizeVec.begin(); iSize!= sizeVec.end(); ++iSize){

    typedef EnthalpyTransport<CellField> EnthalpyTransport;
    typedef SpeciesTransport <CellField> SpeciesTransport;

    Expr::ExpressionFactory initFactory;
    Expr::ExpressionFactory execFactory;

    std::set<Expr::ExpressionID> initIDs;
    EnthalpyTransport* hTrans = new EnthalpyTransport( execFactory, tagMgr );
    initIDs.insert( hTrans->initial_condition( initFactory, tMax, tBase, tMean, tDev) );

    std::list<SpeciesTransport*> specEqns;
    for( size_t n=0; n<(nSpec-1); ++n ){
      SpeciesTransport* specTrans = new SpeciesTransport( execFactory, tagMgr, n );
      initIDs.insert( specTrans->initial_condition( initFactory, yi[n], rho0 ) );
      specEqns.push_back( specTrans );
    }
    specEqns.front()->register_one_time_expressions( execFactory );

    Expr::ExpressionTree initTree( initIDs, initFactory, 0 );

    IntVec gridSize = *iSize;
    Expr::ExprPatch patch( gridSize[0], gridSize[1], gridSize[2] );
    SO::Grid grid( gridSize, SO::DoubleVec( length, length, 1 ) );
    SO::GhostData ghosts( IntVec( 1, 1, 0), IntVec( 1, 1, 0 ) ); // 1 on +-x and +- y and 0 on z

    Expr::TimeStepper timeIntegrator( execFactory, Expr::FORWARD_EULER, patch.id() );
    for( std::list<SpeciesTransport*>::iterator iEqn=specEqns.begin(); iEqn!=specEqns.end(); ++iEqn ){
      timeIntegrator.add_equation<CellField>( (*iEqn)->solution_variable_name(), (*iEqn)->get_rhs_tag(), ghosts );
      (*iEqn)->setup_boundary_conditions( grid, execFactory );
    }
    timeIntegrator.add_equation<CellField>( hTrans->solution_variable_name(), hTrans->get_rhs_tag(), ghosts );
    hTrans->setup_boundary_conditions( grid, execFactory );

    Expr::FieldManagerList& fml = patch.field_manager_list();
#   ifdef ENABLE_CUDA
    initTree.set_device_index( GPU_INDEX, fml );
#   endif
    fml.allocate_fields( patch.field_info() );
    SO::OperatorDatabase& opDB = patch.operator_database();
    SO::build_stencils( grid, opDB );

    timeIntegrator.finalize( fml, patch.operator_database(), patch.field_info() );
    {
      std::ofstream out( "rxn_example.dot" );
      timeIntegrator.get_tree()->write_tree( out );
    }
    timeIntegrator.get_tree()->lock_fields(fml);

    initTree.register_fields( fml );
    initTree.bind_fields( fml );
    initTree.lock_fields( fml );

    CellField& xcoord = fml.field_ref< CellField >( tagMgr[XCOORD] );
    CellField& ycoord = fml.field_ref< CellField >( tagMgr[YCOORD] );
    grid.set_coord<SO::XDIR>( xcoord );
    grid.set_coord<SO::YDIR>( ycoord );

    initTree.execute_tree();

    if( timings ) std::cout << "\nPoKiTT Reaction Diffusion size " << xcoord.window_with_ghost().glob_npts() << std::endl;
    SpatialOps::Timer timer;
    timer.start();
    for( size_t s = 0; s <= nSteps; ++s ){
      if( s%5000 == 0 && print){
        timer.stop();
        std::cout<<"Fields at time "<< s*dt << "; step " << s << "; simulation run time " << timer.elapsed_time() << std::endl;
        print_fields( fml, tag_list( tagMgr[T], tagMgr.rN( 0 ), tagMgr.rhoYiN( 0 ) ) );
        timer.start();
      }
      timeIntegrator.step( dt );
    }
    timer.stop();
    if( timings ) std::cout << "PoKiTT Reaction Diffusion time per step " << timer.elapsed_time() / (nSteps + 1) << std::endl;

    CellField& t = fml.field_ref< CellField >( tagMgr[T] );
#   ifdef ENABLE_CUDA
    t.set_device_as_active( CPU_INDEX );
#   endif
    const double tMean = SO::nebo_sum_interior( t ) / ( gridSize[0] * gridSize[1] * gridSize[2] );
    status( tMean >= 200 && tMean <= 3500 ); // bounds on NASA polynomials

    fml.deallocate_fields();
  } // number of points

  return status.ok();
}

//==============================================================================

int main( int iarg, char* carg[] )
{
  std::string inputFileName = "h2o2.xml";
  std::string inpGroup;
  bool timings = false;
  bool print   = false;
  size_t nSteps = 5000;
  double dt = 1e-10;

  po::options_description desc("Supported Options");
  desc.add_options()
               ( "help", "print help message" )
               ( "xml-input-file", po::value<std::string>(&inputFileName)->default_value("h2o2.xml"), "Cantera xml input file name" )
               ( "phase", po::value<std::string>(&inpGroup), "name of phase in Cantera xml input file" )
               ( "timings", "Generate comparison timings between Cantera and PoKiTT across several problem sizes" )
               ( "print-fields", "Print field values for Temperature, mass of H2, and mass of OH every 2000 time steps")
               ( "nsteps", po::value<size_t>(&nSteps)->default_value(5000), "How many time steps to take" )
               ( "dt",     po::value<double>(&dt)->default_value(1e-10),    "Size of time steps (s) to take"  );

  // parse the command line options input describing the problem
  try {
    po::variables_map args;
    po::store( po::parse_command_line(iarg,carg,desc), args );
    po::notify(args);

    print = args.count("print-fields") > 0;
#   ifdef ENABLE_CUDA
    if( print )
      throw std::runtime_error("print-fields is not (yet) supported with CUDA");
#   endif

    timings = args.count("timings") > 0;

    if( timings && args["nsteps"].defaulted() )
      nSteps = 5;

    if (args.count("help")) {
      std::cout << desc << "\n";
      return -1;
    }
  }
  catch( std::exception& err ){
    std::cout << "Error parsing input arguments\n" << err.what() << std::endl
        << std::endl << "Usage:\n" << desc << std::endl;
    return -2;
  }

  try{
    const CanteraObjects::Setup setup( "Mix", inputFileName, inpGroup );
    CanteraObjects::setup_cantera( setup );

    TestHelper status( !timings );
    status( driver( timings, print, nSteps, dt ), "Reaction Diffusion" );

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
