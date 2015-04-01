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
#include <numeric>

#include <expression/ExprLib.h>

#include <spatialops/structured/Grid.h>
#include <spatialops/OperatorDatabase.h>
#include <spatialops/structured/FieldComparisons.h>
#include <spatialops/util/TimeLogger.h>

#include <pokitt/CanteraObjects.h>
#include "PoKiTTVersion.h"

#include "test/TestHelper.h"

#include "EnthalpyTransport.h"
#include "SpeciesTransport.h"
#include "TagManager.h"

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>

#include <cantera/kernel/global.h>
#include <cantera/kernel/ctexceptions.h>

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
    SO::print_field( field, std::cout, true );
  }
}

void print_matlab( Expr::FieldManagerList& fml, const TagList& fieldTags, const int s ){
  BOOST_FOREACH( const Tag iTag, fieldTags){

    CellField& field = fml.field_ref< CellField >( iTag );
#ifdef ENABLE_CUDA
   field.validate_device_async( GPU_INDEX );
#endif
    SO::write_matlab( field, iTag.name()+boost::lexical_cast<std::string>(s), false);
  }
}

bool driver( const bool timings,
             const bool print,
             const bool matlab,
             const size_t nSteps,
             const double dt,
             const size_t maxPoints)
{
  TestHelper status( !timings );
  const int nSpec = CanteraObjects::number_species();

  const double rho0 = 1.0;
  const double length = 0.1;
  const double tMax = 1700;
  const double tDev = 0.1*length;
  const double tMean = length/2;
  const double tBase = 800;
  std::vector<double> yi(nSpec, 1e-12);
  std::vector<double> slope( nSpec, 1e-12 );
  yi[28]=0.25;
  slope[28]=0.4;
  yi[3]=0.25;
  slope[3]=-0.4;

  const TagManager tagMgr( nSpec,           // enum values:
      Tag( "rhoYi",     Expr::STATE_N    ), // RHOYI
      Tag( "enthalpy",  Expr::STATE_N    ), // H
      Tag( "yi",        Expr::STATE_NONE ), // YI
      Tag( "temp",      Expr::STATE_NONE ), // T
      Tag( "pressure",  Expr::STATE_NONE ), // P
      Tag( "density",   Expr::STATE_NONE ), // RHO
      Tag( "mixMW ",   Expr::STATE_NONE ), // MMW
      Tag( "r",         Expr::STATE_NONE ), // R
      Tag( "lambda",    Expr::STATE_NONE ), // LAM
      Tag( "D",         Expr::STATE_NONE ), // D
      Tag( "J",         Expr::STATE_NONE ), // J
      Tag( "heatFlux", Expr::STATE_NONE ), // Q
      Tag( "XCoord",    Expr::STATE_NONE ), // XCOORD
      Tag( "YCoord",    Expr::STATE_NONE ) // YCOORD
  );

  size_t maxXDim = (maxPoints) - (maxPoints)%10;
  std::vector< IntVec > sizeVec;
  if( timings ){
    sizeVec.push_back( IntVec( maxXDim/10 - 2,  8,  1 ) );
    if( 1024 * 512 < maxPoints ) sizeVec.push_back( IntVec(1022, 510, 1 ) );
    if( 512  * 512 < maxPoints ) sizeVec.push_back( IntVec(510,  510, 1 ) );
    if( 256  * 256 < maxPoints ) sizeVec.push_back( IntVec(254,  254, 1 ) );
    if( 128  * 128 < maxPoints ) sizeVec.push_back( IntVec(126,  126, 1 ) );
    if( 64   * 64  < maxPoints ) sizeVec.push_back( IntVec( 62,  62,  1 ) );
    if( 32   * 32  < maxPoints ) sizeVec.push_back( IntVec( 30,  30,  1 ) );
  }
  sizeVec.push_back( IntVec( 126,  126,  1 ) );
  for( std::vector< IntVec >::iterator iSize = sizeVec.begin(); iSize!= sizeVec.end(); ++iSize){

    typedef EnthalpyTransport<CellField> EnthalpyTransport;
    typedef SpeciesTransport <CellField> SpeciesTransport;

    Expr::ExpressionFactory initFactory;
    Expr::ExpressionFactory execFactory;

    std::set<Expr::ExpressionID> initIDs;
    EnthalpyTransport* hTrans = new EnthalpyTransport( execFactory, tagMgr );
    initIDs.insert( hTrans->initial_condition( initFactory, tMax, tBase, tMean, tDev) );

    typedef typename Expr::ConstantExpr   <CellField>::Builder Rho;
    initIDs.insert( initFactory.register_expression( new Rho( tagMgr[RHO],    rho0 ) ) );

    std::list<SpeciesTransport*> specEqns;
    for( size_t n=0; n<(nSpec-1); ++n ){
      SpeciesTransport* specTrans = new SpeciesTransport( execFactory, tagMgr, n );
      yi[n] = fabs(yi[n] - 0.5 * slope[n] ) + 1e-12;
      initIDs.insert( specTrans->initial_condition( initFactory, yi[n], slope[n]/length, rho0 ) );
      specEqns.push_back( specTrans );
    }
    specEqns.front()->register_one_time_expressions( execFactory );

    Expr::ExpressionTree initTree( initIDs, initFactory, 0 );

    IntVec gridSize = *iSize;
    Expr::ExprPatch patch( gridSize[0], gridSize[1], gridSize[2] );
    SO::Grid grid( gridSize, SO::DoubleVec( length, length, 1 ) );
    SO::GhostData ghosts( IntVec( 1, 1, 0), IntVec( 1, 1, 0 ) ); // 1 on +-x and +- y and 0 on z

    Expr::TimeStepper timeIntegrator( execFactory, Expr::FORWARD_EULER, "timestepper", patch.id() );
    timeIntegrator.request_timings();

    for( std::list<SpeciesTransport*>::iterator iEqn=specEqns.begin(); iEqn!=specEqns.end(); ++iEqn ){
      timeIntegrator.add_equation<CellField>( (*iEqn)->solution_variable_name(), (*iEqn)->get_rhs_tag(), ghosts );
      (*iEqn)->setup_boundary_conditions( grid, execFactory );
    }
    timeIntegrator.add_equation<CellField>( hTrans->solution_variable_name(), hTrans->get_rhs_tag(), ghosts );
    hTrans->setup_boundary_conditions( grid, execFactory );

    Expr::FieldManagerList& fml = patch.field_manager_list();

    fml.allocate_fields( patch.field_info() );
    SO::OperatorDatabase& opDB = patch.operator_database();
    SO::build_stencils( grid, opDB );

    timeIntegrator.finalize( fml, patch.operator_database(), patch.field_info() );
//    {
//      std::ofstream out( "rxn_example.dot" );
//      timeIntegrator.get_tree()->write_tree( out );
//    }
    timeIntegrator.get_tree()->lock_fields(fml);

    initTree.register_fields( fml );
    initTree.bind_fields( fml );
    initTree.lock_fields( fml );

    CellField& xcoord = fml.field_ref< CellField >( tagMgr[XCOORD] );
    CellField& ycoord = fml.field_ref< CellField >( tagMgr[YCOORD] );
    grid.set_coord<SO::XDIR>( xcoord );
    grid.set_coord<SO::YDIR>( ycoord );
#   ifdef ENABLE_CUDA
    initTree.set_device_index( GPU_INDEX, fml );
    timeIntegrator.get_tree()->set_device_index( GPU_INDEX, fml );
#   endif
    initTree.execute_tree();
    int s = -1;
    SpatialOps::Timer timer;
    std::vector< double > times;
    if( timings ){
      std::cout << "\nPoKiTT Reaction Diffusion size " << xcoord.window_with_ghost().glob_npts() << std::endl;
      timeIntegrator.request_timings();
    }
    try{
      if( timings ) timeIntegrator.step( dt ); // sets high water mark on memory
      for( s = 0; s <= nSteps; ++s ){
        if( s%5000 == 0 && print){
          std::cout<<"Fields at time "<< s*dt << "; step " << s << std::endl;
          print_fields( fml, tag_list( tagMgr[T], tagMgr[R_N][28], tagMgr[RHOYI_N][28], tagMgr[RHOYI_N][3] ) );
        }
        if( s%10000 == 0 && matlab ){
          std::cout<<"Fields at time "<< s*dt << "; step " << s << std::endl;
          print_matlab( fml, tag_list( tagMgr[T], tagMgr[R_N][28], tagMgr[RHOYI_N][28], tagMgr[RHOYI_N][3] ), s );
        }
        timer.reset();
        timeIntegrator.step( dt );
        times.push_back( timer.stop() );
      }
    }
    catch( std::runtime_error& err ){
      std::ostringstream msg;
      msg << " \n Error occured during iteration s = " << s << std::endl
          << err.what() << std::endl;
      throw std::runtime_error( msg.str() );
    }

    if( timings ){
      std::sort( times.begin(), times.end() );
      const int chop = floor(nSteps/4);
      const double avgTime = std::accumulate( times.begin() + chop, times.end()-chop, 0.0 )/( nSteps + 1 - 2*chop);
      std::cout << "PoKiTT Reaction Diffusion time per step " << avgTime << std::endl;
    }

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
  bool matlab = false;
  size_t nSteps;
  size_t maxPoints;
  double dt;

  po::options_description desc("Supported Options");
  desc.add_options()
               ( "help", "print help message" )
               ( "xml-input-file,i", po::value<std::string>(&inputFileName)->default_value("methanol.xml"), "Cantera xml input file name" )
               ( "phase", po::value<std::string>(&inpGroup), "name of phase in Cantera xml input file" )
               ( "timings,t", "Generate comparison timings between Cantera and PoKiTT across several problem sizes" )
               ( "print-fields,p", "Print field values for Temperature, mass of H2, and mass of OH every 5000 time steps")
               ( "matlab", "Save field values for Temperature, mass of H2, and mass of OH every 10000 time steps")
               ( "nsteps,n", po::value<size_t>(&nSteps)->default_value(5000), "How many time steps to take" )
               ( "dt",     po::value<double>(&dt)->default_value(1e-15),    "Size of time steps (s) to take"  )
               ( "max-points",     po::value<size_t>(&maxPoints)->default_value(1024*1024),    "Run problems at this size or lower (rounded down to nearest 10)"  );

  // parse the command line options input describing the problem
  try {
    po::variables_map args;
    po::store( po::parse_command_line(iarg,carg,desc), args );
    po::notify(args);

    print = args.count("print-fields");
#   ifdef ENABLE_CUDA
    if( print )
      throw std::runtime_error("print-fields is not (yet) supported with CUDA");
#   endif

    timings = args.count("timings");
    matlab = args.count("matlab");
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

    {
      std::cout << "Reaction-diffusion example\n"
          <<"\n--------------------------------------------------------------\n"
          << "BUILD INFORMATION:\n"
          << "\n\t PoKiTT     git hash: \t" << PoKiTTVersionHash
          << "\n\t ExprLib    git hash: \t" << EXPR_REPO_HASH
          << "\n\t SpatialOps git hash: \t" << SOPS_REPO_HASH
          <<"\n--------------------------------------------------------------\n"
          << "MECHANISM INFORMATION\n"
          << "\n\tXML input file: " << inputFileName
          << "\n\tPhase name    : " << CanteraObjects::phase_name()
          << "\n\t# species     : " << CanteraObjects::number_species()
          << "\n\t# reactions   : " << CanteraObjects::number_rxns()
          <<"\n--------------------------------------------------------------\n\n";
    }

    TestHelper status( !timings );
    status( driver( timings, print, matlab, nSteps, dt, maxPoints ), "Reaction Diffusion" );

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
