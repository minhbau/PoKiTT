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
#include "SpeciesRHS.h"
#include "EnthalpyRHS.h"
#include "HeatFlux.h"
#include "MassFractions.h"
#include "SpeciesDiffusion.h"
#include <pokitt/thermo/Temperature.h>
#include <pokitt/thermo/Pressure.h>
#include <pokitt/thermo/Density.h>
#include <pokitt/MixtureMolWeight.h>
#include <pokitt/kinetics/ReactionRates.h>
#include <pokitt/thermo/Enthalpy.h>
#include <pokitt/transport/DiffusionCoeffMix.h>
#include <pokitt/transport/ThermalCondMix.h>

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


  const Tag xTag  ( "XCoord",               Expr::STATE_NONE );
  const Tag tTag  ( "Temp",                 Expr::STATE_NONE );
  const Tag pTag  ( "Pressure",             Expr::STATE_NONE);
  const Tag rhoTag( "Density",              Expr::STATE_NONE);
  const Tag mmwTag( "mix mol weight",       Expr::STATE_NONE);
  const Tag hTag  ( "enthalpy",             Expr::STATE_N);
  const Tag lamTag( "thermal conductivity", Expr::STATE_NONE);
  const Tag qXTag ( "Heat Flux X",          Expr::STATE_NONE);
  const Tag qYTag ( "Heat Flux Y",          Expr::STATE_NONE);
  const Tag hRHSTag( "enthalpy RHS",        Expr::STATE_N);
  TagList hiTags;
  TagList jXTags;
  TagList jYTags;
  TagList dTags;
  TagList yiTags;
  TagList rhoYiTags;
  TagList specRHSTags;
  TagList rTags;
  for( size_t n=0; n<nSpec; ++n ){
    std::string spec = boost::lexical_cast<std::string>(n);
    yiTags.push_back( Tag( "yi_" + spec, Expr::STATE_NONE ) );
    hiTags.push_back( Tag( "hi_" + spec, Expr::STATE_NONE ) );
    rTags.push_back(  Tag( "ri_" + spec, Expr::STATE_NONE ) );
    dTags.push_back(  Tag( "di_" + spec, Expr::STATE_NONE ) );
    jXTags.push_back(  Tag( "Jxi_" + spec, Expr::STATE_NONE ) );
    jYTags.push_back(  Tag( "Jyi_" + spec, Expr::STATE_NONE ) );
    if( n != (nSpec-1) ){
      specRHSTags.push_back( Tag( "specRHS_" + spec, Expr::STATE_N ) );
      rhoYiTags.push_back( Tag( "rhoYi_" + spec, Expr::STATE_N ) );
    }
  }

  const double rho0 = 0.3;
  const double t0 = 800;
  const double length = 1e-3;

  Expr::ExpressionFactory initFactory;
  std::set<Expr::ExpressionID> initIDs;
  Expr::ExpressionID hID; // temperature
  {
    typedef Expr::PlaceHolder   <CellField>::Builder XCoord;
    typedef Expr::ConstantExpr  <CellField>::Builder Density;
    typedef Expr::LinearFunction<CellField>::Builder rhoYi;
    typedef Expr::LinearFunction<CellField>::Builder Temperature;
    typedef Enthalpy            <CellField>::Builder Enthalpy;
    typedef MassFractions       <CellField>::Builder MassFracs;

    initFactory.register_expression( new XCoord( xTag )    );
    initFactory.register_expression( new Temperature( tTag, xTag, 0.0, t0 )    );
    initFactory.register_expression( new Density( rhoTag, rho0 )    );
    initFactory.register_expression( new rhoYi ( rhoYiTags[0], rhoTag, 0.01, 0.0  ) );
    initFactory.register_expression( new rhoYi ( rhoYiTags[3], rhoTag, 0.20, 0.0  ) );
    for( size_t n = 0; n < (nSpec-1); ++n ){
      if( n != 0 && n != 3 )
        initFactory.register_expression( new rhoYi ( rhoYiTags[n], rhoTag, 1e-8, 0.0  ) );
    }
    initFactory.register_expression( new MassFracs( yiTags, rhoTag, rhoYiTags) );
    hID = initFactory.register_expression( new Enthalpy ( hTag , tTag, yiTags )    );
    initIDs.insert( hID );
  }
  Expr::ExpressionTree initTree( initIDs, initFactory, 0 );

  Expr::ExpressionFactory execFactory;
  {
    typedef Expr::PlaceHolder<CellField>::Builder Density;

    typedef MassFractions       <CellField>::Builder MassFracs;
    typedef MixtureMolWeight    <CellField>::Builder MixMolWeight;
    typedef Temperature         <CellField>::Builder Temperature;
    typedef Pressure            <CellField>::Builder Pressure;
    typedef SpeciesEnthalpy     <CellField>::Builder SpecEnthalpy;
    typedef ThermalConductivity <CellField>::Builder ThermCond;
    typedef ReactionRates       <CellField>::Builder ReactionRates;
    typedef DiffusionCoeff      <CellField>::Builder DiffusionCoeffs;

    typedef SpeciesDiffFlux <XFluxT>::Builder MassFluxX;
    typedef SpeciesDiffFlux <YFluxT>::Builder MassFluxY;
    typedef HeatFlux        <XFluxT>::Builder HeatFluxX;
    typedef HeatFlux        <YFluxT>::Builder HeatFluxY;

    typedef EnthalpyRHS <CellField>::Builder EnthalpyRHS;
    typedef SpeciesRHS  <CellField>::Builder SpeciesRHS;

    execFactory.register_expression( new Density( rhoTag ) );
    execFactory.register_expression( new MassFracs( yiTags, rhoTag, rhoYiTags) );
    execFactory.register_expression( new MixMolWeight ( mmwTag, yiTags )           );
    execFactory.register_expression( new Temperature ( tTag, yiTags, hTag ) );
    execFactory.register_expression( new Pressure ( pTag, tTag, rhoTag, mmwTag) );
    for( size_t n = 0; n < nSpec; ++n ){
      execFactory.register_expression( new SpecEnthalpy ( hiTags[n], tTag, n ) );
    }
    execFactory.register_expression( new ThermCond ( lamTag, tTag, yiTags, mmwTag ) );
    execFactory.register_expression( new ReactionRates(rTags, tTag, pTag, yiTags, mmwTag) );
    execFactory.register_expression( new DiffusionCoeffs(dTags, tTag, pTag, yiTags, mmwTag) );

    execFactory.register_expression( new MassFluxX( jXTags, yiTags,rhoTag, mmwTag, dTags ) );
    execFactory.register_expression( new MassFluxY( jYTags, yiTags,rhoTag, mmwTag, dTags ) );
    execFactory.register_expression( new HeatFluxX( qXTag, tTag, lamTag, hiTags, jXTags ) );
    execFactory.register_expression( new HeatFluxY( qYTag, tTag, lamTag, hiTags, jYTags ) );

    execFactory.register_expression( new EnthalpyRHS ( hRHSTag, rhoTag, tag_list( qXTag, qYTag ) ) );
    for( size_t n = 0; n < (nSpec-1); ++n ){
      execFactory.register_expression( new SpeciesRHS ( specRHSTags[n], rTags[n], tag_list( jXTags[n], jYTags[n] ) ) );
    }
  }

  std::vector< IntVec > sizeVec;
  if( timings ){
    sizeVec.push_back( IntVec(126,126, 1 ) );
    sizeVec.push_back( IntVec( 62, 62, 1 ) );
    sizeVec.push_back( IntVec( 30, 30, 1 ) );
    sizeVec.push_back( IntVec( 14, 14, 1 ) );
    sizeVec.push_back( IntVec(  6,  6, 1 ) );
  }
  else{
    sizeVec.push_back( IntVec( 10,  5,  1) );
  }

  for( std::vector< IntVec >::iterator iSize = sizeVec.begin(); iSize!= sizeVec.end(); ++iSize){

    IntVec gridSize = *iSize;
    Expr::ExprPatch patch( gridSize[0], gridSize[1], gridSize[2] );

    Expr::FieldManagerList& fml = patch.field_manager_list();
    initTree.register_fields( fml );
    initTree.bind_fields( fml );
    initTree.lock_fields( fml );

    fml.allocate_fields( patch.field_info() );
    Expr::TimeStepper timeIntegrator( execFactory, Expr::FORWARD_EULER, patch.id() );

    SO::GhostData ghosts( IntVec( 1, 1, 0), IntVec( 1, 1, 0 ) ); // 1 on +-x and +- y and 0 on z
    for( size_t n = 0; n < (nSpec-1); ++n ){
      timeIntegrator.add_equation<CellField>( (rhoYiTags[n].name()), specRHSTags[n], ghosts );
    }
    timeIntegrator.add_equation<CellField>( hTag.name(), hRHSTag, ghosts );

    SO::Grid grid( gridSize, SO::DoubleVec( length, length, 1 ) );
    SO::OperatorDatabase& opDB = patch.operator_database();
    SO::build_stencils( grid, opDB );

    setup_bcs( execFactory, grid, ghosts, rhoYiTags[0], 0, SO::NEUMANN );
    setup_bcs( execFactory, grid, ghosts, rhoYiTags[3], 0, SO::NEUMANN );
    for( size_t n = 0; n < (nSpec-1); ++n ){
      if( n != 0 && n != 3 )
        setup_bcs( execFactory, grid, ghosts, rhoYiTags[n], 0, SO::NEUMANN );
    }
    setup_bcs( execFactory, grid, ghosts, tTag, t0+500, SO::DIRICHLET );

    timeIntegrator.finalize( fml, patch.operator_database(), patch.field_info() );
    {
      std::ofstream out( "rxn_example.dot" );
      timeIntegrator.get_tree()->write_tree( out );
    }
    timeIntegrator.get_tree()->lock_fields(fml);

    CellField& xcoord = fml.field_ref< CellField >( xTag );
    grid.set_coord<SO::XDIR>( xcoord );
#   ifdef ENABLE_CUDA
    xcoord.set_device_as_active( GPU_INDEX );
#   endif

    initTree.execute_tree();
    CellField& t = fml.field_ref< CellField >( tTag );
    CellField& h = fml.field_ref< CellField >( hTag );
    Timer timer;
    timer.start();
    for( size_t s = 0; s < nSteps; ++s ){
      if( s%5000 == 0 ){
        {
          int n=0;
          std::cout<<"n = " << n << "\n";
          SO::interior_print_field( fml.field_ref< CellField >( rhoYiTags[n] ), std::cout );
          n=3;
          std::cout<<"n = " << n << "\n";
          SO::interior_print_field( fml.field_ref< CellField >( rhoYiTags[n] ), std::cout );
          n=5;
          std::cout<<"n = " << n << "\n";
          SO::interior_print_field( fml.field_ref< CellField >( rhoYiTags[n] ), std::cout );
        }
        std::cout<<"h\n";
        SO::interior_print_field( h, std::cout );
        std::cout<<"t\n";
        SO::interior_print_field( t, std::cout );
      }
      timeIntegrator.step( dt );
    }
    timer.stop();

    const double tMean = SO::nebo_sum_interior( t ) / ( gridSize[0] * gridSize[1] * gridSize[2] );
    status( tMean >= 300 && tMean < 5000 );

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
