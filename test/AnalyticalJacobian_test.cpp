/*
 * The MIT License
 *
 * Copyright (c) 2016 The University of Utah
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

/**
 *  \file   AnalyticalJacobian_test.cpp
 *  \date   May 13, 2016
 *  \author mike
 */


#include <numeric>
#include <iostream>
#include <vector>
#include <cmath>

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include "TestHelper.h"

#include <spatialops/Nebo.h>
#include <spatialops/structured/FieldComparisons.h>
#include <spatialops/structured/Grid.h>
#include <spatialops/structured/MatVecFields.h>
#include <spatialops/structured/MatVecOps.h>
#include <spatialops/util/TimeLogger.h>
#include <spatialops/util/MemoryUsage.h>

#include <expression/ExprLib.h>
#include <expression/Expression.h>

#include <cantera/IdealGasMix.h>

#include <pokitt/CanteraObjects.h>
#include <pokitt/kinetics/AnalyticalJacobian.h>
#include <pokitt/transformations/ConservedPrimitive.h>
#include "LinearMassFracs.h"
#include <pokitt/MixtureMolWeight.h>
#include <pokitt/kinetics/ReactionRates.h>
#include <pokitt/thermo/Density.h>
#include <pokitt/thermo/InternalEnergy.h>
#include <pokitt/thermo/HeatCapacity_Cv.h>

using namespace pokitt;
namespace so = SpatialOps;
namespace po = boost::program_options;

#ifdef ENABLE_CUDA
# define LOCATION GPU_INDEX
#else
# define LOCATION CPU_INDEX
#endif

namespace SO = SpatialOps;
typedef SO::SVolField  CellField;
typedef SO::SpatFldPtr<CellField> CellFieldPtrT;


void get_timings( size_t n, size_t numberOfRepeats, std::string inputFileName )
{
  // ---------------- build the grid, setup a timer ----------------
  const SpatialOps::IntVec    npts   ( n, n, 1 );
  const SpatialOps::DoubleVec lengths( 1, 1, 1 );
  const so::Grid              grid( npts, lengths );
  const so::DoubleVec         spacings = grid.spacing();
  const so::GhostData         nghost(0);
  const so::BoundaryCellInfo  bcInfo = so::BoundaryCellInfo::build<CellField>( npts[0]>1, npts[1]>1, npts[2]>1 );
  const so::MemoryWindow      window( get_window_with_ghost( npts, nghost, bcInfo ) );

  CellField xCoordinate( window, bcInfo, nghost, NULL, so::InternalStorage, LOCATION );
  CellField yCoordinate( window, bcInfo, nghost, NULL, so::InternalStorage, LOCATION );
  CellField zCoordinate( window, bcInfo, nghost, NULL, so::InternalStorage, LOCATION );

  grid.set_coord<so::XDIR>(xCoordinate);
  grid.set_coord<so::YDIR>(yCoordinate);
  grid.set_coord<so::ZDIR>(zCoordinate);


#   ifdef ENABLE_CUDA
    xCoordinate.set_device_as_active( GPU_INDEX );
    yCoordinate.set_device_as_active( GPU_INDEX );
    zCoordinate.set_device_as_active( GPU_INDEX );
#   endif

  const std::string nrStr = boost::lexical_cast<std::string>( numberOfRepeats );
  const std::string nStr = boost::lexical_cast<std::string>( n );
  const std::string repeatAndGridStr = "--" + nrStr + "-repeats--" + nStr + "x" + nStr + "-grid";
  SpatialOps::TimeLogger timeLogger( "Analytical-Jacobian-timings" + repeatAndGridStr + "-" + inputFileName );

  // ---------------- initialize the domain ----------------

  const size_t ns = CanteraObjects::number_species();

  std::vector< double > Msp = CanteraObjects::molecular_weights();

  std::vector< double > invMsp( ns );
  for( size_t s=0; s<ns; ++s ){
    invMsp[s] = 1.0 / Msp[s];
  }
  const double Ru    = CanteraObjects::gas_constant();
  const double invRu = 1.0/Ru;

  CellField T ( window, bcInfo, nghost, NULL, so::InternalStorage, LOCATION ); // temperature

  so::SpatFldPtr<CellField> pPtr    = so::SpatialFieldStore::get<CellField>( T ); CellField& p = *pPtr;       // pressure
  so::SpatFldPtr<CellField> rhoPtr  = so::SpatialFieldStore::get<CellField>( T ); CellField& rho  = *rhoPtr;  // density
  so::SpatFldPtr<CellField> MmixPtr = so::SpatialFieldStore::get<CellField>( T ); CellField& Mmix = *MmixPtr; // mixture molecular weight
  so::SpatFldPtr<CellField> RmixPtr = so::SpatialFieldStore::get<CellField>( T ); CellField& Rmix = *RmixPtr; // specific mixture gas constant

  std::vector< so::SpatFldPtr<CellField> > YPtr;  // mass fractions
  std::vector< so::SpatFldPtr<CellField> > XPtr;  // mole fractions
  for( size_t i=0; i<ns; ++i ){
    YPtr.push_back( so::SpatialFieldStore::get<CellField>( T ) ); *YPtr[i] <<= 0;
    XPtr.push_back( so::SpatialFieldStore::get<CellField>( T ) ); *XPtr[i] <<= 0;
  }

  T <<= 1000;   // system temperature (K)
  p <<= 101325; // system pressure (Pa)

  for( size_t i=0; i<ns; ++i ) *XPtr[i] <<= 1.0; // equimolar mixture

  Mmix <<= 0;
  for( size_t i=0; i<ns; ++i ){
    *XPtr[i] <<= *XPtr[i] / ( (double) ns );
    Mmix <<= Mmix + *XPtr[i] * Msp[i];
  }

  Rmix <<= Ru / Mmix;
  for( size_t i=0; i<ns; ++i ){
    *YPtr[i] <<= *XPtr[i] * Msp[i] / Mmix;
  }

  rho <<= p / T / Rmix;

  // ---------------- calculate net species mass production rates ----------------
  Expr::ExprPatch rxnPatch( npts[0], npts[1], npts[2] );
  Expr::FieldManagerList& rxnFml = rxnPatch.field_manager_list();
  Expr::ExpressionFactory rxnFactory;
  std::set<Expr::ExpressionID> rxnRoots;

  Expr::TagList wTags;
  for( size_t i=0; i<ns; ++i ){
    wTags.push_back( Expr::Tag( "w_" + boost::lexical_cast<std::string>( i ), Expr::STATE_NONE ) );
  }
  const Expr::Tag rhoTag( "rho", Expr::STATE_NONE );
  const Expr::Tag mmwTag( "Mmix", Expr::STATE_NONE );
  const Expr::Tag tempTag( "temperature", Expr::STATE_NONE );
  Expr::TagList yTags;
  for( size_t i=0; i<ns; ++i ){
    yTags.push_back( Expr::Tag( "y_" + boost::lexical_cast<std::string>( i ), Expr::STATE_NONE ) );
  }

  for( size_t i=0; i<ns; ++i ){
    const Expr::ExpressionID yIDrxn = rxnFactory.register_expression( new Expr::PlaceHolder<CellField>::Builder( yTags[i] ) );
    rxnRoots.insert( yIDrxn );
  }
  const Expr::ExpressionID rhoIDrxn = rxnFactory.register_expression( new Expr::PlaceHolder<CellField>::Builder( rhoTag ) );
  rxnRoots.insert( rhoIDrxn );
  const Expr::ExpressionID mmwIDrxn = rxnFactory.register_expression( new Expr::PlaceHolder<CellField>::Builder( mmwTag ) );
  rxnRoots.insert( mmwIDrxn );
  const Expr::ExpressionID tempIDrxn = rxnFactory.register_expression( new Expr::PlaceHolder<CellField>::Builder( tempTag ) );
  rxnRoots.insert( tempIDrxn );
  const Expr::ExpressionID specRIDrxn = rxnFactory.register_expression( new pokitt::ReactionRates<CellField>::Builder( wTags, tempTag, rhoTag, yTags, mmwTag ) );
  rxnRoots.insert( specRIDrxn );

  Expr::ExpressionTree rxnTree( rxnRoots, rxnFactory, rxnPatch.id() );
  rxnTree.register_fields( rxnFml );
  rxnTree.bind_fields( rxnFml );
  rxnFml.allocate_fields( rxnPatch.field_info() );

  rxnFml.field_ref<CellField>( tempTag ) <<= T;
  rxnFml.field_ref<CellField>( rhoTag )  <<= rho;
  rxnFml.field_ref<CellField>( mmwTag )  <<= Mmix;
  for( size_t i=0; i<ns; ++i ){
    rxnFml.field_ref<CellField>( yTags[i] ) <<= *YPtr[i];
  }

  timeLogger.start( "production rates" );
  for( size_t nTimes=0; nTimes<numberOfRepeats; ++nTimes )
    rxnTree.execute_tree();
  timeLogger.stop( "production rates" );

  // ---------------- calculate net species mass production rates and sensitivities ----------------
  so::FieldVector<CellField> productionRates       ( ns+1, window, bcInfo, nghost, LOCATION );
  const std::size_t memoryBeforeJacobian = SpatialOps::get_memory_usage();
  so::FieldMatrix<CellField> primitiveSensitivities( ns+1, window, bcInfo, nghost, LOCATION );
  const std::size_t memoryAfterJacobian = SpatialOps::get_memory_usage();
  std::cout << std::endl << "Memory usage of rates tree on " << n << "x" << n << " grid: " << memoryBeforeJacobian / 1000 << "kB" << std::endl
                         << "Rates tree plus Jacobian matrix: " << memoryAfterJacobian / 1000 << "kB" << std::endl;

#   ifdef ENABLE_CUDA
  for( size_t prodRateIdx=0; prodRateIdx<ns; ++prodRateIdx ){
    productionRates(prodRateIdx).set_device_as_active( GPU_INDEX );
  }
#   endif
#   ifdef ENABLE_CUDA
  for( size_t prodRateIdx=0; prodRateIdx<ns+1; ++prodRateIdx ){
    for( size_t primVarIdx=0; primVarIdx<ns+1; ++primVarIdx ){
      primitiveSensitivities(prodRateIdx,primVarIdx).set_device_as_active( GPU_INDEX );
    }
  }
#   endif

  ChemicalSourceJacobian csj;

  std::vector< const CellField* > YPtrForJac;
  for( size_t i=0; i<YPtr.size(); ++i )
    YPtrForJac.push_back( &( *YPtr[i] ) );

  timeLogger.start( "production rates + Jacobian" );
  for( size_t nTimes=0; nTimes<numberOfRepeats; ++nTimes )
    csj.evaluate_rates_and_jacobian( primitiveSensitivities, productionRates, T, rho, YPtrForJac, Mmix );
  timeLogger.stop( "production rates + Jacobian" );

  // ---------------- form the chemical rhs and primitive Jacobian for the conserved variables ----------------
  so::FieldVector<CellField> chemRhs( ns+1, window, bcInfo, nghost, LOCATION );
  so::FieldMatrix<CellField> chemPrimitiveJac( ns+1, window, bcInfo, nghost, LOCATION );

#   ifdef ENABLE_CUDA
  for( size_t prodRateIdx=0; prodRateIdx<ns+1; ++prodRateIdx ){
    chemRhs(prodRateIdx).set_device_as_active( GPU_INDEX );
    for( size_t primVarIdx=0; primVarIdx<ns+1; ++primVarIdx ){
      chemPrimitiveJac(prodRateIdx,primVarIdx).set_device_as_active( GPU_INDEX );
    }
  }
#   endif

  chemRhs(0) <<= 0.0;
  chemRhs(1) <<= 0.0;
  for( size_t i=0; i<ns-1; ++i )
    chemRhs(i) <<= productionRates(i);

  for( size_t j=0; j<ns+1; ++j ){
    chemPrimitiveJac(0,j) <<= 0.0;
    chemPrimitiveJac(1,j) <<= 0.0;
  }
  for( size_t i=2; i<ns+1; ++i ){
    for( size_t j=0; j<ns+1; ++j )
      chemPrimitiveJac(i,j) <<= primitiveSensitivities(i-2,j);
  }

  // ---------------- use a state transformation to get the conserved Jacobian for the conserved variables ----------------
  so::FieldMatrix<CellField> dVdU( ns+1, window, bcInfo, nghost, LOCATION );
  so::FieldMatrix<CellField> chemJac( ns+1, window, bcInfo, nghost, LOCATION );

#   ifdef ENABLE_CUDA
  for( size_t prodRateIdx=0; prodRateIdx<ns+1; ++prodRateIdx ){
    for( size_t primVarIdx=0; primVarIdx<ns+1; ++primVarIdx ){
      dVdU(prodRateIdx,primVarIdx).set_device_as_active( GPU_INDEX );
      chemJac(prodRateIdx,primVarIdx).set_device_as_active( GPU_INDEX );
    }
  }
#   endif

  pokitt::ConservedPrimitiveTransformer pct;

  // compute additional fields necessary for transformation
  const Expr::Tag egyTag( "egy", Expr::STATE_NONE );
  const Expr::Tag cvTag( "cv", Expr::STATE_NONE );

  Expr::TagList massTags, specEgyTags;
  for( size_t i=0; i<ns; ++i ){
    massTags.push_back( Expr::Tag( "Y_" + CanteraObjects::species_name( i ), Expr::STATE_NONE ) );
    specEgyTags.push_back( Expr::Tag( "e_" + CanteraObjects::species_name( i ), Expr::STATE_NONE ) );
  }

  Expr::ExpressionFactory factory;
  std::set<Expr::ExpressionID> roots;

  typedef Expr::PlaceHolder<FieldT>::Builder PlaceHolderT;
  typedef pokitt::InternalEnergy<FieldT>::Builder EnergyT;
  typedef pokitt::SpeciesInternalEnergy<FieldT>::Builder SpecEnergyT;
  typedef pokitt::HeatCapacity_Cv<FieldT>::Builder HeatCapT;

  const Expr::ExpressionID tempID = factory.register_expression( new PlaceHolderT( tempTag ) );
  roots.insert( tempID );
  for( size_t i=0; i<ns; ++i ){
    const Expr::ExpressionID massID = factory.register_expression( new PlaceHolderT( massTags[i] ) );
    roots.insert( massID );
    const Expr::ExpressionID specEgyID = factory.register_expression( new SpecEnergyT( specEgyTags[i], tempTag, i ) );
    roots.insert( specEgyID );
  }
  const Expr::ExpressionID cvID = factory.register_expression( new HeatCapT( cvTag, tempTag, massTags ) );
  roots.insert( cvID );
  const Expr::ExpressionID egyID = factory.register_expression( new EnergyT( egyTag, tempTag, massTags ) );
  roots.insert( egyID );

  Expr::ExpressionTree tree( roots , factory, 0 );
  Expr::FieldManagerList fml;
  tree.register_fields( fml );
  tree.bind_fields( fml );
  tree.lock_fields( fml );
  fml.allocate_fields( Expr::FieldAllocInfo( npts, 0, 0, false, false, false ) );

  fml.field_ref<FieldT>( tempTag ) <<= T;
  for( size_t i=0; i<ns; ++i ){
    fml.field_ref<FieldT>( massTags[i] ) <<= *YPtr[i];
  }
  tree.execute_tree();

  const FieldT& cv = fml.field_ref<FieldT>( cvTag );
  const FieldT& etotal = fml.field_ref<FieldT>( egyTag );
  std::vector< const FieldT* > speciesEnergyPtr;
  for( size_t i=0; i<YPtr.size(); ++i ){
    speciesEnergyPtr.push_back( &( fml.field_ref<FieldT>( specEgyTags[i] ) ) );
  }

  timeLogger.start( "state transformation matrix" );
  for( size_t nTimes=0; nTimes<numberOfRepeats; ++nTimes )
    pct.sens_of_primitive_matrix( dVdU, T, rho, YPtrForJac, cv, speciesEnergyPtr );
  timeLogger.stop( "state transformation matrix" );

  timeLogger.start( "dense state transformation" );
  for( size_t nTimes=0; nTimes<numberOfRepeats; ++nTimes )
    chemJac = chemPrimitiveJac * dVdU;
  timeLogger.stop( "dense state transformation" );

  // in-place transformation
  for( size_t i=0; i<chemJac.elements(); ++i ){
    chemJac(i) <<= chemPrimitiveJac(i);
  }
  timeLogger.start( "sparse state transformation" );
  for( size_t nTimes=0; nTimes<numberOfRepeats; ++nTimes )
    pct.right_multiply_to_conserved( chemJac, T, rho, YPtrForJac, cv, speciesEnergyPtr );
  timeLogger.stop( "sparse state transformation" );

  // ---------------- do a direct linear solve with the Jacobian ----------------
  so::FieldVector<CellField> xVector( ns+1, window, bcInfo, nghost, LOCATION );

#   ifdef ENABLE_CUDA
  for( size_t prodRateIdx=0; prodRateIdx<ns+1; ++prodRateIdx ){
    xVector(prodRateIdx).set_device_as_active( GPU_INDEX );
  }
#   endif

  timeLogger.start( "linear solve" );
  for( size_t nTimes=0; nTimes<numberOfRepeats; ++nTimes )
    xVector = chemJac.solve( chemRhs );
  timeLogger.stop( "linear solve" );

  // ---------------- calculate eigenvalues of the Jacobian ----------------
  so::FieldVector<CellField> eigenvalues( ns+1, window, bcInfo, nghost, LOCATION );

#   ifdef ENABLE_CUDA
  for( size_t prodRateIdx=0; prodRateIdx<ns+1; ++prodRateIdx ){
    eigenvalues(prodRateIdx).set_device_as_active( GPU_INDEX );
  }
#   endif

  timeLogger.start( "eigendecomposition" );
  for( size_t nTimes=0; nTimes<numberOfRepeats; ++nTimes )
    eigenvalues = chemJac.eigen_values();
  timeLogger.stop( "eigendecomposition" );

  // ---------------- find maximum local eigenvalue of the Jacobian ----------------
  double maxRealPosEig = 0.0;

  so::SpatFldPtr<CellField> dsPtr     = so::SpatialFieldStore::get<CellField>( T ); CellField& ds     = *dsPtr;
  so::SpatFldPtr<CellField> maxEigPtr = so::SpatialFieldStore::get<CellField>( T ); CellField& maxEig = *maxEigPtr;

#   ifdef ENABLE_CUDA
    ds.set_device_as_active( GPU_INDEX );
    maxEig.set_device_as_active( GPU_INDEX );
#   endif

  const double gesatSafety = 0.1;
  const double gesatRamp   = 1.1;
  const double gesatMax    = 1.0;

  const double eps = std::sqrt( std::numeric_limits<double>::epsilon() );

  timeLogger.start( "GESAT dual time step size" );
  for( size_t nTimes=0; nTimes<numberOfRepeats; ++nTimes ){
    maxEig <<= 0;
    for( size_t i=0; i<eigenvalues.elements(); ++i ){
      maxEig <<= max( maxEig, eigenvalues(i) );
    }
    ds <<= min( min( gesatMax, gesatSafety / maxEig ), gesatRamp * ds );
  }
  timeLogger.stop( "GESAT dual time step size" );

  std::cout << "Problem size: " << n << " x " << n << " x " << 1 << " = " << n*n << " cells, avg over " << numberOfRepeats << " runs" << std::endl
            << "--------------------------------------------------" << std::endl;
  std::cout << " Timings in microseconds per cell" << std::endl
            << " - Production rates: \t " << timeLogger.timer("production rates" ).elapsed_time() * 1.0e6 / n / n / numberOfRepeats << std::endl
            << " - Rates + Jacobian: \t " << timeLogger.timer("production rates + Jacobian" ).elapsed_time() * 1.0e6 / n / n / numberOfRepeats << std::endl
            << " - Dense transform: \t " << ( timeLogger.timer("state transformation matrix" ).elapsed_time() + \
                                              timeLogger.timer("dense state transformation" ).elapsed_time() ) * 1.0e6 / n / n / numberOfRepeats << std::endl
            << " - Sparse transform: \t " << ( timeLogger.timer("state transformation matrix" ).elapsed_time() + \
                                               timeLogger.timer("sparse state transformation" ).elapsed_time() ) * 1.0e6 / n / n / numberOfRepeats << std::endl
            << " - Linear solve: \t " << timeLogger.timer("linear solve" ).elapsed_time() * 1.0e6 / n / n / numberOfRepeats << std::endl
            << " - Eigenvalues: \t " << timeLogger.timer("eigendecomposition" ).elapsed_time() * 1.0e6 / n / n / numberOfRepeats << std::endl
            << " - Gesat step size: \t " << timeLogger.timer("GESAT dual time step size" ).elapsed_time() * 1.0e6 / n / n / numberOfRepeats << std::endl
            << "--------------------------------------------------" << std::endl;
}


//==============================================================================

std::string primitive_variable_string( size_t idx )
{
  if     ( idx == 0 ) return "rho";
  else if( idx == 1 ) return "T";
  else                return "Y_" + CanteraObjects::species_name( idx-2 );
}

//==============================================================================

void generate_jacobian_for_accuracy_test( const std::string& inputFile )
{
  // equimolar mixture at four temperatures, four pressures
  const size_t nT = 4; // number of temperatures
  const size_t nP = 4; // number of pressures

  const double Tmin = 300;  // minimum temperature (K)
  const double Tmax = 3000; // maximum temperature (K)
  const double pmin = 0.01  * 101325;  // minimum pressure (Pa)
  const double pmax = 100.0 * 101325; // maximum pressure (Pa)

  // ---------------- build the grid ----------------
  const SpatialOps::IntVec    npts   ( nT, nP, 1 );
  const SpatialOps::DoubleVec lengths( 1 , 1 , 1 );
  const so::Grid              grid( npts, lengths );
  const so::DoubleVec         spacings = grid.spacing();
  const so::GhostData         nghost(0);
  const so::BoundaryCellInfo  bcInfo = so::BoundaryCellInfo::build<CellField>( npts[0]>1, npts[1]>1, npts[2]>1 );
  const so::MemoryWindow      window( get_window_with_ghost( npts, nghost, bcInfo ) );

  CellField xCoordinate( window, bcInfo, nghost, NULL, so::InternalStorage, LOCATION );
  CellField yCoordinate( window, bcInfo, nghost, NULL, so::InternalStorage, LOCATION );

  grid.set_coord<so::XDIR>(xCoordinate);
  grid.set_coord<so::YDIR>(yCoordinate);

#   ifdef ENABLE_CUDA
  xCoordinate.set_device_as_active( GPU_INDEX );
  yCoordinate.set_device_as_active( GPU_INDEX );
#   endif

  // ---------------- initialize the domain ----------------
  const size_t ns = CanteraObjects::number_species();

  std::vector< double > Msp = CanteraObjects::molecular_weights();

  std::vector< double > invMsp( ns );
  for( size_t s=0; s<ns; ++s ){
    invMsp[s] = 1.0 / Msp[s];
  }
  const double Ru    = CanteraObjects::gas_constant();
  const double invRu = 1.0/Ru;

  CellField T ( window, bcInfo, nghost, NULL, so::InternalStorage, LOCATION ); // temperature

  so::SpatFldPtr<CellField> pPtr    = so::SpatialFieldStore::get<CellField>( T ); CellField& p = *pPtr;       // pressure
  so::SpatFldPtr<CellField> rhoPtr  = so::SpatialFieldStore::get<CellField>( T ); CellField& rho  = *rhoPtr;  // density
  so::SpatFldPtr<CellField> MmixPtr = so::SpatialFieldStore::get<CellField>( T ); CellField& Mmix = *MmixPtr; // mixture molecular weight
  so::SpatFldPtr<CellField> RmixPtr = so::SpatialFieldStore::get<CellField>( T ); CellField& Rmix = *RmixPtr; // specific mixture gas constant

  std::vector< so::SpatFldPtr<CellField> > YPtr;  // mass fractions
  std::vector< so::SpatFldPtr<CellField> > XPtr;  // mole fractions
  for( size_t i=0; i<ns; ++i ){
    YPtr.push_back( so::SpatialFieldStore::get<CellField>( T ) ); *YPtr[i] <<= 0;
    XPtr.push_back( so::SpatialFieldStore::get<CellField>( T ) ); *XPtr[i] <<= 0;
  }

  T <<= Tmin + ( Tmax - Tmin ) * xCoordinate;
  p <<= pmin + ( pmax - pmin ) * yCoordinate;

  for( size_t i=0; i<ns; ++i ) *XPtr[i] <<= 1.0; // equimolar mixture

  Mmix <<= 0;
  for( size_t i=0; i<ns; ++i ){
    *XPtr[i] <<= *XPtr[i] / ( (double) ns );
    Mmix <<= Mmix + *XPtr[i] * Msp[i];
  }

  Rmix <<= Ru / Mmix;
  for( size_t i=0; i<ns; ++i ){
    *YPtr[i] <<= *XPtr[i] * Msp[i] / Mmix;
  }

  rho <<= p / T / Rmix;

  // ---------------- calculate net species mass production rates and sensitivities ----------------
  so::FieldVector<CellField> productionRates       ( ns+1, window, bcInfo, nghost, LOCATION );
  so::FieldMatrix<CellField> primitiveSensitivities( ns+1, window, bcInfo, nghost, LOCATION );

#   ifdef ENABLE_CUDA
  for( size_t prodRateIdx=0; prodRateIdx<ns; ++prodRateIdx ){
    productionRates(prodRateIdx).set_device_as_active( GPU_INDEX );
  }
#   endif
#   ifdef ENABLE_CUDA
  for( size_t prodRateIdx=0; prodRateIdx<ns+1; ++prodRateIdx ){
    for( size_t primVarIdx=0; primVarIdx<ns+1; ++primVarIdx ){
      primitiveSensitivities(prodRateIdx,primVarIdx).set_device_as_active( GPU_INDEX );
    }
  }
#   endif

  ChemicalSourceJacobian csj;

  std::vector< const CellField* > YPtrForJac;
  for( size_t i=0; i<YPtr.size(); ++i )
    YPtrForJac.push_back( &( *YPtr[i] ) );
  csj.evaluate_rates_and_jacobian( primitiveSensitivities, productionRates, T, rho, YPtrForJac, Mmix );

  // ---------------- write the Jacobian to disk ----------------

#ifdef ENABLE_CUDA
  for( size_t prodRateIdx=0; prodRateIdx<ns+1; ++prodRateIdx ){
    for( size_t primVarIdx=0; primVarIdx<ns+1; ++primVarIdx ){
      primitiveSensitivities(prodRateIdx,primVarIdx).add_device( CPU_INDEX );
    }
  }
#endif

  std::vector<std::string> speciesList;
  boost::split( speciesList, inputFile, boost::is_any_of( "." ) );
  std::string inputNoExt = speciesList[0];
  std::string outputFile = "jacobian-" + inputNoExt + ".dat";

  std::ofstream outputStream;

  // clean the file
  outputStream.open( outputFile.c_str() , std::ofstream::out | std::ofstream::trunc );
  outputStream.close();

  // write the Jacobian
  outputStream.open( outputFile.c_str() , std::ofstream::out | std::ofstream::app );

  for( size_t prodRateIdx=0; prodRateIdx<ns; ++prodRateIdx ){
    for( size_t primVarIdx=0; primVarIdx<ns+1; ++primVarIdx ){
      outputStream << "partial w_" << CanteraObjects::species_name( prodRateIdx ) << " / partial " << primitive_variable_string( primVarIdx ) << std::endl;
      print_field( primitiveSensitivities(prodRateIdx,primVarIdx), outputStream );
    }
  }
  std::cout << "- Wrote Jacobian output to " << outputFile << std::endl;
}

//==============================================================================

const std::vector< std::vector<double> >
extract_mass_fracs( const Expr::TagList yiTags, Expr::FieldManagerList& fml ){
  CellField& yi0 = fml.field_ref< CellField >(yiTags[0]);
  const size_t nPts = yi0.window_with_ghost().glob_npts();
  const int nSpec = yiTags.size();

  std::vector< std::vector<double> > massFracs;
  for( size_t i=0; i<nPts; ++i ){
    std::vector<double> massFrac( nSpec, 0.0 );
    massFracs.push_back(massFrac);
  }
  for( size_t n=0; n<nSpec; ++n ){
    CellField& yi = fml.field_ref< CellField >( yiTags[n] );
#   ifdef ENABLE_CUDA
    yi.set_device_as_active( CPU_INDEX );
#   endif
    size_t i=0;
    for( CellField::const_iterator iY = yi.begin(); iY != yi.end(); ++iY, ++i ){
      massFracs[i][n] = *iY;
    }
  }
  return massFracs;
}

//==============================================================================

const std::vector< CellFieldPtrT >
get_cantera_results( Cantera::IdealGasMix& gasMix,
                     Expr::FieldManagerList& fml,
                     const Expr::Tag& tTag,
                     const Expr::TagList& yiTags,
                     const Expr::Tag& pTag){

  const std::vector<double>& molecularWeights = gasMix.molecularWeights();
  const int nSpec = gasMix.nSpecies();

  const std::vector< std::vector<double> > massFracs = extract_mass_fracs( yiTags, fml );

  CellField& temp   = fml.field_ref< CellField >( tTag );
  CellField& press  = fml.field_ref< CellField >( pTag );
# ifdef ENABLE_CUDA
  temp.set_device_as_active  ( CPU_INDEX );
  press.set_device_as_active ( CPU_INDEX );
# endif

  std::vector< CellFieldPtrT > canteraResults;
  for( size_t n=0; n < nSpec; ++n){
    canteraResults.push_back( SO::SpatialFieldStore::get<CellField>(temp) );
  }

  std::vector< std::vector<double> >::const_iterator iMass;
  CellField::const_iterator                          iTemp;
  CellField::const_iterator                          iPress;
  CellField::iterator                                iCant;

  std::vector<double> rResult(nSpec,0.0);
  iPress = press.begin();
  iMass  = massFracs.begin();
  size_t i = 0;
  for( iTemp = temp.begin(); iTemp != temp.end(); ++iTemp, ++iMass, ++i){
    gasMix.setMassFractions_NoNorm( &(*iMass)[0] );
    gasMix.setState_TP( *iTemp, *iPress );
    gasMix.getNetProductionRates(&rResult[0]);
    for( size_t n=0; n<nSpec; ++n){
      (*canteraResults[n])[i] = rResult[n] * molecularWeights[n];
    }
  }

  return canteraResults;
}

bool driver( const double pressure = 101325 )
{
  TestHelper status;
  Cantera::IdealGasMix* const gasMix = CanteraObjects::get_gasmix();
  const int nSpec=gasMix->nSpecies();

  const Expr::Tag xTag  ( "XCoord",      Expr::STATE_NONE );
  const Expr::Tag tTag  ( "Temperature", Expr::STATE_NONE);
  const Expr::Tag pTag  ( "Pressure",    Expr::STATE_NONE);
  const Expr::Tag mmwTag( "mmw",         Expr::STATE_NONE);
  const Expr::Tag rhoTag( "rho",         Expr::STATE_NONE);
  Expr::TagList yiTags;
  Expr::TagList rTags;
  for( size_t n=0; n<nSpec; ++n ){
    yiTags.push_back( Expr::Tag( "yi_" + boost::lexical_cast<std::string>(n), Expr::STATE_NONE ) );
    rTags.push_back(  Expr::Tag( "ri" + boost::lexical_cast<std::string>(n), Expr::STATE_NONE ) );
  }

  // we use an initialization tree to avoid recalculations when timing the execution
  Expr::ExpressionFactory initFactory;
  std::set<Expr::ExpressionID> initIDs;
  Expr::ExpressionID tID; // temperature
  Expr::ExpressionID pID; // pressure
  Expr::ExpressionID mID; // mixture molecular weight
  Expr::ExpressionID yID; // mass fractions
  {
    typedef Expr::PlaceHolder     <CellField>::Builder XCoord;
    typedef       LinearMassFracs <CellField>::Builder MassFracs;
    typedef Expr::ConstantExpr    <CellField>::Builder Pressure;
    typedef       MixtureMolWeight<CellField>::Builder MixMolWeight;
    typedef Expr::LinearFunction  <CellField>::Builder Temperature;
    typedef       Density         <CellField>::Builder Density;

    initFactory.register_expression(       new XCoord       ( xTag )                     );
    yID = initFactory.register_expression( new MassFracs    ( yiTags, xTag )             );
    initFactory.register_expression(       new Pressure     ( pTag, pressure )           );
    mID = initFactory.register_expression( new MixMolWeight ( mmwTag, yiTags, MASS )     );
    tID = initFactory.register_expression( new Temperature  ( tTag, xTag, 2000, 500 )    );
    pID = initFactory.register_expression( new Density      ( rhoTag, tTag, pTag, mmwTag));
    initIDs.insert( tID );
    initIDs.insert( pID );
    initIDs.insert( mID );
    initIDs.insert( yID );
  }

  Expr::ExpressionFactory execFactory;
  std::set<Expr::ExpressionID> execIDs;
  {
    typedef Expr::PlaceHolder   <CellField>::Builder MassFracs;
    typedef Expr::PlaceHolder   <CellField>::Builder Temperature;
    typedef Expr::PlaceHolder   <CellField>::Builder Density;
    typedef Expr::PlaceHolder   <CellField>::Builder MixMolWeight;

    const Expr::ExpressionID yId = execFactory.register_expression( new MassFracs   ( yiTags ) );
    execIDs.insert( yId );
    const Expr::ExpressionID tId = execFactory.register_expression( new Temperature ( tTag   ) );
    execIDs.insert( tId );
    const Expr::ExpressionID dId = execFactory.register_expression( new Density     ( rhoTag ) );
    execIDs.insert( dId );
    const Expr::ExpressionID mId = execFactory.register_expression( new MixMolWeight( mmwTag ) );
    execIDs.insert( mId );
  }

  Expr::ExpressionTree initTree( initIDs, initFactory, 0 );
  Expr::ExpressionTree execTree( execIDs , execFactory, 0 );

  Expr::FieldManagerList fml;
  initTree.register_fields( fml );
  initTree.bind_fields( fml );
  initTree.lock_fields( fml );

  execTree.register_fields( fml );
  execTree.bind_fields( fml );
  execTree.lock_fields( fml );

  SO::IntVec gridSize = SO::IntVec( 20,  1,  1);
  fml.allocate_fields( Expr::FieldAllocInfo( gridSize, 0, 0, false, false, false ) );
  SO::Grid grid( gridSize, SO::DoubleVec(1,1,1) );

  CellField& xcoord = fml.field_ref< CellField >( xTag );
  grid.set_coord<SO::XDIR>( xcoord );
  const int nPoints = xcoord.window_with_ghost().glob_npts();
#   ifdef ENABLE_CUDA
  xcoord.set_device_as_active( GPU_INDEX );
#   endif
  initTree.execute_tree();
  execTree.execute_tree();

  // get fields for aj method

  const so::GhostData         nghost(1);
  const so::BoundaryCellInfo  bcInfo = so::BoundaryCellInfo::build<CellField>( gridSize[0]>1, gridSize[1]>1, gridSize[2]>1 );
  const so::MemoryWindow      window( get_window_with_ghost( gridSize, nghost, bcInfo ) );
  so::FieldVector<CellField> productionRates       ( nSpec+1, window, bcInfo, nghost, LOCATION );
  so::FieldMatrix<CellField> primitiveSensitivities( nSpec+1, window, bcInfo, nghost, LOCATION );

#   ifdef ENABLE_CUDA
  for( size_t prodRateIdx=0; prodRateIdx<nSpec; ++prodRateIdx ){
    productionRates(prodRateIdx).set_device_as_active( GPU_INDEX );
  }
#   endif
#   ifdef ENABLE_CUDA
  for( size_t prodRateIdx=0; prodRateIdx<nSpec+1; ++prodRateIdx ){
    for( size_t primVarIdx=0; primVarIdx<nSpec+1; ++primVarIdx ){
      primitiveSensitivities(prodRateIdx,primVarIdx).set_device_as_active( GPU_INDEX );
    }
  }
#   endif

  CellField T ( window, bcInfo, nghost, NULL, so::InternalStorage, LOCATION );                                // temperature
  so::SpatFldPtr<CellField> pPtr    = so::SpatialFieldStore::get<CellField>( T ); CellField& p = *pPtr;       // pressure
  so::SpatFldPtr<CellField> rhoPtr  = so::SpatialFieldStore::get<CellField>( T ); CellField& rho  = *rhoPtr;  // density
  so::SpatFldPtr<CellField> MmixPtr = so::SpatialFieldStore::get<CellField>( T ); CellField& Mmix = *MmixPtr; // mixture molecular weight
  std::vector< so::SpatFldPtr<CellField> > YPtr;  // mass fractions

  p <<= fml.field_ref<CellField>( pTag );
  T <<= fml.field_ref<CellField>( tTag );
  rho <<= fml.field_ref<CellField>( rhoTag );
  Mmix <<= fml.field_ref<CellField>( mmwTag );

  for( size_t i=0; i<nSpec; ++i ){
    YPtr.push_back( so::SpatialFieldStore::get<CellField>( T ) );
    *YPtr[i] <<= fml.field_ref<CellField>( yiTags[i] );
  }

  ChemicalSourceJacobian csj;

  std::vector< const CellField* > YPtrForJac;
  for( size_t i=0; i<YPtr.size(); ++i )
    YPtrForJac.push_back( &( *YPtr[i] ) );

  csj.evaluate_rates_and_jacobian( primitiveSensitivities, productionRates, T, rho, YPtrForJac, Mmix );

  const std::vector< CellFieldPtrT > canteraResults = get_cantera_results( *gasMix,
                                                                           fml,
                                                                           tTag,
                                                                           yiTags,
                                                                           pTag );
#ifdef ENABLE_CUDA
  for( size_t prodRateIdx=0; prodRateIdx<nSpec; ++prodRateIdx ){
    productionRates( prodRateIdx ).add_device( CPU_INDEX );
  }
#endif

  std::vector< CellFieldPtrT >::const_iterator iCantera = canteraResults.begin();
  size_t i=0;
  BOOST_FOREACH( const Expr::Tag& rTag, rTags ){
    status( field_equal( productionRates(i), **iCantera, 1e-8 ) || field_equal_abs( productionRates(i), **iCantera, 1e-10 ), rTag.name() );
    ++iCantera;
    ++i;
  }

  fml.deallocate_fields();

  return status.ok();
}

//==============================================================================


int main( int iarg, char* carg[] )
{

  std::string inputFileName = "h2-burke.xml"; // mechanism input in xml format
  bool doTimings = false;
  bool doJacVsGs = false;
  bool doProdRates = false;
  double pressure = 101325;
  try {
    po::options_description desc("Supported Options");
    desc.add_options()
      ( "help", "print help message" )
      ( "xml-input-file", po::value<std::string>(&inputFileName)->default_value(inputFileName), "Cantera xml input file name" )
      ( "timings" , "provide timings" )
      ( "jacobian-vs-gs", "test Jacobian vs gold standards" )
      ( "rates-vs-cantera", "test production rates vs Cantera" )
      ( "pressure", po::value<double>(&pressure), "system pressure (for rates-vs-cantera test only)" );

    po::variables_map args;
    po::store( po::parse_command_line(iarg,carg,desc), args );
    po::notify(args);

    doTimings = args.count( "timings" );
    doJacVsGs = args.count( "jacobian-vs-gs" );
    doProdRates = args.count( "rates-vs-cantera" );

    if (args.count("help")) {
      std::cout << desc << "\n";
      return -1;
    }
  }
  catch( std::exception& err ){
    std::cout << "Error parsing input arguments\n" << err.what() << std::endl;
    return -2;
  }



  const CanteraObjects::Setup setup( "Mix", inputFileName );
  CanteraObjects::setup_cantera( setup );
  const size_t ns = CanteraObjects::number_species();
  const size_t nr = CanteraObjects::number_rxns();


  if( doJacVsGs ) generate_jacobian_for_accuracy_test( inputFileName );


  if( doProdRates ){
    TestHelper status;
    status( driver( pressure ), "Production rates vs cantera" );

    if( status.ok() ){
      std::cout << "\nPASS\n";
      return 0;
    }
  }


  if( doTimings ){
    std::cout << std::endl
              << "Timings for prod. rates, Jacobian matrix, state transforms, and linear algebra" << std::endl
              << "------------------------------------------------------------------------------" << std::endl
              << " - mechanism          : " << inputFileName << std::endl
              << " - number of species  : " << ns << std::endl
              << " - number of reactions: " << nr << std::endl << std::endl;

    get_timings( 4   , 512, inputFileName );
    get_timings( 8   , 256, inputFileName );
    get_timings( 16  , 128, inputFileName );
    get_timings( 32  , 64 , inputFileName );
    get_timings( 64  , 32 , inputFileName );
    get_timings( 128 , 16 , inputFileName );
    get_timings( 256 , 8  , inputFileName );
    get_timings( 512 , 4  , inputFileName );
  }


}

