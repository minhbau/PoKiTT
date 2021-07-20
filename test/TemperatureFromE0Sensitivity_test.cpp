/*
 * The MIT License
 *
 * Copyright (c) 2016-2017 The University of Utah
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

#include <spatialops/Nebo.h>
#include <spatialops/structured/FieldComparisons.h>
#include <spatialops/structured/Grid.h>

#include <expression/ExprLib.h>

#include <pokitt/CanteraObjects.h>
#include <pokitt/SpeciesN.h>
#include <pokitt/thermo/Temperature.h>
#include <pokitt/thermo/HeatCapacity_Cv.h>
#include <pokitt/thermo/InternalEnergy.h>

#include "TestHelper.h"

namespace so = SpatialOps;


using CellFieldT = so::SVolField;
using CellFieldPtrT = so::SpatFldPtr<CellFieldT>;

using PlaceHolderT   = Expr::PlaceHolder<CellFieldT>::Builder;
using ConstantT      = Expr::ConstantExpr<CellFieldT>::Builder;
using LinearT        = Expr::LinearFunction<CellFieldT>::Builder;

using SpeciesNT = pokitt::SpeciesN<CellFieldT>::Builder;
using TemperatureFromE0T = pokitt::TemperatureFromE0<CellFieldT>::Builder;
using HeatCapacity_CvT = pokitt::HeatCapacity_Cv<CellFieldT>::Builder;
using SpecEgyT = pokitt::SpeciesInternalEnergy<CellFieldT>::Builder;


/**
 * !!!! GPU SUPPORT WILL BE DISABLED UNTIL SENSITIVITIES WORK ON GPU !!!!
 */
#ifdef ENABLE_CUDA
# define LOCATION GPU_INDEX
#else
# define LOCATION CPU_INDEX
#endif
/**
 * !!!! GPU SUPPORT WILL BE DISABLED UNTIL SENSITIVITIES WORK ON GPU !!!!
 */


int main()
{
  std::string inputFileName = "h2-burke.yaml";

  const CanteraObjects::Setup setup( "Mix", inputFileName );
  CanteraObjects::setup_cantera( setup );

  const size_t nSpec = CanteraObjects::number_species();
  const std::vector<std::string>& spNames = CanteraObjects::species_names();

  std::vector<double> massFracs( nSpec, 0.0 );
  for( int i=0; i<nSpec; ++i ){
    massFracs[i] = 1.0 / static_cast<double>( nSpec );
  }

  const double offset = 1;
  const double pressure = 101325;
  const double leftEgy = 2e7;
  const double rightEgy = 3e7;

  Expr::ExpressionFactory factory;
  Expr::ExpressionTree tree( factory, 0, "tree" );
  Expr::FieldManagerList fml;

  const Expr::Tag gridTag( "x", Expr::STATE_NONE );
  const Expr::TagList massTags = Expr::tag_list( spNames, Expr::STATE_NONE, "Y_" );
  const Expr::TagList specEgyTags = Expr::tag_list( spNames, Expr::STATE_NONE, "e_" );
  const Expr::Tag egyTag( "e", Expr::STATE_NONE );
  const Expr::Tag offsetEgyTag( "e+de", Expr::STATE_NONE );
  const Expr::Tag tempTag( "T", Expr::STATE_NONE );
  const Expr::Tag offsetTempTag( "T+dT", Expr::STATE_NONE );
  const Expr::Tag cvTag( "cv", Expr::STATE_NONE );
  const Expr::Tag keTag( "ke", Expr::STATE_NONE );
  const Expr::Tag tempGuessTag( "T_guess", Expr::STATE_NONE );

  tree.insert_tree( factory.register_expression( new ConstantT( tempGuessTag, 1200 ) ) );
  tree.insert_tree( factory.register_expression( new PlaceHolderT( gridTag ) ) );
  tree.insert_tree( factory.register_expression( new LinearT( egyTag, gridTag, rightEgy-leftEgy, leftEgy ) ) );
  tree.insert_tree( factory.register_expression( new LinearT( offsetEgyTag, gridTag, rightEgy-leftEgy, leftEgy+offset ) ) );

  for( int i=0; i<nSpec-1; ++i ){
    tree.insert_tree( factory.register_expression( new ConstantT( massTags[i], massFracs[i] ) ) );
  }
  tree.insert_tree( factory.register_expression( new SpeciesNT( massTags[nSpec-1], massTags ) ) );

  tree.insert_tree( factory.register_expression( new ConstantT( keTag, 0.0 ) ) );
  tree.insert_tree( factory.register_expression( new TemperatureFromE0T( tempTag, massTags, egyTag, keTag, tempGuessTag, 1e-3, 5000, 20 )));
  tree.insert_tree( factory.register_expression( new TemperatureFromE0T( offsetTempTag, massTags, offsetEgyTag, keTag, tempGuessTag, 1e-3, 5000, 20 )));
  tree.insert_tree( factory.register_expression( new HeatCapacity_CvT( cvTag, tempTag, massTags ) ) );
  for( int i=0; i<nSpec; ++i ){
    tree.insert_tree( factory.register_expression( new SpecEgyT( specEgyTags[i], tempTag, i ) ) );
  }

  Expr::TagList sensFxnTags = tag_list( tempTag );
  Expr::TagList sensVarTags = tag_list( egyTag );
  for (int i=0; i<nSpec-1; ++i){
    sensVarTags.push_back( massTags[i] );
  }

  tree.compute_sensitivities( sensFxnTags, sensVarTags );

  tree.register_fields( fml );
  tree.bind_fields( fml );
  tree.lock_fields( fml );

  so::IntVec gridSize = so::IntVec( 100,  1,  1 );
  fml.allocate_fields( Expr::FieldAllocInfo( gridSize, 0, 0, false, false, false ) );
  so::Grid grid( gridSize, so::DoubleVec( 1, 1, 1 ) );

  CellFieldT& xcoord = fml.field_ref<CellFieldT>( gridTag );
  grid.set_coord<so::XDIR>( xcoord );
  const int nPoints = xcoord.window_with_ghost().glob_npts();
#   ifdef ENABLE_CUDA
xcoord.set_device_as_active( GPU_INDEX );
#   endif
  tree.execute_tree();

#   ifdef ENABLE_CUDA
fml.field_ref<CellFieldT>( gridTag ).add_device( CPU_INDEX );
fml.field_ref<CellFieldT>( tempTag ).add_device( CPU_INDEX );
fml.field_ref<CellFieldT>( offsetTempTag ).add_device( CPU_INDEX );
for( int i=0; i<nSpec; ++i ){
  fml.field_ref<CellFieldT>( specEnthTags[i] ).add_device( CPU_INDEX );
}
#   endif
// need to do add_device on all sensitivity fields but this won't work until GPU sensitivity support is in ExprLib

  CellFieldT& T = fml.field_ref<CellFieldT>( tempTag );
  CellFieldT& TpdT = fml.field_ref<CellFieldT>( offsetTempTag );
  CellFieldT& eMix = fml.field_ref<CellFieldT>( egyTag );
  CellFieldT& epdTMix = fml.field_ref<CellFieldT>( offsetEgyTag );
  CellFieldT& cv = fml.field_ref<CellFieldT>( cvTag );

  TestHelper fullTest( true );

  TestHelper mixture( false );
  CellFieldPtrT dTdePtr = so::SpatialFieldStore::get<CellFieldT>( T );
  *dTdePtr <<= ( TpdT - T ) / ( epdTMix -eMix );
  mixture( so::field_equal( *dTdePtr, fml.field_ref<CellFieldT>( Expr::sens_tag( tempTag, egyTag ) ), 1e-4 ), "T_sens_e" );

  for( int i=0; i<nSpec-1; ++i ){
    CellFieldPtrT dTdYi = so::SpatialFieldStore::get<CellFieldT>( T );
    *dTdYi <<= 0.0;

    const Expr::Tag sensTag = Expr::sens_tag( tempTag, massTags[i] );
    mixture( so::field_equal( *dTdYi, fml.field_ref<CellFieldT>( sensTag ), 1e-5 ), sensTag.name() );
  }
  fullTest( mixture.ok(), "mixture enthalpy" );

  if( fullTest.ok() ){
    std::cout << "\nPASS\n";
    return 0;
  }

  std::cout << "\nFAIL\n";
  return -1;
}
