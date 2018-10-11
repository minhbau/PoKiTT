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
#include <pokitt/thermo/HeatCapacity_Cp.h>
#include <pokitt/thermo/Enthalpy.h>

#include "TestHelper.h"

namespace so = SpatialOps;


using CellFieldT = so::SVolField;
using CellFieldPtrT = so::SpatFldPtr<CellFieldT>;

using PlaceHolderT   = Expr::PlaceHolder<CellFieldT>::Builder;
using ConstantT      = Expr::ConstantExpr<CellFieldT>::Builder;
using LinearT        = Expr::LinearFunction<CellFieldT>::Builder;

using SpeciesNT = pokitt::SpeciesN<CellFieldT>::Builder;
using TemperatureT = pokitt::Temperature<CellFieldT>::Builder;
using HeatCapacity_CpT = pokitt::HeatCapacity_Cp<CellFieldT>::Builder;
using SpecEnthT = pokitt::SpeciesEnthalpy<CellFieldT>::Builder;


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
  std::string inputFileName = "h2-burke.xml";

  const CanteraObjects::Setup setup( "Mix", inputFileName );
  CanteraObjects::setup_cantera( setup );

  const size_t nSpec = CanteraObjects::number_species();
  const std::vector<std::string>& spNames = CanteraObjects::species_names();

  std::vector<double> massFracs( nSpec, 0.0 );
  for( int i=0; i<nSpec; ++i ){
    massFracs[i] = 1.0 / static_cast<double>( nSpec );
  }

  const double offset = 1;
  const double leftEnth = 2e7;
  const double rightEnth = 3e7;
  const double pressure = 101325;

  Expr::ExpressionFactory factory;
  Expr::ExpressionTree tree( factory, 0, "tree" );
  Expr::FieldManagerList fml;

  const Expr::Tag gridTag( "x", Expr::STATE_NONE );
  const Expr::TagList massTags = Expr::tag_list( spNames, Expr::STATE_NONE, "Y_" );
  const Expr::TagList specEnthTags = Expr::tag_list( spNames, Expr::STATE_NONE, "h_" );
  const Expr::Tag enthTag( "h", Expr::STATE_NONE );
  const Expr::Tag offsetEnthTag( "h+dh", Expr::STATE_NONE );
  const Expr::Tag tempTag( "T", Expr::STATE_NONE );
  const Expr::Tag offsetTempTag( "T+dT", Expr::STATE_NONE );
  const Expr::Tag cpTag( "cp", Expr::STATE_NONE );
  const Expr::Tag tempGuessTag( "T_guess", Expr::STATE_NONE );

  tree.insert_tree( factory.register_expression( new ConstantT( tempGuessTag, 1200 ) ) );
  tree.insert_tree( factory.register_expression( new PlaceHolderT( gridTag ) ) );
  tree.insert_tree( factory.register_expression( new LinearT( enthTag, gridTag, rightEnth-leftEnth, leftEnth ) ) );
  tree.insert_tree( factory.register_expression( new LinearT( offsetEnthTag, gridTag, rightEnth-leftEnth, leftEnth+offset ) ) );

  for( int i=0; i<nSpec-1; ++i ){
    tree.insert_tree( factory.register_expression( new ConstantT( massTags[i], massFracs[i] ) ) );
  }
  tree.insert_tree( factory.register_expression( new SpeciesNT( massTags[nSpec-1], massTags ) ) );
  tree.insert_tree( factory.register_expression( new TemperatureT( tempTag, massTags, enthTag, tempGuessTag, 1e-3, 5000, 20 )));
  tree.insert_tree( factory.register_expression( new TemperatureT( offsetTempTag, massTags, offsetEnthTag, tempGuessTag, 1e-3, 5000, 20 )));
  tree.insert_tree( factory.register_expression( new HeatCapacity_CpT( cpTag, tempTag, massTags ) ) );
  for( int i=0; i<nSpec; ++i ){
    tree.insert_tree( factory.register_expression( new SpecEnthT( specEnthTags[i], tempTag, i ) ) );
  }

  {
    std::ofstream tGraph( "temperature_sens.dot" );
    tree.write_tree( tGraph,false,true );
  }

  Expr::TagList sensFxnTags = tag_list( tempTag );
  Expr::TagList sensVarTags = tag_list( enthTag);
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
  CellFieldT& hMix = fml.field_ref<CellFieldT>( enthTag );
  CellFieldT& hpdTMix = fml.field_ref<CellFieldT>( offsetEnthTag );
  CellFieldT& cp = fml.field_ref<CellFieldT>( cpTag );

  TestHelper fullTest( true );

  TestHelper mixture( false );
  CellFieldPtrT dTdHPtr = so::SpatialFieldStore::get<CellFieldT>( T );
  *dTdHPtr <<= ( TpdT - T ) / ( hpdTMix - hMix );
  mixture( so::field_equal( *dTdHPtr, fml.field_ref<CellFieldT>( Expr::sens_tag( tempTag, enthTag ) ), 1e-4 ), "T_sens_h" );

  for( int i=0; i<nSpec-1; ++i ){
    CellFieldPtrT dTdYi = so::SpatialFieldStore::get<CellFieldT>( T );
    *dTdYi <<= -(fml.field_ref<CellFieldT>( specEnthTags[i] ) - fml.field_ref<CellFieldT>( specEnthTags[nSpec-1] ))/cp;

    const Expr::Tag sensTag = Expr::sens_tag( tempTag, massTags[i] );

    mixture( so::field_equal( *dTdYi, fml.field_ref<CellFieldT>( sensTag ), 1e-5 ), sensTag.name() );
  }
  fullTest( mixture.ok(), "mixture enthalpy" );


  if( fullTest.ok() ){
    std::cout << "\nPASS\n";
    return 0;
  }
  else{
    std::cout << "\nFAIL\n";
    return -1;
  }

}





