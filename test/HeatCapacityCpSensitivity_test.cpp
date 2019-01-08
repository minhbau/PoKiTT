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
#include <pokitt/thermo/HeatCapacity_Cp.h>

#include "TestHelper.h"

namespace so = SpatialOps;


using CellFieldT = so::SVolField;
using CellFieldPtrT = so::SpatFldPtr<CellFieldT>;

using PlaceHolderT   = Expr::PlaceHolder<CellFieldT>::Builder;
using ConstantT      = Expr::ConstantExpr<CellFieldT>::Builder;
using LinearT        = Expr::LinearFunction<CellFieldT>::Builder;

using CpT = pokitt::HeatCapacity_Cp<CellFieldT>::Builder;
using CpiT = pokitt::SpeciesHeatCapacity_Cp<CellFieldT>::Builder;
using SpeciesNT = pokitt::SpeciesN<CellFieldT>::Builder;


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

  const double offset = 1e-6;
  const double leftTemp = 400;
  const double rightTemp = 3000;
  const double pressure = 101325;

  Expr::ExpressionFactory factory;
  Expr::ExpressionTree tree( factory, 0, "tree" );
  Expr::FieldManagerList fml;

  const Expr::Tag gridTag( "x", Expr::STATE_NONE );
  const Expr::TagList massTags = Expr::tag_list( spNames, Expr::STATE_NONE, "Y_" );
  const Expr::TagList specCpTags = Expr::tag_list( spNames, Expr::STATE_NONE, "Cp_" );
  const Expr::TagList offsetSpecCpTags = Expr::tag_list( spNames, Expr::STATE_NONE, "Cp_offset_" );
  const Expr::Tag CpTag( "Cp", Expr::STATE_NONE );
  const Expr::Tag offsetCpTag( "Cp+dCp", Expr::STATE_NONE );
  const Expr::Tag tempTag( "T", Expr::STATE_NONE );
  const Expr::Tag offsetTempTag( "T+dT", Expr::STATE_NONE );

  tree.insert_tree( factory.register_expression( new PlaceHolderT( gridTag ) ) );
  tree.insert_tree( factory.register_expression( new LinearT( tempTag, gridTag, rightTemp-leftTemp, leftTemp ) ) );
  tree.insert_tree( factory.register_expression( new LinearT( offsetTempTag, gridTag, rightTemp-leftTemp, leftTemp+offset ) ) );
  for( int i=0; i<nSpec-1; ++i ){
    tree.insert_tree( factory.register_expression( new ConstantT( massTags[i], massFracs[i] ) ) );
  }
  tree.insert_tree( factory.register_expression( new SpeciesNT( massTags[nSpec-1], massTags ) ) );
  tree.insert_tree( factory.register_expression( new CpT( CpTag, tempTag, massTags ) ) );
  tree.insert_tree( factory.register_expression( new CpT( offsetCpTag, offsetTempTag, massTags ) ) );
  for( int i=0; i<nSpec; ++i ){
    tree.insert_tree( factory.register_expression( new CpiT( specCpTags[i]      , tempTag      , i ) ) );
    tree.insert_tree( factory.register_expression( new CpiT( offsetSpecCpTags[i], offsetTempTag, i ) ) );
  }

  Expr::TagList sensFxnTags = specCpTags;
  sensFxnTags.push_back( CpTag );

  Expr::TagList sensVarTags;
  for( auto i=0; i<nSpec-1; ++i ) sensVarTags.push_back( massTags[i] );
  sensVarTags.push_back( tempTag );
  tree.compute_sensitivities( sensFxnTags, sensVarTags );

  {
    std::ofstream fout("cp_sens.dot");
    tree.write_tree(fout,false,true);
  }

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
  fml.field_ref<CellFieldT>( specCpTags[i] ).add_device( CPU_INDEX );
  fml.field_ref<CellFieldT>( offsetSpecCpTags[i] ).add_device( CPU_INDEX );
}
#   endif
// need to do add_device on all sensitivity fields but this won't work until GPU sensitivity support is in ExprLib

  CellFieldT& T        = fml.field_ref<CellFieldT>( tempTag );
  CellFieldT& TpdT     = fml.field_ref<CellFieldT>( offsetTempTag );
  CellFieldT& cPMix    = fml.field_ref<CellFieldT>( CpTag );
  CellFieldT& cPpdTMix = fml.field_ref<CellFieldT>( offsetCpTag );

  TestHelper fullTest( true );

  TestHelper species( false );
  for( int i=0; i<nSpec; ++i ){
    CellFieldT& cP = fml.field_ref<CellFieldT>( specCpTags[i] );
    CellFieldT& cPpdT = fml.field_ref<CellFieldT>( offsetSpecCpTags[i] );
    CellFieldPtrT dCPdTPtr = so::SpatialFieldStore::get<CellFieldT>( T );
    *dCPdTPtr <<= ( cPpdT - cP ) / ( TpdT - T );

    const Expr::Tag sensTag = Expr::sens_tag( specCpTags[i], tempTag );
    species( so::field_equal( *dCPdTPtr, fml.field_ref<CellFieldT>( sensTag ), 1e-2 ), sensTag.name() );
  }
  fullTest( species.ok(), "species heat capacities" );

  TestHelper mixture( false );
  CellFieldPtrT dHdTPtr = so::SpatialFieldStore::get<CellFieldT>( T );
  *dHdTPtr <<= ( cPpdTMix - cPMix ) / ( TpdT - T );
  mixture( so::field_equal( *dHdTPtr, fml.field_ref<CellFieldT>( Expr::sens_tag( CpTag, tempTag ) ), 1e-4 ), "Cp_sens_T" );


  { // check some intermediate sensitivities
    TestHelper status(false);
    for( auto i=0; i<nSpec-1; ++i ){
      const Expr::Tag& yi = massTags[i];
      status( SpatialOps::field_equal(  1.0, fml.field_ref<CellFieldT>( Expr::sens_tag( yi, yi ) ) ), "d" + yi.name() + "/d" + yi.name() );
      status( SpatialOps::field_equal( -1.0, fml.field_ref<CellFieldT>( Expr::sens_tag( massTags[nSpec-1], yi ) ) ), "d" + massTags[nSpec-1].name() + "/d" + yi.name() );
    }
    fullTest( status.ok(), "computed fields" );
  }

  const CellFieldT& cpn = fml.field_ref<CellFieldT>( specCpTags[nSpec-1] );
  for( int i=0; i<nSpec-1; ++i ){
    CellFieldPtrT dCPdYi = so::SpatialFieldStore::get<CellFieldT>( T );
    *dCPdYi <<= fml.field_ref<CellFieldT>( specCpTags[i] ) - cpn;

    const Expr::Tag sensTag = Expr::sens_tag( CpTag, massTags[i] );
    mixture( so::field_equal( *dCPdYi, fml.field_ref<CellFieldT>( sensTag ), 1e-8 ), sensTag.name() );
//    std::cout << "\t Expected: " << (*dCPdYi)[0] << ", found: " << fml.field_ref<CellFieldT>( sensTag )[0] << std::endl;
  }
  fullTest( mixture.ok(), "mixture heat capacity" );

  if( fullTest.ok() ){
    std::cout << "\nPASS\n";
    return 0;
  }
  std::cout << "\nFAIL\n";
  return -1;
}
