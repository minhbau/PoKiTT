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
#include <pokitt/MixtureMolWeight.h>
#include <pokitt/SpeciesN.h>
#include <pokitt/transport/ThermalCondMix.h>
#include "TestHelper.h"

namespace so = SpatialOps;

/*
 * Molar weight of mixture is a function of `Y_i`.
 */


using CellFieldT = so::SVolField;
using CellFieldPtrT = so::SpatFldPtr<CellFieldT>;

using PlaceHolderT   = Expr::PlaceHolder<CellFieldT>::Builder;
using ConstantT      = Expr::ConstantExpr<CellFieldT>::Builder;
using LinearT        = Expr::LinearFunction<CellFieldT>::Builder;

using MixtureMWT   = pokitt::MixtureMolWeight<CellFieldT>::Builder;
using SpeciesNT    = pokitt::SpeciesN<CellFieldT>::Builder;
using ThermalCondT = pokitt::ThermalConductivity<CellFieldT>::Builder;


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

  const std::vector<double> mwVec = CanteraObjects::molecular_weights();

  const double offset = 1e-6;

  std::vector<double> massFracs( nSpec, 0.0 );
  std::vector<double> massFracsOffset( nSpec, 0.0 );
  std::vector<double> moleFracs( nSpec, 0.0 );
  std::vector<double> moleFracsOffset( nSpec, 0.0 );
  for( int i=0; i<nSpec; ++i ){
    massFracs[i] = 1.0 / static_cast<double>( nSpec );
    moleFracs[i] = 1.0 / static_cast<double>( nSpec );
  }
  for( int i=0; i<nSpec-1; ++i ){
    massFracsOffset[i] = massFracs[i] + offset;
    moleFracsOffset[i] = moleFracs[i] + offset;
  }

  TestHelper fullTest( true );

  /* check dMw/dYi */
  for( int i = 0;i < nSpec - 1;++i )
  {

    Expr::ExpressionFactory factory;
    Expr::ExpressionTree tree( factory, 0, "tree" );
    Expr::FieldManagerList fml;

    const Expr::Tag gridTag( "x", Expr::STATE_NONE );
    const Expr::TagList massTags = Expr::tag_list( spNames, Expr::STATE_NONE, "Y_" );
    const Expr::TagList massOffsetTags = Expr::tag_list( spNames, Expr::STATE_NONE, "Y_", "_offset" );
    const Expr::Tag mmwTag( "mmw", Expr::STATE_NONE );
    const Expr::Tag mmwOffsetTag( "mmw_offset", Expr::STATE_NONE );

    tree.insert_tree( factory.register_expression( new PlaceHolderT( gridTag )));
    tree.insert_tree( factory.register_expression( new ConstantT( massTags[i], massFracs[i] )));
    tree.insert_tree( factory.register_expression( new ConstantT( massOffsetTags[i], massFracsOffset[i] )));
    for ( int j = 0; j < nSpec - 1;++j){
      if( i != j ){
        tree.insert_tree( factory.register_expression( new ConstantT( massTags[j], massFracs[j] )));
        tree.insert_tree( factory.register_expression( new ConstantT( massOffsetTags[j], massFracs[j] )));
      }
    }
    tree.insert_tree( factory.register_expression( new SpeciesNT( massTags[nSpec - 1], massTags )));
    tree.insert_tree( factory.register_expression( new SpeciesNT( massOffsetTags[nSpec - 1], massOffsetTags )));
    tree.insert_tree( factory.register_expression( new MixtureMWT( mmwTag, massTags, pokitt::MASS )));
    tree.insert_tree( factory.register_expression( new MixtureMWT( mmwOffsetTag, massOffsetTags, pokitt::MASS )));

    Expr::TagList sensVarTags = tag_list( massTags[i] );
    Expr::TagList sensFxnTags = tag_list( mmwTag );

    tree.compute_sensitivities( sensFxnTags, sensVarTags );

    tree.register_fields( fml );
    tree.bind_fields( fml );
    tree.lock_fields( fml );

    so::IntVec gridSize = so::IntVec( 20, 1, 1 );
    fml.allocate_fields( Expr::FieldAllocInfo( gridSize, 0, 0, false, false, false ));
    so::Grid grid( gridSize, so::DoubleVec( 1, 1, 1 ));

    CellFieldT& xcoord = fml.field_ref<CellFieldT>( gridTag );
    grid.set_coord<so::XDIR>( xcoord );
    const int nPoints = xcoord.window_with_ghost().glob_npts();
#   ifdef ENABLE_CUDA
    xcoord.set_device_as_active( GPU_INDEX );
#   endif
    tree.execute_tree();

#   ifdef ENABLE_CUDA
    fml.field_ref<CellFieldT>( gridTag ).add_device( CPU_INDEX );
      fml.field_ref<CellFieldT>( mmwTag  ).add_device( CPU_INDEX );
      fml.field_ref<CellFieldT>( mmwOffsetTag  ).add_device( CPU_INDEX );
    for( int i=0; i<nSpec; ++i ){
      fml.field_ref<CellFieldT>( massTags[i] ).add_device( CPU_INDEX );
      fml.field_ref<CellFieldT>( massOffsetTags[i] ).add_device( CPU_INDEX );
    }
#   endif
    // need to do add_device on all sensitivity fields but this won't work until GPU sensitivity support is in ExprLib

    TestHelper MassFractionSens( true );

    CellFieldT& mmwField = fml.field_ref<CellFieldT>( mmwTag );
    CellFieldT& mmwOffsetField = fml.field_ref<CellFieldT>( mmwOffsetTag );

    CellFieldT& massFracField = fml.field_ref<CellFieldT>( massTags[i] );
    CellFieldT& massFracOffsetField = fml.field_ref<CellFieldT>( massOffsetTags[i] );
    CellFieldPtrT dmmwdYi = so::SpatialFieldStore::get<CellFieldT>( massFracField );
    *dmmwdYi <<= ( mmwOffsetField - mmwField ) / ( massFracOffsetField - massFracField);

    const std::string sensStr = mmwTag.name() + "_sens_" + massTags[i].name();
    const Expr::Tag sensTag( sensStr, Expr::STATE_NONE );
    MassFractionSens( so::field_equal( *dmmwdYi, fml.field_ref<CellFieldT>( sensTag ), 1e-3 ), sensStr );
    fullTest( MassFractionSens.ok(), "mixture molar weight sensitivity to species mass fraction" );
  }

  /* check dMw/dXi */
  for( int i = 0;i < nSpec - 1;++i )
  {

    Expr::ExpressionFactory factory;
    Expr::ExpressionTree tree( factory, 0, "tree" );
    Expr::FieldManagerList fml;

    const Expr::Tag gridTag( "x", Expr::STATE_NONE );
    const Expr::TagList moleTags = Expr::tag_list( spNames, Expr::STATE_NONE, "X_" );
    const Expr::TagList moleOffsetTags = Expr::tag_list( spNames, Expr::STATE_NONE, "X_", "_offset" );
    const Expr::Tag mmwTag( "mmw", Expr::STATE_NONE );
    const Expr::Tag mmwOffsetTag( "mmw_offset", Expr::STATE_NONE );

    tree.insert_tree( factory.register_expression( new PlaceHolderT( gridTag )));
    tree.insert_tree( factory.register_expression( new ConstantT( moleTags[i], moleFracs[i] )));
    tree.insert_tree( factory.register_expression( new ConstantT( moleOffsetTags[i], moleFracsOffset[i] )));
    for ( int j = 0; j < nSpec - 1;++j){
      if( i != j ){
        tree.insert_tree( factory.register_expression( new ConstantT( moleTags[j], moleFracs[j] )));
        tree.insert_tree( factory.register_expression( new ConstantT( moleOffsetTags[j], moleFracs[j] )));
      }
    }
    tree.insert_tree( factory.register_expression( new SpeciesNT( moleTags[nSpec - 1], moleTags )));
    tree.insert_tree( factory.register_expression( new SpeciesNT( moleOffsetTags[nSpec - 1], moleOffsetTags )));
    tree.insert_tree( factory.register_expression( new MixtureMWT( mmwTag, moleTags, pokitt::MOLE )));
    tree.insert_tree( factory.register_expression( new MixtureMWT( mmwOffsetTag, moleOffsetTags, pokitt::MOLE )));

    Expr::TagList sensVarTags = tag_list( moleTags[i] );
    Expr::TagList sensFxnTags = tag_list( mmwTag );

    tree.compute_sensitivities( sensFxnTags, sensVarTags );

    tree.register_fields( fml );
    tree.bind_fields( fml );
    tree.lock_fields( fml );

    so::IntVec gridSize = so::IntVec( 20, 1, 1 );
    fml.allocate_fields( Expr::FieldAllocInfo( gridSize, 0, 0, false, false, false ));
    so::Grid grid( gridSize, so::DoubleVec( 1, 1, 1 ));

    CellFieldT& xcoord = fml.field_ref<CellFieldT>( gridTag );
    grid.set_coord<so::XDIR>( xcoord );
    const int nPoints = xcoord.window_with_ghost().glob_npts();
#   ifdef ENABLE_CUDA
    xcoord.set_device_as_active( GPU_INDEX );
#   endif
    tree.execute_tree();

#   ifdef ENABLE_CUDA
    fml.field_ref<CellFieldT>( gridTag ).add_device( CPU_INDEX );
      fml.field_ref<CellFieldT>( mmwTag  ).add_device( CPU_INDEX );
      fml.field_ref<CellFieldT>( mmwOffsetTag  ).add_device( CPU_INDEX );
    for( int i=0; i<nSpec; ++i ){
      fml.field_ref<CellFieldT>( moleTags[i] ).add_device( CPU_INDEX );
      fml.field_ref<CellFieldT>( moleOffsetTags[i] ).add_device( CPU_INDEX );
    }
#   endif
    // need to do add_device on all sensitivity fields but this won't work until GPU sensitivity support is in ExprLib

    TestHelper MoleFractionSens( true );

    CellFieldT& mmwField = fml.field_ref<CellFieldT>( mmwTag );
    CellFieldT& mmwOffsetField = fml.field_ref<CellFieldT>( mmwOffsetTag );

    CellFieldT& moleFracField = fml.field_ref<CellFieldT>( moleTags[i] );
    CellFieldT& moleFracOffsetField = fml.field_ref<CellFieldT>( moleOffsetTags[i] );
    CellFieldPtrT dmmwdYi = so::SpatialFieldStore::get<CellFieldT>( moleFracField );
    *dmmwdYi <<= ( mmwOffsetField - mmwField ) / ( moleFracOffsetField - moleFracField);

    const std::string sensStr = mmwTag.name() + "_sens_" + moleTags[i].name();
    const Expr::Tag sensTag( sensStr, Expr::STATE_NONE );
    MoleFractionSens( so::field_equal( *dmmwdYi, fml.field_ref<CellFieldT>( sensTag ), 1e-3 ), sensStr );
    fullTest( MoleFractionSens.ok(), "mixture molar weight sensitivity to species mole fraction" );
  }

  if( fullTest.ok() ){
    std::cout << "\nPASS\n";
    return 0;
  }
  else{
    std::cout << "\nFAIL\n";
    return -1;
  }

}
