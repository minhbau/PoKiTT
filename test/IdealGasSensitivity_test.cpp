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
#include <pokitt/thermo/Density.h>
#include <pokitt/thermo/Pressure.h>

#include "TestHelper.h"

namespace so = SpatialOps;


using CellFieldT = so::SVolField;
using CellFieldPtrT = so::SpatFldPtr<CellFieldT>;

using PlaceHolderT   = Expr::PlaceHolder<CellFieldT>::Builder;
using ConstantT      = Expr::ConstantExpr<CellFieldT>::Builder;
using LinearT        = Expr::LinearFunction<CellFieldT>::Builder;

using MixtureMWT = pokitt::MixtureMolWeight<CellFieldT>::Builder;
using DensityT   = pokitt::Density<CellFieldT>::Builder;
using PressureT  = pokitt::Pressure<CellFieldT>::Builder;
using SpeciesNT  = pokitt::SpeciesN<CellFieldT>::Builder;


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

  std::vector<double> massFracs( nSpec, 0.0 );
  for( int i=0; i<nSpec; ++i ){
    massFracs[i] = 1.0 / static_cast<double>( nSpec );
  }


  const double leftTemp = 300;
  const double rightTemp = 2000;


  TestHelper allVarSets( true );

  /* set pressure */
  {
    const double leftPres = 10132.5;
    const double rightPres = 1013250;

    Expr::ExpressionFactory factory;
    Expr::ExpressionTree tree( factory, 0, "tree" );
    Expr::FieldManagerList fml;

    const Expr::Tag gridTag( "x", Expr::STATE_NONE );
    const Expr::TagList massTags = Expr::tag_list( spNames, Expr::STATE_NONE, "Y_" );
    const Expr::Tag tempTag( "T", Expr::STATE_NONE );
    const Expr::Tag presTag( "p", Expr::STATE_NONE );
    const Expr::Tag rhoTag( "rho", Expr::STATE_NONE );
    const Expr::Tag mmwTag( "mmw", Expr::STATE_NONE );

    tree.insert_tree( factory.register_expression( new PlaceHolderT( gridTag ) ) );
    tree.insert_tree( factory.register_expression( new LinearT( tempTag, gridTag, rightTemp-leftTemp, leftTemp ) ) );
    tree.insert_tree( factory.register_expression( new LinearT( presTag, gridTag, rightPres-leftPres, leftPres ) ) );
    for( int i=0; i<nSpec-1; ++i ){
      tree.insert_tree( factory.register_expression( new ConstantT( massTags[i], massFracs[i] ) ) );
    }
    tree.insert_tree( factory.register_expression( new SpeciesNT( massTags[nSpec-1], massTags ) ) );
    tree.insert_tree( factory.register_expression( new MixtureMWT( mmwTag, massTags, pokitt::MASS ) ) );
    tree.insert_tree( factory.register_expression( new DensityT( rhoTag, tempTag, presTag, mmwTag ) ) );

    Expr::TagList sensVarTags = {tempTag, presTag};
    Expr::TagList sensFxnTags = sensVarTags;
    for( int i=0; i<nSpec-1; ++i ){
      sensFxnTags.push_back( massTags[i] );
      sensVarTags.push_back( massTags[i] );
    }
    Expr::TagList computedTags = {rhoTag, mmwTag, massTags[nSpec-1]};
    for( const auto& t : computedTags ){
      sensFxnTags.push_back( t );
    }

    tree.compute_sensitivities( sensFxnTags, sensVarTags );

    tree.register_fields( fml );
    tree.bind_fields( fml );
    tree.lock_fields( fml );

    so::IntVec gridSize = so::IntVec( 20,  1,  1 );
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
    fml.field_ref<CellFieldT>( presTag ).add_device( CPU_INDEX );
    fml.field_ref<CellFieldT>( rhoTag  ).add_device( CPU_INDEX );
    fml.field_ref<CellFieldT>( mmwTag  ).add_device( CPU_INDEX );
    for( int i=0; i<nSpec; ++i ){
      fml.field_ref<CellFieldT>( massTags[i] ).add_device( CPU_INDEX );
    }
#   endif
    // need to do add_device on all sensitivity fields but this won't work until GPU sensitivity support is in ExprLib

    TestHelper allSens( true );

    TestHelper varToVar( false );

    for( const auto& f : sensVarTags ){
      for( const auto& v : sensVarTags ){
        const Expr::Tag ddTag = sens_tag( f, v );

        if( f == v ){
          varToVar( so::field_equal( 1.0, fml.field_ref<CellFieldT>( ddTag ) ), ddTag.name() + " = 1.0" );
        }
        else{
          varToVar( so::field_equal( 0.0, fml.field_ref<CellFieldT>( ddTag ) ), ddTag.name() + " = 0.0" );
        }
      }
    }
    allSens( varToVar.ok(), "variable-to-variable identity matrix" );


    CellFieldT& rhoField = fml.field_ref<CellFieldT>( rhoTag );
    CellFieldT& pField   = fml.field_ref<CellFieldT>( presTag );
    CellFieldT& TField   = fml.field_ref<CellFieldT>( tempTag );
    CellFieldT& mmwField = fml.field_ref<CellFieldT>( mmwTag );
    CellFieldT& YnsField = fml.field_ref<CellFieldT>( massTags[nSpec-1] );
    CellFieldPtrT drhodp = so::SpatialFieldStore::get<CellFieldT>( rhoField );
    CellFieldPtrT drhodT = so::SpatialFieldStore::get<CellFieldT>( rhoField );

    TestHelper nthSpecies( false );
    const Expr::Tag& yNsTag = massTags[nSpec-1];
    const std::string& yNsName = yNsTag.name();
    for( int i=0; i<nSpec-1; ++i ){
      const Expr::Tag yNsSensTag = Expr::sens_tag( yNsTag, massTags[i] );

      CellFieldPtrT dYnsdYi = so::SpatialFieldStore::get<CellFieldT>( rhoField );
      *dYnsdYi <<= -1.0;

      nthSpecies( so::field_equal( *dYnsdYi, fml.field_ref<CellFieldT>( yNsSensTag ), 1e-8 ), yNsSensTag.name() );
    }
    nthSpecies( so::field_equal( 0.0, fml.field_ref<CellFieldT>( Expr::sens_tag( yNsTag, presTag ) ), 1e-8 ), yNsName + "__sens_p" );
    nthSpecies( so::field_equal( 0.0, fml.field_ref<CellFieldT>( Expr::sens_tag( yNsTag, tempTag ) ), 1e-8 ), yNsName + "_sens_T" );
    allSens( nthSpecies.ok(), "n-th species" );


    TestHelper mixMolWeight( false );
    for( int i=0; i<nSpec-1; ++i ){
      const Expr::Tag mmwSensTag = Expr::sens_tag( mmwTag, massTags[i] );

      CellFieldPtrT dmmwdYi = so::SpatialFieldStore::get<CellFieldT>( rhoField );
      *dmmwdYi <<= - ( 1.0 / mwVec[i] - 1.0 / mwVec[nSpec-1] ) * mmwField * mmwField;

      mixMolWeight( so::field_equal( *dmmwdYi, fml.field_ref<CellFieldT>( mmwSensTag ), 1e-8 ), mmwSensTag.name() );
    }
    mixMolWeight( so::field_equal( 0.0, fml.field_ref<CellFieldT>( Expr::sens_tag( mmwTag, presTag ) ), 1e-8 ), "mmw_sens_p" );
    mixMolWeight( so::field_equal( 0.0, fml.field_ref<CellFieldT>( Expr::sens_tag( mmwTag, tempTag ) ), 1e-8 ), "mmw_sens_T" );
    allSens( mixMolWeight.ok(), "mixture molecular weight" );


    TestHelper density( false );
    for( int i=0; i<nSpec-1; ++i ){
      const Expr::Tag rhoSensTag = Expr::sens_tag( rhoTag, massTags[i] );
      const Expr::Tag mmwSensTag = Expr::sens_tag( mmwTag, massTags[i] );

      CellFieldPtrT drhodYi = so::SpatialFieldStore::get<CellFieldT>( rhoField );
      *drhodYi <<= rhoField / mmwField * fml.field_ref<CellFieldT>( mmwSensTag );

      density( so::field_equal( *drhodYi, fml.field_ref<CellFieldT>( rhoSensTag ), 1e-8 ), rhoSensTag.name() );
    }
    *drhodp <<= rhoField / pField;
    *drhodT <<= - rhoField / TField;
    density( so::field_equal( *drhodp, fml.field_ref<CellFieldT>( Expr::sens_tag( rhoTag, presTag ) ), 1e-8 ), "rho_sens_p" );
    density( so::field_equal( *drhodT, fml.field_ref<CellFieldT>( Expr::sens_tag( rhoTag, tempTag ) ), 1e-8 ), "rho_sens_T" );
    allSens( density.ok(), "density" );
    allVarSets( allSens.ok(), "set pressure" );
  }

  /* set density */
  {
    const double leftRho = 0.1;
    const double rightRho = 10.0;

    Expr::ExpressionFactory factory;
    Expr::ExpressionTree tree( factory, 0, "tree" );
    Expr::FieldManagerList fml;

    const Expr::Tag gridTag( "x", Expr::STATE_NONE );
    const Expr::TagList massTags = Expr::tag_list( spNames, Expr::STATE_NONE, "Y_" );
    const Expr::Tag tempTag( "T", Expr::STATE_NONE );
    const Expr::Tag presTag( "p", Expr::STATE_NONE );
    const Expr::Tag rhoTag( "rho", Expr::STATE_NONE );
    const Expr::Tag mmwTag( "mmw", Expr::STATE_NONE );

    tree.insert_tree( factory.register_expression( new PlaceHolderT( gridTag ) ) );
    tree.insert_tree( factory.register_expression( new LinearT( tempTag, gridTag, leftTemp, rightTemp-leftTemp ) ) );
    tree.insert_tree( factory.register_expression( new LinearT( rhoTag, gridTag, leftRho, rightRho-leftRho ) ) );
    for( int i=0; i<nSpec-1; ++i ){
      tree.insert_tree( factory.register_expression( new ConstantT( massTags[i], massFracs[i] ) ) );
    }
    tree.insert_tree( factory.register_expression( new SpeciesNT( massTags[nSpec-1], massTags ) ) );
    tree.insert_tree( factory.register_expression( new MixtureMWT( mmwTag, massTags, pokitt::MASS ) ) );
    tree.insert_tree( factory.register_expression( new PressureT( presTag, tempTag, rhoTag, mmwTag ) ) );

    Expr::TagList sensVarTags = {tempTag, rhoTag};
    Expr::TagList sensFxnTags = sensVarTags;
    for( int i=0; i<nSpec-1; ++i ){
      sensFxnTags.push_back( massTags[i] );
      sensVarTags.push_back( massTags[i] );
    }
    Expr::TagList computedTags = {presTag, mmwTag, massTags[nSpec-1]};
    for( const auto& t : computedTags ){
      sensFxnTags.push_back( t );
    }

    tree.compute_sensitivities( sensFxnTags, sensVarTags );

    tree.register_fields( fml );
    tree.bind_fields( fml );
    tree.lock_fields( fml );

    so::IntVec gridSize = so::IntVec( 20,  1,  1 );
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
    fml.field_ref<CellFieldT>( presTag ).add_device( CPU_INDEX );
    fml.field_ref<CellFieldT>( rhoTag  ).add_device( CPU_INDEX );
    fml.field_ref<CellFieldT>( mmwTag  ).add_device( CPU_INDEX );
    for( int i=0; i<nSpec; ++i ){
      fml.field_ref<CellFieldT>( massTags[i] ).add_device( CPU_INDEX );
    }
#   endif
    // need to do add_device on all sensitivity fields but this won't work until GPU sensitivity support is in ExprLib

    TestHelper allSens( true );

    TestHelper varToVar( false );

    for( const auto& f : sensVarTags ){
      for( const auto& v : sensVarTags ){
        const Expr::Tag ddTag = Expr::sens_tag( f, v  );

        if( f == v ){
          varToVar( so::field_equal( 1.0, fml.field_ref<CellFieldT>( ddTag ) ), ddTag.name() + " = 1.0" );
        }
        else{
          varToVar( so::field_equal( 0.0, fml.field_ref<CellFieldT>( ddTag ) ), ddTag.name() + " = 0.0" );
        }
      }
    }
    allSens( varToVar.ok(), "variable-to-variable identity matrix" );

    CellFieldT& rhoField = fml.field_ref<CellFieldT>( rhoTag );
    CellFieldT& pField   = fml.field_ref<CellFieldT>( presTag );
    CellFieldT& TField   = fml.field_ref<CellFieldT>( tempTag );
    CellFieldT& mmwField = fml.field_ref<CellFieldT>( mmwTag );
    CellFieldT& YnsField = fml.field_ref<CellFieldT>( massTags[nSpec-1] );
    CellFieldPtrT dpdrho = so::SpatialFieldStore::get<CellFieldT>( rhoField );
    CellFieldPtrT dpdT   = so::SpatialFieldStore::get<CellFieldT>( rhoField );

    TestHelper nthSpecies( false );
    const Expr::Tag& yNsTag = massTags[nSpec-1];
    for( int i=0; i<nSpec-1; ++i ){
      CellFieldPtrT dYnsdYi = so::SpatialFieldStore::get<CellFieldT>( rhoField );
      *dYnsdYi <<= -1.0;
      const Expr::Tag& yNsSensyi = Expr::sens_tag( yNsTag, massTags[i] );
      nthSpecies( so::field_equal( *dYnsdYi, fml.field_ref<CellFieldT>( yNsSensyi ), 1e-8 ), yNsSensyi.name() );
    }
    nthSpecies( so::field_equal( 0.0, fml.field_ref<CellFieldT>( Expr::sens_tag( yNsTag, rhoTag  ) ), 1e-8 ), yNsTag.name() + "_sens_rho" );
    nthSpecies( so::field_equal( 0.0, fml.field_ref<CellFieldT>( Expr::sens_tag( yNsTag, tempTag ) ), 1e-8 ), yNsTag.name() + "_sens_T" );
    allSens( nthSpecies.ok(), "n-th species" );


    TestHelper mixMolWeight( false );
    for( int i=0; i<nSpec-1; ++i ){
      CellFieldPtrT dmmwdYi = so::SpatialFieldStore::get<CellFieldT>( rhoField );
      *dmmwdYi <<= - ( 1.0 / mwVec[i] - 1.0 / mwVec[nSpec-1] ) * mmwField * mmwField;
      const Expr::Tag mwSensYi = Expr::sens_tag( mmwTag, massTags[i] );
      mixMolWeight( so::field_equal( *dmmwdYi, fml.field_ref<CellFieldT>( mwSensYi ), 1e-8 ), mwSensYi.name() );
    }
    mixMolWeight( so::field_equal( 0.0, fml.field_ref<CellFieldT>( Expr::sens_tag( mmwTag, rhoTag ) ), 1e-8 ), "mmw_sens_rho" );
    mixMolWeight( so::field_equal( 0.0, fml.field_ref<CellFieldT>( Expr::sens_tag( mmwTag, tempTag) ), 1e-8 ), "mmw_sens_T" );
    allSens( mixMolWeight.ok(), "mixture molecular weight" );


    TestHelper pressure( false );
    for( int i=0; i<nSpec-1; ++i ){
      const Expr::Tag pSensTag   = Expr::sens_tag( presTag, massTags[i] );
      const Expr::Tag mmwSensTag = Expr::sens_tag( mmwTag,  massTags[i] );

      CellFieldPtrT dpdYi = so::SpatialFieldStore::get<CellFieldT>( rhoField );
      *dpdYi <<= -pField / mmwField * fml.field_ref<CellFieldT>( mmwSensTag );

      pressure( so::field_equal( *dpdYi, fml.field_ref<CellFieldT>( pSensTag ), 1e-8 ), pSensTag.name() );
    }
    *dpdrho <<= pField / rhoField;
    *dpdT <<= pField / TField;
    pressure( so::field_equal( *dpdrho, fml.field_ref<CellFieldT>( Expr::sens_tag( presTag, rhoTag ) ), 1e-8 ), "p_sens_rho" );
    pressure( so::field_equal( *dpdT  , fml.field_ref<CellFieldT>( Expr::sens_tag( presTag, tempTag) ), 1e-8 ), "p_sens_T" );
    allSens( pressure.ok(), "pressure" );
    allVarSets( allSens.ok(), "set density" );
  }

  if( allVarSets.ok() ){
    std::cout << "\nPASS\n";
    return 0;
  }
  std::cout << "\nFAIL\n";
  return -1;
}
