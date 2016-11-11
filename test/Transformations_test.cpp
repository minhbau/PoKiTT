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
 * THE SOFTQARE IS PROVIDED "AS IS", QITHOUT QARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE QARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, QHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERQISE, ARISING
 * FROM, OUT OF OR IN CONNECTION QITH THE SOFTQARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTQARE.
 */

/**
 *  \file   Transformations_test.cpp
 *  \date   May 1, 2016
 *  \author mike
 */

#include <numeric>
#include <iostream>
#include <vector>
#include <cmath>

#include "TestHelper.h"

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include <spatialops/Nebo.h>
#include <spatialops/structured/Grid.h>
#include <spatialops/structured/MatVecFields.h>
#include <spatialops/structured/MatVecOps.h>
#include <spatialops/structured/FieldComparisons.h>

#include <expression/ExprLib.h>
#include <expression/Expression.h>

#include <pokitt/CanteraObjects.h>

#include <pokitt/transformations/ConservedPrimitive.h>
#include <pokitt/transformations/PTYPrimitive.h>
#include <pokitt/thermo/InternalEnergy.h>
#include <pokitt/thermo/HeatCapacity_Cv.h>

namespace pokitt{};
using namespace pokitt;
namespace po = boost::program_options;

#ifdef ENABLE_CUDA
# define LOCATION GPU_INDEX
#else
# define LOCATION CPU_INDEX
#endif

typedef SpatialOps::SVolField FieldT;

// --------------------------------
template< typename FieldT >
bool test_identity_matrix( const SpatialOps::FieldMatrix<FieldT>& fieldMatrix )
{
  const size_t ncols = std::sqrt( fieldMatrix.elements() );
  bool pass = true;
  for( size_t i=0; i<ncols; ++i ){
    for( size_t j=0; j<ncols; ++j ){
      if( i==j ) pass = pass && field_equal( 1.0, fieldMatrix(i,j), 1.0, 1e-8 );
      else       pass = pass && field_equal( 0.0, fieldMatrix(i,j), 1.0, 1e-8 );
    }
  }
  return pass;
}
// --------------------------------
template< typename FieldT >
void equate_matrices(       SpatialOps::FieldMatrix<FieldT>& dst,
                      const SpatialOps::FieldMatrix<FieldT>& src )
{
  bool pass = true;
  const size_t nelem = dst.elements();
  assert( dst.elements() == src.elements() );
  for( size_t i=0; i<nelem; ++i ){
      dst(i) <<= src(i);
  }
}
// --------------------------------


int main( int iarg, char* carg[] )
{
  // ---------------- parse options ----------------

  std::string inputFileName = "h2-burke.xml"; // mechanism input in xml format

  double temperature = 1000;
  double pressure = 101325;

  try {
    po::options_description desc("Supported Options");
    desc.add_options()
      ( "help", "print help message" )
      ( "xml-input-file", po::value<std::string>(&inputFileName)->default_value(inputFileName), "Cantera xml input file name" );

    po::variables_map args;
    po::store( po::parse_command_line(iarg,carg,desc), args );
    po::notify(args);

    if (args.count("help")) {
      std::cout << desc << "\n";
      return -1;
    }
  }
  catch( std::exception& err ){
    std::cout << "Error parsing input arguments\n" << err.what() << std::endl;
    return -2;
  }

  // ---------------- set up Cantera ----------------
  const CanteraObjects::Setup setup( "Mix", inputFileName );
  CanteraObjects::setup_cantera( setup );
  const size_t ns = CanteraObjects::number_species();

  // ---------------- build the grid ----------------
  const SpatialOps::IntVec    npts   ( 1, 1, 1 );
  const SpatialOps::DoubleVec lengths( 1, 1, 1 );
  const SpatialOps::Grid              grid( npts, lengths );
  const SpatialOps::DoubleVec         spacings = grid.spacing();
  const SpatialOps::GhostData         nghost(1);
  const SpatialOps::BoundaryCellInfo  bcInfo = SpatialOps::BoundaryCellInfo::build<FieldT>( npts[0]>1, npts[1]>1, npts[2]>1 );
  const SpatialOps::MemoryWindow      window( get_window_with_ghost( npts, nghost, bcInfo ) );

  // ---------------- initialize an equimolar mixture ----------------

  std::vector< double > Msp = CanteraObjects::molecular_weights();

  std::vector< double > invMsp( ns );
  for( size_t s=0; s<ns; ++s ){
    invMsp[s] = 1.0 / Msp[s];
  }
  const double Ru    = CanteraObjects::gas_constant();
  const double invRu = 1.0/Ru;

  FieldT T( window, bcInfo, nghost, NULL, SpatialOps::InternalStorage, LOCATION ); // temperature

  SpatialOps::SpatFldPtr<FieldT> pPtr    = SpatialOps::SpatialFieldStore::get<FieldT>( T );
  SpatialOps::SpatFldPtr<FieldT> rhoPtr  = SpatialOps::SpatialFieldStore::get<FieldT>( T );
  SpatialOps::SpatFldPtr<FieldT> MmixPtr = SpatialOps::SpatialFieldStore::get<FieldT>( T );
  SpatialOps::SpatFldPtr<FieldT> RmixPtr = SpatialOps::SpatialFieldStore::get<FieldT>( T );

  FieldT& p = *pPtr;       // pressure
  FieldT& rho = *rhoPtr;   // density
  FieldT& Mmix = *MmixPtr; // mixture molecular weight
  FieldT& Rmix = *RmixPtr; // specific mixture gas constant

  std::vector< SpatialOps::SpatFldPtr<FieldT> > YPtr;  // mass fractions
  std::vector< SpatialOps::SpatFldPtr<FieldT> > XPtr;  // mole fractions
  for( size_t i=0; i<ns; ++i ){
    YPtr.push_back( SpatialOps::SpatialFieldStore::get<FieldT>( T ) ); *YPtr[i] <<= 0;
    XPtr.push_back( SpatialOps::SpatialFieldStore::get<FieldT>( T ) ); *XPtr[i] <<= 0;
  }

  T <<= temperature;
  p <<= pressure;

  Mmix <<= 0;
  for( size_t i=0; i<ns; ++i ){
    *XPtr[i] <<= 1.0 / ( (double) ns );
    Mmix <<= Mmix + *XPtr[i] * Msp[i];
  }

  Rmix <<= Ru / Mmix;
  for( size_t i=0; i<ns; ++i ){
    *YPtr[i] <<= *XPtr[i] * Msp[i] / Mmix;
  }

  rho <<= p / T / Rmix;

  // ---------------- run the tests ----------------

  TestHelper status;

  pokitt::ConservedPrimitiveTransformer stateTransformUV;
  pokitt::PTYPrimitiveTransformer  stateTransformQV;

  SpatialOps::FieldMatrix<FieldT> dUdV ( ns+1, window, bcInfo, nghost, LOCATION );
  SpatialOps::FieldMatrix<FieldT> dVdU ( ns+1, window, bcInfo, nghost, LOCATION );
  SpatialOps::FieldMatrix<FieldT> dQdV ( ns+1, window, bcInfo, nghost, LOCATION );
  SpatialOps::FieldMatrix<FieldT> dVdQ ( ns+1, window, bcInfo, nghost, LOCATION );

  SpatialOps::FieldMatrix<FieldT> UV_dense ( ns+1, window, bcInfo, nghost, LOCATION );
  SpatialOps::FieldMatrix<FieldT> VU_dense ( ns+1, window, bcInfo, nghost, LOCATION );
  SpatialOps::FieldMatrix<FieldT> UV_sparse( ns+1, window, bcInfo, nghost, LOCATION );
  SpatialOps::FieldMatrix<FieldT> VU_sparse( ns+1, window, bcInfo, nghost, LOCATION );
  SpatialOps::FieldMatrix<FieldT> VQ_dense ( ns+1, window, bcInfo, nghost, LOCATION );
  SpatialOps::FieldMatrix<FieldT> QV_dense ( ns+1, window, bcInfo, nghost, LOCATION );
  SpatialOps::FieldMatrix<FieldT> QV_sparse( ns+1, window, bcInfo, nghost, LOCATION );
  SpatialOps::FieldMatrix<FieldT> VQ_sparse( ns+1, window, bcInfo, nghost, LOCATION );

  std::vector< const FieldT* > YPtrForMat;
  for( size_t i=0; i<YPtr.size(); ++i ){
    YPtrForMat.push_back( &( *YPtr[i] ) );
  }

  // compute additional fields necessary for transformation
  const Expr::Tag tempTag( "T", Expr::STATE_NONE );
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

  Expr::ExpressionTree tree( roots, factory, 0 );
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

#   ifdef ENABLE_CUDA
  for( size_t rowIdx=0; rowIdx<ns+1; ++rowIdx ){
    for( size_t colIdx=0; colIdx<ns+1; ++colIdx ){
      dVdU(rowIdx,colIdx).set_device_as_active( GPU_INDEX );
      dUdV(rowIdx,colIdx).set_device_as_active( GPU_INDEX );
      dVdQ(rowIdx,colIdx).set_device_as_active( GPU_INDEX );
      dQdV(rowIdx,colIdx).set_device_as_active( GPU_INDEX );
    }
  }
#   endif

  // test matrices
  stateTransformUV.sens_of_primitive_matrix( dVdU, T, rho, YPtrForMat, cv, speciesEnergyPtr );
  stateTransformUV.sens_to_primitive_matrix( dUdV, T, rho, YPtrForMat, cv, etotal, speciesEnergyPtr );
  stateTransformQV.sens_of_primitive_matrix( dVdQ, T, rho, YPtrForMat, Mmix );
  stateTransformQV.sens_to_primitive_matrix( dQdV, T, rho, YPtrForMat, Mmix );

  UV_dense = dUdV * dVdU;
  VU_dense = dVdU * dUdV;
  equate_matrices( UV_sparse, dUdV );
  equate_matrices( VU_sparse, dVdU );
  stateTransformUV.right_multiply_to_conserved( UV_sparse, T, rho, YPtrForMat, cv, speciesEnergyPtr );
  stateTransformUV.right_multiply_to_primitive( VU_sparse, T, rho, YPtrForMat, cv, etotal, speciesEnergyPtr );
  status( test_identity_matrix( UV_dense  ), "dense  dUdV * dVdU = I" );
  status( test_identity_matrix( VU_dense  ), "dense  dVdU * dUdV = I" );
  status( test_identity_matrix( UV_sparse ), "sparse dUdV * dVdU = I" );
  status( test_identity_matrix( VU_sparse ), "sparse dVdU * dUdV = I" );
  std::cout << std::endl;

  QV_dense = dQdV * dVdQ;
  VQ_dense = dVdQ * dQdV;
  equate_matrices( QV_sparse, dQdV );
  equate_matrices( VQ_sparse, dVdQ );
  stateTransformQV.right_multiply_to_pressure ( QV_sparse, T, rho, YPtrForMat, Mmix );
  stateTransformQV.right_multiply_to_primitive( VQ_sparse, T, rho, YPtrForMat, Mmix );
  status( test_identity_matrix( QV_dense  ), "dense  dQdV * dVdQ = I" );
  status( test_identity_matrix( VQ_dense  ), "dense  dVdQ * dQdV = I" );
  status( test_identity_matrix( QV_sparse ), "sparse dQdV * dVdQ = I" );
  status( test_identity_matrix( VQ_sparse ), "sparse dVdQ * dQdV = I" );
  std::cout << std::endl;


}


