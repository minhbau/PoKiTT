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
 *  \file   PressurePrimitivePrimitive.h
 *  \date   Apr 29, 2016
 *  \author mike
 */

#ifndef PRESSUREPRIMITIVEPRIMITIVE_H_
#define PRESSUREPRIMITIVEPRIMITIVE_H_

#include <numeric>
#include <iostream>
#include <vector>

#include <cmath>

#include <pokitt/CanteraObjects.h>

#include <spatialops/Nebo.h>
#include <spatialops/structured/MatVecFields.h>
#include <spatialops/structured/MatVecOps.h>

namespace pokitt{

  /**
   *  \class  PTYPrimitiveTransformer
   *  \author Mike Hansen
   *  \date 2016
   *
   *  \brief Calculates transformation matrices and performs sparse transformations
   *  between the primitives with pressure instead of rho and the standard primitives.
   *
   *  PressurePrimitive = (p  , T, {Y})
   *  Primitive         = (rho, T, {Y})
   *
   *  Currently only sparse right-multiplication is supported.
   *
   *  See the writeup for details.
   */
struct PTYPrimitiveTransformer
{
  PTYPrimitiveTransformer()
  : ns( CanteraObjects::number_species() ),\
    Msp( CanteraObjects::molecular_weights() ),
    Ru( CanteraObjects::gas_constant() ),
    invRu( 1.0 / CanteraObjects::gas_constant() ),
    pref( CanteraObjects::reference_pressure() )
  {
    for( size_t s=0; s<ns; ++s ){
      invMsp.push_back( 1.0 / Msp[s] );
    }
  } // constructor

  const size_t ns;

  const std::vector< double > Msp;
  std::vector< double > invMsp;

  const double Ru;
  const double invRu;

  const double pref;



  /**
   *  @brief Obtain the transformation matrix for \pder{Primitive}{PressurePrimitive}
   *  @param transformationMatrix a FieldMatrix for the transformation
   *  @param T a field containing the temperature
   *  @param rho a field containing the density
   *  @param YPtr a vector of SpatFldPtr containing the mass fraction fields
   *  @param Mmix a field containing the mixture molecular weight
   *
   *  sens_of_primitive means sensitivities of the primitive variables to the pressure-primitive variables
   *
   */
  template< typename FieldT >
  void sens_of_primitive_matrix( SpatialOps::FieldMatrix<FieldT>& transformationMatrix,
                                 const FieldT& T,
                                 const FieldT& rho,
                                 const std::vector< const FieldT* >& YPtr,
                                 const FieldT& Mmix ) const;

  /**
   *  @brief Obtain the transformation matrix for \pder{PressurePrimitive}{Primitive}
   *  @param transformationMatrix a FieldMatrix for the transformation
   *  @param T a field containing the temperature
   *  @param rho a field containing the density
   *  @param YPtr a vector of SpatFldPtr containing the mass fraction fields
   *  @param Mmix a field containing the mixture molecular weight
   *
   *  sens_to_primitive means sensitivities of the pressure-primitive variables to the primitive variables
   *
   */
  template< typename FieldT >
  void sens_to_primitive_matrix( SpatialOps::FieldMatrix<FieldT>& transformationMatrix,
                                 const FieldT& T,
                                 const FieldT& rho,
                                 const std::vector< const FieldT* >& YPtr,
                                 const FieldT& Mmix ) const;

  /**
   *  @brief Perform sparse matrix-matrix right-multiplication by \pder{Primitive}{PressurePrimitive}
   *  @param matrixToTransform a FieldMatrix to be transformed in place
   *  @param T a field containing the temperature
   *  @param rho a field containing the density
   *  @param YPtr a vector of SpatFldPtr containing the mass fraction fields
   *  @param Mmix a field containing the mixture molecular weight
   *
   *  right_multiply_to_conserved means sensitivities computed with respect to primitives are
   *  converted to those with respect to pressure-primitive, as
   *
   *  \pder{Function}{Primitive} * \pder{Primitive}{PressurePrimitive}
   *
   */
  template< typename FieldT >
  void right_multiply_to_pressure( SpatialOps::FieldMatrix<FieldT>& matrix,
                                   const FieldT& T,
                                   const FieldT& rho,
                                   const std::vector< const FieldT* >& YPtr,
                                   const FieldT& Mmix ) const;

  /**
   *  @brief Perform sparse matrix-matrix right-multiplication by \pder{PressurePrimitive}{Primitive}
   *  @param matrixToTransform a FieldMatrix to be transformed in place
   *  @param T a field containing the temperature
   *  @param rho a field containing the density
   *  @param YPtr a vector of SpatFldPtr containing the mass fraction fields
   *  @param Mmix a field containing the mixture molecular weight
   *
   *  right_multiply_to_primitive means sensitivities computed with respect to pressure-primitive are
   *  converted to those with respect to primitives, as
   *
   *  \pder{Function}{PressurePrimitive} * \pder{PressurePrimitive}{Primitive}
   *
   */
  template< typename FieldT >
  void right_multiply_to_primitive( SpatialOps::FieldMatrix<FieldT>& matrix,
                                    const FieldT& T,
                                    const FieldT& rho,
                                    const std::vector< const FieldT* >& YPtr,
                                    const FieldT& Mmix ) const;

};

template< typename FieldT >
void
PTYPrimitiveTransformer::
sens_of_primitive_matrix( SpatialOps::FieldMatrix<FieldT>& transformationMatrix,
                          const FieldT& T,
                          const FieldT& rho,
                          const std::vector< const FieldT* >& YPtr,
                          const FieldT& Mmix ) const
{
  for( size_t i=0; i<(ns+1)*(ns+1); ++i )
    transformationMatrix(i) <<= 0.0;
  for( size_t i=0; i<ns+1; ++i )
    transformationMatrix(i,i) <<= 1.0;

  transformationMatrix(0,0) <<= Mmix / ( T * Ru );
  transformationMatrix(0,1) <<= - rho / T;
  for( size_t s=0; s<ns-1; ++s ){
    transformationMatrix(0,2+s) <<= - rho * Mmix * ( invMsp[s] - invMsp[ns-1] );
  }
}

template< typename FieldT >
void
PTYPrimitiveTransformer::
sens_to_primitive_matrix( SpatialOps::FieldMatrix<FieldT>& transformationMatrix,
                          const FieldT& T,
                          const FieldT& rho,
                          const std::vector< const FieldT* >& YPtr,
                          const FieldT& Mmix ) const
{

  for( size_t i=0; i<(ns+1)*(ns+1); ++i )
    transformationMatrix(i) <<= 0.0;
  for( size_t i=0; i<ns+1; ++i )
    transformationMatrix(i,i) <<= 1.0;

  SpatialOps::SpatFldPtr<FieldT> rMixPtr = SpatialOps::SpatialFieldStore::get<FieldT>( T );
  FieldT& Rmix = *rMixPtr; // mixture specific gas constant

  Rmix <<= Ru / Mmix;

  transformationMatrix(0,0) <<= Rmix * T;
  transformationMatrix(0,1) <<= Rmix * rho;
  for( size_t s=0; s<ns-1; ++s ){
    transformationMatrix(0,2+s) <<= Ru * rho * T * ( invMsp[s] - invMsp[ns-1] );
  }
}

template< typename FieldT >
void
PTYPrimitiveTransformer::
right_multiply_to_pressure( SpatialOps::FieldMatrix<FieldT>& matrix,
                            const FieldT& T,
                            const FieldT& rho,
                            const std::vector< const FieldT* >& YPtr,
                            const FieldT& Mmix ) const
{
  SpatialOps::SpatFldPtr<FieldT> rMixPtr = SpatialOps::SpatialFieldStore::get<FieldT>( T );
  FieldT& Rmix = *rMixPtr; // mixture specific gas constant

  Rmix <<= Ru / Mmix;

  std::vector< SpatialOps::SpatFldPtr<FieldT> > temporaryRowPtr;
  for( size_t i=0; i<ns+1; ++i ){
    temporaryRowPtr.push_back( SpatialOps::SpatialFieldStore::get<FieldT>( T ) );
  }
  SpatialOps::FieldVector<FieldT> temporaryRow( temporaryRowPtr );

  // transform each row of matrixToTransform to resultantMatrix with a sparse matvec without the transformation matrix
  for( size_t row = 0; row<ns+1; ++row ){
    for( size_t col=0; col<ns+1; ++col ){
      temporaryRow(col) <<= matrix(row,col);
    }
    matrix(row,0) <<= temporaryRow(0) / ( T * Rmix );
    matrix(row,1) <<= - rho / T * temporaryRow(0) + temporaryRow(1);
    for( size_t s = 0; s<ns-1; ++s )
      matrix(row,2+s) <<= -rho * Mmix * ( invMsp[s] - invMsp[ns-1] ) * temporaryRow(0) + temporaryRow(2+s);
  }
}

template< typename FieldT >
void
PTYPrimitiveTransformer::
right_multiply_to_primitive( SpatialOps::FieldMatrix<FieldT>& matrix,
                             const FieldT& T,
                             const FieldT& rho,
                             const std::vector< const FieldT* >& YPtr,
                             const FieldT& Mmix ) const
{
  SpatialOps::SpatFldPtr<FieldT> rMixPtr = SpatialOps::SpatialFieldStore::get<FieldT>( T );
  FieldT& Rmix = *rMixPtr; // mixture specific gas constant

  Rmix <<= Ru / Mmix;

  std::vector< SpatialOps::SpatFldPtr<FieldT> > temporaryRowPtr;
  for( size_t i=0; i<ns+1; ++i ){
    temporaryRowPtr.push_back( SpatialOps::SpatialFieldStore::get<FieldT>( T ) );
  }
  SpatialOps::FieldVector<FieldT> temporaryRow( temporaryRowPtr );

  // transform each row of matrixToTransform to resultantMatrix with a sparse matvec without the transformation matrix
  for( size_t row = 0; row<ns+1; ++row ){
    for( size_t col=0; col<ns+1; ++col ){
      temporaryRow(col) <<= matrix(row,col);
    }
    matrix(row,0) <<= Rmix * T * temporaryRow(0);
    matrix(row,1) <<= Rmix * rho * temporaryRow(0) + temporaryRow(1);
    for( size_t s = 0; s<ns-1; ++s )
      matrix(row,2+s) <<= Ru * T * rho * ( invMsp[s] - invMsp[ns-1] ) * temporaryRow(0) + temporaryRow(2+s);
  }
}

}



#endif /* PRESSUREPRIMITIVEPRIMITIVE_H_ */
