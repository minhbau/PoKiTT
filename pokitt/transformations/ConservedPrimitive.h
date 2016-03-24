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
 *  \file   ConservedPrimitive.h
 *  \date   Apr 29, 2016
 *  \author mike
 */

#ifndef CONSERVEDPRIMITIVE_H_
#define CONSERVEDPRIMITIVE_H_

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
   *  \class  ConservedPrimitiveTransformer
   *  \author Mike Hansen
   *  \date 2016
   *
   *  \brief Calculates transformation matrices and performs sparse transformations
   *  between the conserved variables and primitives.
   *
   *  Conserved = (rho, rho e_0, {rho Y})
   *  Primitive = (rho, T      , {Y})
   *
   *  Currently only sparse right-multiplication is supported.
   *
   *  See the writeup for details.
   */
struct ConservedPrimitiveTransformer
{
  ConservedPrimitiveTransformer()
  : ns( CanteraObjects::number_species() ),\
    Msp( CanteraObjects::molecular_weights() ),
    Ru( CanteraObjects::gas_constant() ),
    invRu( 1.0 / CanteraObjects::gas_constant() ),
    pref( CanteraObjects::reference_pressure() )
  {
    for( size_t s=0; s<ns; ++s ){
      invMsp.push_back( 1.0 / Msp[s] );
    }
  }

  const size_t ns;

  const std::vector< double > Msp;
  std::vector< double > invMsp;

  const double Ru;
  const double invRu;

  const double pref;

  /**
   *  @brief Obtain the transformation matrix for \pder{Primitive}{Conserved}
   *  @param transformationMatrix a FieldMatrix for the transformation
   *  @param T a field containing the temperature
   *  @param rho a field containing the density
   *  @param YPtr a vector of const FieldT* containing the mass fraction fields
   *  @param cv mixture specific heat capacity
   *  @param speciesEnergyPtr a vector of const FieldT* containing the species energy fields
   *
   *  sens_of_primitive means sensitivities of the primitive variables to the conserved variables
   *
   */
  template< typename FieldT >
  void sens_of_primitive_matrix( SpatialOps::FieldMatrix<FieldT>& transformationMatrix,
                                 const FieldT& T,
                                 const FieldT& rho,
                                 const std::vector< const FieldT* >& YPtr,
                                 const FieldT& cv,
                                 const std::vector< const FieldT*>& speciesEnergyPtr ) const;

  /**
   *  @brief Obtain the transformation matrix for \pder{Conserved}{Primitive}
   *  @param transformationMatrix a FieldMatrix for the transformation
   *  @param T a field containing the temperature
   *  @param rho a field containing the density
   *  @param YPtr a vector of const FieldT* containing the mass fraction fields
   *  @param cv mixture specific heat capacity
   *  @param etotal mixture specific total energy
   *  @param speciesEnergyPtr a vector of const FieldT* containing the species energy fields
   *
   *  sens_to_primitive means sensitivities of the conserved variables to the primitive variables
   *
   */
  template< typename FieldT >
  void sens_to_primitive_matrix( SpatialOps::FieldMatrix<FieldT>& transformationMatrix,
                                 const FieldT& T,
                                 const FieldT& rho,
                                 const std::vector< const FieldT* >& YPtr,
                                 const FieldT& cv,
                                 const FieldT& etotal,
                                 const std::vector< const FieldT*>& speciesEnergyPtr ) const;

  /**
   *  @brief Perform sparse matrix-matrix right-multiplication by \pder{Primitive}{Conserved}
   *  @param matrix a FieldMatrix to be transformed in place
   *  @param T a field containing the temperature
   *  @param rho a field containing the density
   *  @param YPtr a vector of const FieldT* containing the mass fraction fields
   *  @param cv mixture specific heat capacity
   *  @param etotal mixture specific total energy
   *  @param speciesEnergyPtr a vector of const FieldT* containing the species energy fields
   *
   *  right_multiply_to_conserved means sensitivities computed with respect to primitives are
   *  converted to those with respect to conserved, as
   *
   *  \pder{Function}{Primitive} * \pder{Primitive}{Conserved}
   *
   */
  template< typename FieldT >
  void right_multiply_to_conserved( SpatialOps::FieldMatrix<FieldT>& matrix,
                                    const FieldT& T,
                                    const FieldT& rho,
                                    const std::vector< const FieldT* >& YPtr,
                                    const FieldT& cv,
                                    const std::vector< const FieldT*>& speciesEnergyPtr ) const;

  /**
   *  @brief Perform sparse matrix-matrix right-multiplication by \pder{Conserved}{Primitive}
   *  @param matrix a FieldMatrix to be transformed in place
   *  @param T a field containing the temperature
   *  @param rho a field containing the density
   *  @param YPtr a vector of const FieldT* containing the mass fraction fields
   *  @param cv mixture specific heat capacity
   *  @param speciesEnergyPtr a vector of const FieldT* containing the species energy fields
   *
   *  right_multiply_to_primitive means sensitivities computed with respect to conserved are
   *  converted to those with respect to primitives, as
   *
   *  \pder{Function}{Conserved} * \pder{Conserved}{Primitive}
   *
   */
  template< typename FieldT >
  void right_multiply_to_primitive( SpatialOps::FieldMatrix<FieldT>& matrix,
                                    const FieldT& T,
                                    const FieldT& rho,
                                    const std::vector< const FieldT* >& YPtr,
                                    const FieldT& cv,
                                    const FieldT& etotal,
                                    const std::vector< const FieldT*>& speciesEnergyPtr ) const;

};

template< typename FieldT >
void
ConservedPrimitiveTransformer::
sens_of_primitive_matrix( SpatialOps::FieldMatrix<FieldT>& transformationMatrix,
                          const FieldT& T,
                          const FieldT& rho,
                          const std::vector< const FieldT* >& YPtr,
                          const FieldT& cv,
                          const std::vector< const FieldT*>& speciesEnergyPtr ) const
{
  std::vector< SpatialOps::SpatFldPtr<FieldT> > specOffEgyPtr; // offset species specific energies, theta_i = e_i - e_{n_s}, a vector of ns quantities
  for( size_t i=0; i<ns; ++i ){
    specOffEgyPtr.push_back( SpatialOps::SpatialFieldStore::get<FieldT>( T ) );
  }

  for( size_t n=0; n<ns; ++n ){
    *specOffEgyPtr[n] <<= *speciesEnergyPtr[n] - *speciesEnergyPtr[ns-1];
  }

  for( size_t i=0; i<(ns+1)*(ns+1); ++i )
    transformationMatrix(i) <<= 0.0;

  transformationMatrix(0,0) <<= 1.0;
  transformationMatrix(1,0) <<= -*speciesEnergyPtr[ns-1] / rho / cv;
  transformationMatrix(1,1) <<= 1 / rho / cv;
  for( size_t s=0; s<ns-1; ++s ){
    transformationMatrix(1,2+s)   <<= -*specOffEgyPtr[s] / rho / cv;
    transformationMatrix(2+s,0)   <<= -*YPtr[s] / rho;
    transformationMatrix(2+s,2+s) <<= 1 / rho;
  }
}

template< typename FieldT >
void
ConservedPrimitiveTransformer::
sens_to_primitive_matrix( SpatialOps::FieldMatrix<FieldT>& transformationMatrix,
                       const FieldT& T,
                       const FieldT& rho,
                       const std::vector< const FieldT* >& YPtr,
                       const FieldT& cv,
                       const FieldT& etotal,
                       const std::vector< const FieldT*>& speciesEnergyPtr ) const
{
  std::vector< SpatialOps::SpatFldPtr<FieldT> > specOffEgyPtr; // offset species specific energies, theta_i = e_i - e_{n_s}, a vector of ns quantities
  for( size_t i=0; i<ns; ++i ){
    specOffEgyPtr.push_back( SpatialOps::SpatialFieldStore::get<FieldT>( T ) );
  }

  for( size_t n=0; n<ns; ++n ){
    *specOffEgyPtr[n] <<= *speciesEnergyPtr[n] - *speciesEnergyPtr[ns-1];
  }

  for( size_t i=0; i<(ns+1)*(ns+1); ++i )
    transformationMatrix(i) <<= 0.0;

  transformationMatrix(0,0) <<= 1.0;
  transformationMatrix(1,0) <<= etotal;
  transformationMatrix(1,1) <<= rho * cv;
  for( size_t s=0; s<ns-1; ++s ){
    transformationMatrix(1,2+s)   <<= rho * *specOffEgyPtr[s];
    transformationMatrix(2+s,0)   <<= *YPtr[s];
    transformationMatrix(2+s,2+s) <<= rho;
  }
}

template< typename FieldT >
void
ConservedPrimitiveTransformer::
right_multiply_to_conserved( SpatialOps::FieldMatrix<FieldT>& matrix,
                             const FieldT& T,
                             const FieldT& rho,
                             const std::vector< const FieldT* >& YPtr,
                             const FieldT& cv,
                             const std::vector< const FieldT*>& speciesEnergyPtr ) const
{
  std::vector< SpatialOps::SpatFldPtr<FieldT> > specOffEgyPtr; // offset species specific energies, theta_i = e_i - e_{n_s}, a vector of ns quantities
  for( size_t i=0; i<ns; ++i ){
    specOffEgyPtr.push_back( SpatialOps::SpatialFieldStore::get<FieldT>( T ) );
  }

  for( size_t n=0; n<ns; ++n ){
    *specOffEgyPtr[n] <<= *speciesEnergyPtr[n] - *speciesEnergyPtr[ns-1];
  }

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

    matrix(row,0) <<= rho * temporaryRow(0) - *speciesEnergyPtr[ns-1] / cv * temporaryRow(1);
    for( size_t s = 0; s<ns-1; ++s ){
      matrix(row,0) <<= matrix(row,0) - *YPtr[s] * temporaryRow(2+s);
    }
    matrix(row,0) <<= matrix(row,0) / rho;

    matrix(row,1) <<= temporaryRow(1) / ( rho * cv );

    for( size_t s = 0; s<ns-1; ++s ){
      matrix(row,2+s) <<= ( temporaryRow(2+s) - temporaryRow(1) * *specOffEgyPtr[s] / cv ) / rho;
    }
  }
}

template< typename FieldT >
void
ConservedPrimitiveTransformer::
right_multiply_to_primitive( SpatialOps::FieldMatrix<FieldT>& matrix,
                             const FieldT& T,
                             const FieldT& rho,
                             const std::vector< const FieldT* >& YPtr,
                             const FieldT& cv,
                             const FieldT& etotal,
                             const std::vector< const FieldT*>& speciesEnergyPtr ) const
{
  std::vector< SpatialOps::SpatFldPtr<FieldT> > specOffEgyPtr; // offset species specific energies, theta_i = e_i - e_{n_s}, a vector of ns quantities
  for( size_t i=0; i<ns; ++i ){
    specOffEgyPtr.push_back( SpatialOps::SpatialFieldStore::get<FieldT>( T ) );
  }

  for( size_t n=0; n<ns; ++n ){
    *specOffEgyPtr[n] <<= *speciesEnergyPtr[n] - *speciesEnergyPtr[ns-1];
  }

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

    matrix(row,0) <<= temporaryRow(0) + etotal * temporaryRow(1);
    for( size_t s = 0; s<ns-1; ++s )
      matrix(row,0) <<= matrix(row,0) + *YPtr[s] * temporaryRow(2+s);

    matrix(row,1) <<= temporaryRow(1) * rho * cv;

    for( size_t s = 0; s<ns-1; ++s )
      matrix(row,2+s) <<= *specOffEgyPtr[s] * rho * temporaryRow(1) + rho * temporaryRow(2+s);
  }
}

}


#endif /* CONSERVEDPRIMITIVE_H_ */
