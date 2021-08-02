#pragma once
#ifndef QTENSOR_H
#define QTENSOR_H

#include <array>
#include <algorithm>
#include <vector>
#include <numeric> 
#include <cmath>
#include "tools/CommonTypes.h"
#include "tools/Assert.h"

namespace Qtensor{
    // Scalars
  	using Real    = CommonTypes::Real;

	// Arrays
	using Real3    = CommonTypes::Real3;

	// Vectors
	using VectorReal     = CommonTypes::VectorReal;
	using VectorReal3    = CommonTypes::VectorReal3;

	// Matrices/tensors
	using Matrix = CommonTypes::Matrix; 

	const Real PI=3.14159265359;
	const int  N_DIM = 3;
	
	// Eigenvalue stuff

	// Function that calculates the eigenvector of a symmetric 3x3 matrix
	// using the Cayley-Hamilton theorem from https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3.C3.973_matrices
	//
	// Args:
	// 	Q(Qtensor::Matrix): The Q matrix, it has to be symmetric
	//
	// Returns:
	// 	eigenvalue(Qtensor::Real3): A array of 3 eigenvalues
	Real3 eigval_Qtensor(const Matrix&);

	// NOTE!!! This only works for 3x3 matrix (have only tested on symmetric traceless 3x3 matrices)
	//
	// Find the eigenvector of a corresponding eigenvalue of the matrix, eigenvalues are passed in as 
	// the original order
	Matrix eigvec_Qtensor(const Matrix&,const Real3&);
	std::pair<Real3,Matrix> eig_Qtensor(const Matrix&);

	// Miscellaneous vector stuff
	Real3  vec_add(const Real3&, const Real3&);
	Real3  vec_add(const Real3&, const Real);
	Real3  vec_add(const Real, const Real3&);
	Real3  vec_sub(const Real3&, const Real3&);
	Real3  vec_mult(const Real3&, const Real3&);
	Real3  vec_mult(const Real3&, const Real);
	Real3  vec_mult(const Real, const Real3&);
	Real   vec_dot(const Real3&, const Real3&);
	Real   vec_norm(const Real3);
	Matrix vec_dyadic(const Qtensor::Real3&, const Qtensor::Real3&);
	Real3	vec_cross(const Qtensor::Real3&, const Qtensor::Real3&);

	// Miscellaneous matrix stuff
	Real   matrix_trace(const Matrix&);
	Matrix matrix_add(const Matrix&, const Matrix&);
	Matrix matrix_sub(const Matrix&, const Matrix&);
	Matrix matrix_Identity();
	Matrix matrix_dot(const Matrix&,const Matrix&);
	void   matrix_mult_inplace(Matrix&,Real);
	void   matrix_zero_inplace(Matrix&);
	void   matrix_accum_inplace(Matrix&, const Matrix&);

	// Order parameter stuff

	// Function that returns p2 & corresponding v1 (eigenvectors that corresponds to second largest eigenvalue of Q matrix)
	std::pair<Real,Real3> OP_Qtensor(const Matrix& Q);

	// Perform argsort of a Real3 array
	std::vector<size_t> argsort(const Real3& array);

	// Returns ordered eigenvector & eigenvalues from largest to smallest according to the absolute value of eigenvalue, 0 is largest, 2 is smallest
	// in absolute language
	std::pair<VectorReal3, Real3> orderedeig_Qtensor(const Matrix& Q);
	Real3 normalize_director(const Real3&);
}

#endif