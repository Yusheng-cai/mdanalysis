#pragma once
#include "tools/CommonTypes.h"
#include "tools/Assert.h"
#include "Eigen/Core"
#include "Eigen/Dense"
#include "tools/Constants.h"
#include "tools/CommonOperations.h"

#include <cmath>
#include <algorithm>
#include <numeric>

namespace LinAlg3x3
{
    using Real = CommonTypes::Real;
    using Matrix = CommonTypes::Matrix;
    using Real3= CommonTypes::Real3;
    static constexpr int N_DIM = 3;

    Real3 CrossProduct(const Real3& v1, const Real3& v2);
    Real DotProduct(const Real3& v1, const Real3& v2);
    Real norm(const Real3& v1);
    void normalize(Real3& v1);
    Matrix dyad(const Real3& v1, const Real3& v2);

    Real3 MatrixDotVector(const Matrix& mat1, const Real3& v1); 

    // This is the rotation matrix that rotates vector 1 onto vector 2
    Matrix GetRotationMatrix(const Real3& v1, const Real3& v2);
    Matrix RotationMatrix(const Real3& v1, const Real3& v2);

	// Rotation matrix around a vector 
	Matrix RodriguesRotationFormula(Real degree, const Real3& vec);

    // calculate 3uiuj - I 
    Matrix LocalQtensor(const Real3& v1);

	// find azimuthal angle 
	Real CalculateAzimuthalAngle(const Real3& v);

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
	Matrix vec_dyadic(const Real3&, const Real3&);
	Real3	vec_cross(const Real3&, const Real3&);

	// Miscellaneous matrix stuff
	Real   matrix_trace(const Matrix&);
	Matrix matrix_add(const Matrix&, const Matrix&);
	Matrix matrix_sub(const Matrix&, const Matrix&);
	Matrix matrix_Identity();
	Matrix matrix_dot(const Matrix&,const Matrix&);
	void   matrix_mult_inplace(Matrix&,Real);
	void   matrix_zero_inplace(Matrix&);
	void   matrix_accum_inplace(Matrix&, const Matrix&);

	// Function that returns p2 & corresponding v1 (eigenvectors that corresponds to second largest eigenvalue of Q matrix)
	std::pair<Real3,Matrix> EigenSolver(const Matrix& Q);
	std::pair<Real3, Matrix> OrderEigenSolver(const Matrix& Q);

	// Perform argsort of a Real3 array
	std::vector<size_t> argsort(const Real3& array);
}