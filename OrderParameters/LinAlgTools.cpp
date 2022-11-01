#include "LinAlgTools.h"

LinAlg3x3::Real3 LinAlg3x3::CrossProduct(const Real3& v1, const Real3& v2)
{
    Real3 ret;
    ret.fill(0);

    ret[0] = v1[1]*v2[2] - v1[2]*v2[1];
    ret[1] = v1[2]*v2[0] - v1[0]*v2[2];
    ret[2] = v1[0]*v2[1] - v1[1]*v2[0];

    return ret;
}

LinAlg3x3::Real LinAlg3x3::CalculateAzimuthalAngle(const Real3& v1){
  // azimuthal angle phi (0, 2 pi)
  // first quadrant
  Real phi = std::atan2(v1[1],v1[0]);
  if (phi < 0){
      phi += 2 * Constants::PI;
  }

  return phi;
}

LinAlg3x3::Matrix LinAlg3x3::RodriguesRotationFormula(Real degree, const Real3& vec){
  Real radian = degree * Constants::PI / 180.0;

  Matrix W = {{{0, -vec[2], vec[1]},{vec[2], 0, -vec[1]},{-vec[1], vec[0], 0}}};

  Matrix W2 = LinAlg3x3::matrix_dot(W,W);
  Real sinphi = std::sin(radian);
  Real sinsq  = std::pow(std::sin(radian/2),2.0);

  Matrix ret = LinAlg3x3::matrix_Identity() +  W * sinphi + W2 * 2.0 * sinsq;

  return ret;
}

LinAlg3x3::Real LinAlg3x3::DotProduct(const Real3& v1, const Real3& v2)
{
    Real ret = 0.0;

    for (int i=0;i<3;i++)
    {
        ret += v1[i]*v2[i];
    }

    return ret;
}

LinAlg3x3::Matrix LinAlg3x3::LocalQtensor(const Real3& v1)
{
    Matrix Qtemp = dyad(v1, v1);

    for (int i=0;i<3;i++)
    {
        for (int j=0;j<3;j++)
        {
            Qtemp[i][j] = 3 * Qtemp[i][j];

            if (i == j)
            {
                Qtemp[i][j] -= 1;
            }
        }
    }

    return Qtemp;
}

LinAlg3x3::Real LinAlg3x3::norm(const Real3& v1)
{
    Real ret = 0.0;

    for (int i=0;i<3;i++)
    {
        ret += v1[i]*v1[i];
    }

    ret = std::sqrt(ret);

    return ret;
}

LinAlg3x3::Matrix LinAlg3x3::dyad(const Real3& v1, const Real3& v2)
{
    Matrix ret;

    for (int i=0;i<3;i++)
    {
        for (int j=0;j<3;j++)
        {
            ret[i][j] = v1[i]*v2[j];
        }
    }

    return ret;
}

void LinAlg3x3::normalize(Real3& v1)
{
    Real norm_ = norm(v1);

    for (int i=0;i<3;i++)
    {
        v1[i] = v1[i]/norm_;
    }
}

LinAlg3x3::Real3 LinAlg3x3::MatrixDotVector(const Matrix& A, const Real3& v1)
{
    Real3 ans;
    ans.fill(0);

    for (int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            ans[i] += A[i][j]*v1[j];
        }
    }

    return ans;
}

LinAlg3x3::Matrix LinAlg3x3::RotationMatrix(const Real3& v1, const Real3& v2)
{
  Matrix rotMat;
	Real rdotz      = DotProduct(v1, v2);
	Real3 cross 	= CrossProduct(v1, v2);
	Real x  = cross[0];
	Real y  = cross[1];
	Real z  = cross[2];

	Real auxk;
	if (std::abs(rdotz + 1) < 1e-7) auxk = 0;
	else auxk = 1/(1+rdotz);

	rotMat[0][0] = (x * x * auxk) + rdotz;
	rotMat[0][1] = (y * x * auxk) - z;
	rotMat[0][2] = (z * x * auxk) + y;
	rotMat[1][0] = (x * y * auxk) + z;
	rotMat[1][1] = (y * y * auxk) + rdotz;
	rotMat[1][2] = (z * y * auxk) - x;
	rotMat[2][0] = (x * z * auxk) - y;
	rotMat[2][1] = (y * z * auxk) + x;
	rotMat[2][2] = (z * z * auxk) + rdotz;

    return rotMat;
}

LinAlg3x3::Matrix LinAlg3x3::GetRotationMatrix(const Real3& v1, const Real3& v2)
{
    Real3 crossProduct = LinAlg3x3::CrossProduct(v1, v2);
    Real cosine = LinAlg3x3::DotProduct(v1, v2);

    Matrix ret;
    Real denom = 1.0 + cosine;
    Real factor;

    if (std::abs(denom) < 1e-7) { 
        factor = 0.0;
        ret.fill({});
        ret[0][0] = -1;
        ret[1][1] = -1;
        ret[2][2] = -1;
    }
    else{factor = 1.0/(1.0 + cosine);}


    ret[0][0] = 1 + factor*(-crossProduct[2]*crossProduct[2] - crossProduct[1]*crossProduct[1]);
    ret[0][1] = -crossProduct[2] + factor*crossProduct[0]*crossProduct[1];
    ret[0][2] = crossProduct[1] + factor*crossProduct[0]*crossProduct[2];
    ret[1][0] = crossProduct[2] + factor*crossProduct[0]*crossProduct[1];
    ret[1][1] = 1 + factor*(-crossProduct[2]*crossProduct[2] - crossProduct[0]*crossProduct[0]);
    ret[1][2] = -crossProduct[0] + factor*crossProduct[1]*crossProduct[2];
    ret[2][0] = -crossProduct[1] + factor*crossProduct[2]*crossProduct[0];
    ret[2][1] = crossProduct[0] + factor*crossProduct[2]*crossProduct[1];
    ret[2][2] = 1 + factor*(-crossProduct[1]*crossProduct[1] - crossProduct[0]*crossProduct[0]);

    return ret;
}

LinAlg3x3::Real3 LinAlg3x3::vec_add(const LinAlg3x3::Real3& v1, const LinAlg3x3::Real3& v2)
{
  // Function that sum a vector with another vector
  Real3 eval;
  for(int i = 0; i < N_DIM; i++){
    eval[i] = v1[i]+v2[i];
  }
  return eval;
}

LinAlg3x3::Real3 LinAlg3x3::vec_add(const Real3& v1, const Real s){
	Real3 eval;
	for(int i=0;i<3;i++){
		eval[i] = v1[i] + s;
	}

	return eval;
}

LinAlg3x3::Real3 LinAlg3x3::vec_add(const Real s, const Real3& v1){
	return LinAlg3x3::vec_add(v1,s);
}

LinAlg3x3::Real3 LinAlg3x3::vec_sub(const Real3& v1, const Real3& v2)
{ 
  Real3 eval;
  for(int i = 0; i < N_DIM; i++){
    eval[i] = v1[i]-v2[i];
  }
  return eval;
}


LinAlg3x3::Real3 LinAlg3x3::vec_mult(const Real3& v1, const Real3& v2)
{
  Real3 eval;
  for(int i = 0; i < N_DIM; i++){
    eval[i] = v1[i]*v2[i];
  }
  return eval;
}

LinAlg3x3::Real3 LinAlg3x3::vec_mult(const Real3& v1, const Real scalar_value)
{
  Real3 eval;
  for(int i = 0; i < N_DIM; i++){
    eval[i] = v1[i]*scalar_value;
  }
  return eval;
}

LinAlg3x3::Real3 LinAlg3x3::vec_mult(const Real scalar_value, const Real3& v1){
	return LinAlg3x3::vec_mult(v1,scalar_value);
}

LinAlg3x3::Real LinAlg3x3::vec_dot(const Real3& v1, const Real3& v2){
	Real dot_product = 0.0;
	for(int i=0;i<3;i++){
		dot_product += v1[i]*v2[i];
	}

	return dot_product;
}

LinAlg3x3::Real LinAlg3x3::vec_norm(Real3 v1)
{
  Real eval=LinAlg3x3::vec_dot(v1,v1);

  return std::sqrt(eval);
}

LinAlg3x3::Matrix LinAlg3x3::vec_dyadic(const Real3& v1, const Real3& v2)
{
  Matrix eval;
  for(int i = 0; i < N_DIM; i++){
    for(int j = 0; j < N_DIM; j++){
    eval[i][j] = v1[i]*v2[j];
    }
  }
  return eval;
}

LinAlg3x3::Real3 LinAlg3x3::vec_cross(const Real3& v1, const Real3& v2){
	Real3 output;
	output.fill(0);

	output[0] 	= v1[1]*v2[2] - v1[2]*v2[1];	
	output[1]   = v1[2]*v2[0] - v1[0]*v2[2];
	output[2]   = v1[0]*v2[1] - v1[1]*v2[0];  

	return output;
}

LinAlg3x3::Real LinAlg3x3::matrix_trace(const Matrix& Q){
	// Function that returns the trace of a 3x3 matrix (Q tensor)
	// 
	// Args:
	// Q (LinAlg3x3::Matrix): 3x3 matrix of the Q tensor
	//
	// returns:
	// trace (LinAlg3x3::Real): the trace of the Q tensor
	
	
	Real trace = 0.0;
	for(int i=0;i<3;i++){
		trace = trace + Q[i][i]; 
	}
	return trace;
}

LinAlg3x3::Matrix LinAlg3x3::matrix_add(const Matrix& m1, const Matrix& m2)
{
  Matrix eval;
  for(int i = 0; i < N_DIM; i++){
    for(int j = 0; j < N_DIM; j++){
      eval[i][j] = m1[i][j]+m2[i][j];
    }
  }
  return eval;
}

void LinAlg3x3::matrix_accum_inplace(Matrix& m1, const Matrix& m2)
{
  for(int i = 0; i < N_DIM; i++){
    for(int j = 0; j < N_DIM; j++){
      m1[i][j] = m1[i][j]+m2[i][j];
    }
  }
  return;
}
LinAlg3x3::Matrix LinAlg3x3::matrix_sub(const Matrix& m1, const Matrix& m2)
{
  // Function that subtracts a matrix from another matrix
  Matrix eval;
  for(int i = 0; i < N_DIM; i++){
    for(int j = 0; j < N_DIM; j++){
      eval[i][j] = m1[i][j]-m2[i][j];
    }
  }
  return eval;
}

LinAlg3x3::Matrix LinAlg3x3::matrix_Identity()
{
  // Function that returns the identity matrix
  Matrix eval;
  for(int i = 0; i < N_DIM; i++){
    for(int j = 0; j < N_DIM; j++){
      if(i == j) eval[i][j] = 1.0;
      else eval[i][j] = 0.0;
    }
  }
  return eval;
}

void LinAlg3x3::matrix_mult_inplace(Matrix& m1, Real scalar)
{
  // Function that multiplies a matrix by a scalar in place (modifies the matrix in place)
  for(int i = 0; i < N_DIM; i++){
    for(int j = 0; j < N_DIM; j++){
      m1[i][j] = m1[i][j]*scalar;
    }
  }
  return;
}

void LinAlg3x3::matrix_zero_inplace(Matrix& m1)
{
	// Function that zeros out a matrix
	m1.fill({});
}

LinAlg3x3::Matrix LinAlg3x3::matrix_dot(const Matrix& Qt1,const Matrix& Qt2){	
	// Function that performs matrix matrix dot product of 3x3 matrices
	//
	// Args:
	// 	Qt1(LinAlg3x3::Matrix -> std::array<std::array<double,3>,3>: Q tensor 1
	// 	Qt2(LinAlg3x3::Matrix -> std::array<std::array<double,3>,3>: Q tensor 2
	//
	// returns:
	// 	res(LinAlg3x3::Matrix -> std::array<std::array<double,3>,3>: Q1*Q2
	
	
	// Initialize matrix to a matrix of 0's
	Matrix res;	
	// default initializes the array
	res.fill({});
		
	// perform matrix matrix multiplication (simpliest N^3 algorithm but it's only 27 = 3x3x3)
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			for(int k=0;k<3;k++){
				res[i][j] += Qt1[i][k] * Qt2[k][j]; 
			}
		}	
	}
	return res;
}

std::pair<LinAlg3x3::Real3,LinAlg3x3::Matrix> LinAlg3x3::EigenSolver(const Matrix& Q){
	Eigen::Map<Eigen::Matrix3f> Qeigen(const_cast<Real*>(Q[0].data()),3,3);
	Eigen::EigenSolver<Eigen::Matrix3f> eigensolver;
	eigensolver.compute(Qeigen);

	Eigen::Vector3f eigenvalue = eigensolver.eigenvalues().real();
	Eigen::Matrix3f eigenvec   = eigensolver.eigenvectors().real();

    // copy eigenvalue
	Real3 eigval;
	for (int i=0;i<3;i++)
	{
		eigval[i] = eigenvalue[i];
	}

    // copy eigenvector
    Matrix eigvec;
    for (int i=0;i<N_DIM;i++)
    {
        for (int j=0;j<N_DIM;j++)
        {
            eigvec[i][j] = eigenvec(i,j);
        }
    }

    // initialize return value
    std::pair<Real3, Matrix> ret;
    ret.first = eigval;
    ret.second = eigvec;

    return ret;
}

std::pair<LinAlg3x3::Real3, LinAlg3x3::Matrix> LinAlg3x3::OrderEigenSolver(const Matrix& Q)
{
  auto res = EigenSolver(Q);
	std::vector<size_t> order = LinAlg3x3::argsort(res.first);
	std::pair<Real3,Matrix> ordered_eig;

    // find the ordered eigenvector
	for(int i=0;i<3;i++){
		ordered_eig.first[i] = res.first[order[i]];
		for(int j=0;j<3;j++){
			ordered_eig.second[j][i] = res.second[j][order[i]];
		}
	}

 	return ordered_eig;	
}

std::vector<size_t> LinAlg3x3::argsort(const Real3& array) {
    std::vector<size_t> indices(array.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(),
              [&array](int left, int right) -> bool {
                  // sort indices according to corresponding array element
                  return array[left] > array[right];
              });

    return indices;
}