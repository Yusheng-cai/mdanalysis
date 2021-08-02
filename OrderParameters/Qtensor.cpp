#include "Qtensor.h"

Qtensor::Real3 Qtensor::vec_add(const Qtensor::Real3& v1, const Qtensor::Real3& v2)
{
  // Function that sum a vector with another vector
  Real3 eval;
  for(int i = 0; i < N_DIM; i++){
    eval[i] = v1[i]+v2[i];
  }
  return eval;
}

Qtensor::Real3 Qtensor::vec_add(const Real3& v1, const Real s){
	Real3 eval;
	for(int i=0;i<3;i++){
		eval[i] = v1[i] + s;
	}

	return eval;
}

Qtensor::Real3 Qtensor::vec_add(const Real s, const Real3& v1){
	return Qtensor::vec_add(v1,s);
}

Qtensor::Real3 Qtensor::vec_sub(const Real3& v1, const Real3& v2)
{ 
  Real3 eval;
  for(int i = 0; i < N_DIM; i++){
    eval[i] = v1[i]-v2[i];
  }
  return eval;
}


Qtensor::Real3 Qtensor::vec_mult(const Real3& v1, const Real3& v2)
{
  Real3 eval;
  for(int i = 0; i < N_DIM; i++){
    eval[i] = v1[i]*v2[i];
  }
  return eval;
}

Qtensor::Real3 Qtensor::vec_mult(const Real3& v1, const Real scalar_value)
{
  Real3 eval;
  for(int i = 0; i < N_DIM; i++){
    eval[i] = v1[i]*scalar_value;
  }
  return eval;
}

Qtensor::Real3 Qtensor::vec_mult(const Real scalar_value, const Real3& v1){
	return Qtensor::vec_mult(v1,scalar_value);
}

Qtensor::Real Qtensor::vec_dot(const Real3& v1, const Real3& v2){
	Real dot_product = 0.0;
	for(int i=0;i<3;i++){
		dot_product += v1[i]*v2[i];
	}

	return dot_product;
}

Qtensor::Real Qtensor::vec_norm(Real3 v1)
{
  Real eval=Qtensor::vec_dot(v1,v1);

  return std::sqrt(eval);
}

Qtensor::Matrix Qtensor::vec_dyadic(const Real3& v1, const Real3& v2)
{
  Matrix eval;
  for(int i = 0; i < N_DIM; i++){
    for(int j = 0; j < N_DIM; j++){
    eval[i][j] = v1[i]*v2[j];
    }
  }
  return eval;
}

Qtensor::Real3 Qtensor::vec_cross(const Real3& v1, const Real3& v2){
	Real3 output;
	output.fill(0);

	output[0] 	= v1[1]*v2[2] - v1[2]*v2[1];	
	output[1]   = v1[2]*v2[0] - v1[0]*v2[2];
	output[2]   = v1[0]*v2[1] - v1[1]*v2[0];  

	return output;
}

Qtensor::Real Qtensor::matrix_trace(const Matrix& Q){
	// Function that returns the trace of a 3x3 matrix (Q tensor)
	// 
	// Args:
	// Q (Qtensor::Matrix): 3x3 matrix of the Q tensor
	//
	// returns:
	// trace (Qtensor::Real): the trace of the Q tensor
	
	
	Real trace = 0.0;
	for(int i=0;i<3;i++){
		trace = trace + Q[i][i]; 
	}
	return trace;
}

Qtensor::Matrix Qtensor::matrix_add(const Matrix& m1, const Matrix& m2)
{
  Matrix eval;
  for(int i = 0; i < N_DIM; i++){
    for(int j = 0; j < N_DIM; j++){
      eval[i][j] = m1[i][j]+m2[i][j];
    }
  }
  return eval;
}

void Qtensor::matrix_accum_inplace(Matrix& m1, const Matrix& m2)
{
  for(int i = 0; i < N_DIM; i++){
    for(int j = 0; j < N_DIM; j++){
      m1[i][j] = m1[i][j]+m2[i][j];
    }
  }
  return;
}
Qtensor::Matrix Qtensor::matrix_sub(const Matrix& m1, const Matrix& m2)
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

Qtensor::Matrix Qtensor::matrix_Identity()
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

void Qtensor::matrix_mult_inplace(Matrix& m1, Real scalar)
{
  // Function that multiplies a matrix by a scalar in place (modifies the matrix in place)
  for(int i = 0; i < N_DIM; i++){
    for(int j = 0; j < N_DIM; j++){
      m1[i][j] = m1[i][j]*scalar;
    }
  }
  return;
}

void Qtensor::matrix_zero_inplace(Matrix& m1)
{
	// Function that zeros out a matrix
	m1.fill({});
}

Qtensor::Matrix Qtensor::matrix_dot(const Matrix& Qt1,const Matrix& Qt2){	
	// Function that performs matrix matrix dot product of 3x3 matrices
	//
	// Args:
	// 	Qt1(Qtensor::Matrix -> std::array<std::array<double,3>,3>: Q tensor 1
	// 	Qt2(Qtensor::Matrix -> std::array<std::array<double,3>,3>: Q tensor 2
	//
	// returns:
	// 	res(Qtensor::Matrix -> std::array<std::array<double,3>,3>: Q1*Q2
	
	
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

Qtensor::Real3 Qtensor::eigval_Qtensor(const Matrix& Q){
	// First check if the matrix is diagonal or if the matrix is the 0 matrix
	Real sum = pow(Q[0][1],2) + pow(Q[0][2],2) + pow(Q[1][2],2);
	Real3 eigenvalue;
	Real tolerance=1e-8;
	eigenvalue.fill(0);

	if( sum <= tolerance){
		// Eigenvalues are the diagonals of the matrix if matrix is diagonal
		eigenvalue = {Q[0][0],Q[1][1],Q[2][2]};
	}
	else{
		// make matrix traceless first
		Matrix tempQ = Q;
		Real trQ = Qtensor::matrix_trace(Q);
		for(int i=0;i<3;i++){
			tempQ[i][i] -= trQ/3.0;
		}
		

		// follow the algorithm of traceless matrix from now on
		// Find Q^2 and Q^3
		Matrix Q2 = Qtensor::matrix_dot(tempQ,tempQ);
		Matrix Q3 = Qtensor::matrix_dot(Q2,tempQ);
		
		// Find tr(Q^2) and tr(Q^3)
		Real trQ2 = Qtensor::matrix_trace(Q2);
		Real trQ3 = Qtensor::matrix_trace(Q3);
		
		// find A and phi in the algorithm
		Real A = 2.0 * std::sqrt(trQ2/6.0);
		Real phi = std::acos((2.0*trQ3)/(A*trQ2));
		
		// analytical formula for traceless symmetric matrices for eigenvalue
		eigenvalue[0] = A*std::cos(phi/3) + trQ/3.0;
		eigenvalue[1] = A*std::cos((phi+2*Qtensor::PI)/3.0) + trQ/3.0;
		eigenvalue[2] = trQ -eigenvalue[1] - eigenvalue[0];
	}

	return eigenvalue;
}



Qtensor::Matrix Qtensor::eigvec_Qtensor(const Matrix& Q, const Real3& eigval){	
	Matrix eigenvector;	

	// first check if the matrix is diagonal
	Real sum = Q[0][1] + Q[0][2] + Q[1][2];	

	if(sum == 0){
		// if matrix is diagonal or the zero matrix, then the eigenvector is the identity matrix
		for(int i=0;i<3;i++){
			eigenvector[i][i] = 1.0;
		}
	}
	else{
		// This algorithm is adapted from Cayley-Hamilton theorem
		// please refer to https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3.C3.973_matrices
		Matrix Q1,Q2,Q3;
		Q1 = Q;
		Q2 = Q;
		Q3 = Q;
		
		// find (Q - lambda_i*I)
		for(int i=0;i<3;i++){
			Q1[i][i] -= eigval[0];
			Q2[i][i] -= eigval[1];
			Q3[i][i] -= eigval[2];
		}
		
		// Calculate the matrix in which the columns must be the corresponding eigenvectors
		Matrix e1 = Qtensor::matrix_dot(Q2,Q3);
		Matrix e2 = Qtensor::matrix_dot(Q1,Q3);
		Matrix e3 = Qtensor::matrix_dot(Q1,Q2);	
		std::vector<Matrix> evec = {e1,e2,e3};
		
		// Calculate the norm of each column of the matrix
		Matrix norm;

		// initializing norm matrix to be 0
		norm.fill({});
			
		for(int j=0;j<3;j++){
			// sum over rows
			for(int i=0;i<3;i++){
				// iterate over number of eigenvectors
				for(int num=0;num<3;num++){
					norm[num][j] += std::pow(evec[num][i][j],2.0);
				}
			}
		}
		
		// sqrt the norms
		for(int j=0;j<3;j++){	
			for(int num=0;num<3;num++){
				norm[num][j] = std::sqrt(norm[num][j]);
			}
		}
		
		// normalize the eigenvectors
		for(int j=0;j<3;j++){
			for(int i=0;i<3;i++){
				for(int num=0;num<3;num++){
					evec[num][i][j] = evec[num][i][j]/norm[num][j];
				}
			}
		}

		// write to eigenvector
		for(int num=0;num<3;num++){
			for(int j=0;j<3;j++){
				if(norm[num][j] != 0){
					for(int k=0;k<3;k++){
						eigenvector[k][num] = evec[num][k][j];
					}
					break;
				}
			}
		}	
	}

	return eigenvector;
}

std::pair<Qtensor::Real3,Qtensor::Matrix> Qtensor::eig_Qtensor(const Matrix& Q){
	Real3 eigval       = Qtensor::eigval_Qtensor(Q);
	Matrix eigenvector = Qtensor::eigvec_Qtensor(Q,eigval);
	
		
	// Create eigenvector/eigenvalue struct 
	std::pair<Real3,Matrix> e;
	e.first = eigval;
	e.second = eigenvector;

	return e;
}

std::pair<Qtensor::Real,Qtensor::Real3> Qtensor::OP_Qtensor(const Matrix& Q){
	std::pair<Real3,Matrix> e = Qtensor::eig_Qtensor(Q);	
	Real3 eigval = e.first;
	Matrix eigvec = e.second;
	std::pair<Real,Real3> p;

	std::vector<size_t> order = Qtensor::argsort(eigval);
	Real P2 = -2.0*eigval[order[1]];
	Real3 v1;

	for(int i=0;i<3;i++){
		v1[i] = eigvec[i][order[1]];
	}

	p.first = P2;
	p.second = v1;

 	return p;	
}

std::pair<Qtensor::VectorReal3, Qtensor::Real3> Qtensor::orderedeig_Qtensor(const Matrix& Q){
	std::pair<Real3,Matrix> eig = Qtensor::eig_Qtensor(Q);	
	Real3 eigval  = eig.first;
	Matrix eigvec = eig.second;
	

	std::vector<size_t> order = Qtensor::argsort(eigval);
	std::pair<VectorReal3,Real3> ordered_eig;

	for(int i=0;i<3;i++){
		ordered_eig.second[i] = eigval[order[i]];
		Real3 temp_eigvec;
		for(int j=0;j<3;j++){
			temp_eigvec[j] = eigvec[j][order[i]];
		}
		ordered_eig.first.push_back(temp_eigvec);
	}

 	return ordered_eig;	
}

std::vector<size_t> Qtensor::argsort(const Real3& array) {
    std::vector<size_t> indices(array.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(),
              [&array](int left, int right) -> bool {
                  // sort indices according to corresponding array element
                  return array[left] > array[right];
              });

    return indices;
}

Qtensor::Real3 Qtensor::normalize_director(const Real3& director){
	// Function that returns a normalized director bond	
	Real3 normalized;
	
	Real norm = 0.0;
	for(int i=0;i<3;i++){
		norm += pow(director[i],2.0);
	}

	norm = std::sqrt(norm);

	for(int i=0;i<3;i++){
		normalized[i] = director[i]/norm;
	}

	return normalized;
}