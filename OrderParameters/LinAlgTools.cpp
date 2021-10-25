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

LinAlg3x3::Real LinAlg3x3::DotProduct(const Real3& v1, const Real3& v2)
{
    Real ret = 0.0;

    for (int i=0;i<3;i++)
    {
        ret += v1[i]*v2[i];
    }

    return ret;
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

    for (int i=0;i<3;i++)
    {
        for (int j=0;j<3;j++)
        {
            std::cout << rotMat[i][j] << "\t";
        }
        std::cout << "\n";
    }

    return rotMat;
}

LinAlg3x3::Matrix LinAlg3x3::GetRotationMatrix(const Real3& v1, const Real3& v2)
{
    Real3 crossProduct = LinAlg3x3::CrossProduct(v1, v2);
    Real norm = LinAlg3x3::norm(crossProduct);
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