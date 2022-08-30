#pragma once 

#include "CommonTypes.h"

using Real3  = CommonTypes::Real3;
using Real   = CommonTypes::Real;
using Matrix = CommonTypes::Matrix;

Real3 operator+(const Real3& v1, const Real3& v2)
{
    Real3 ret;
    for (int i=0;i<3;i++)
    {
        ret[i] = v1[i] + v2[i];
    }

    return ret;
}

Real3 operator*(const Real3& v1, const Real3& v2)
{
    Real3 ret;
    for (int i=0;i<3;i++)
    {
        ret[i] = v1[i] * v2[i];
    }
    return ret;
}

Matrix operator*(const Matrix& mat1, Real value)
{
    Matrix ret;
    for (int i=0;i<3;i++)
    {
        for (int j=0;j<3;j++)
        {
            ret[i][j] = mat1[i][j] * value;
        }
    }

    return ret;
}

Matrix operator+(const Matrix& mat1, Real value)
{
    Matrix ret;
    for (int i=0;i<3;i++)
    {
        for (int j=0;j<3;j++)
        {
            ret[i][j] = mat1[i][j] + value;
        }
    }

    return ret;
}

Matrix operator+(const Matrix& mat1, const Matrix& mat2)
{
    Matrix ret;
    for (int i=0;i<3;i++)
    {
        for (int j=0;j<3;j++)
        {
            ret[i][j] = mat1[i][j] + mat2[i][j];
        }
    }

    return ret;
}