#pragma once 

#include "CommonTypes.h"

#include <iostream>

using Real3  = CommonTypes::Real3;
using Real   = CommonTypes::Real;
using Matrix = CommonTypes::Matrix;


inline std::ostream& operator<<(std::ostream &out, const Real3& v) {
    return out << v[0] << ' ' << v[1] << ' ' << v[2];
}

inline Real3 operator+(const Real3& v1, const Real3& v2)
{
    return {{v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]}};
}

inline Real3 operator*(const Real3& v1, const Real3& v2)
{
    return {{v1[0] * v2[0], v1[1] * v2[1], v1[2] * v2[2]}};
}

inline Real3 operator/(const Real3& v1, const Real3& v2)
{
    return {{v1[0]/v2[0], v1[1]/v2[1], v1[2]/v2[2]}};
}

inline Real3 operator+(const Real3& v1, Real value)
{
    return {{v1[0] + value, v1[1] + value, v1[2] + value}};
}

inline Real3 operator*(const Real3& v1, Real value)
{
    return {{v1[0] * value, v1[1] * value, v1[2] * value}};
}

inline Real3 operator/(const Real3& v1, Real value)
{
    return {{v1[0]/value, v1[1]/value, v1[2]/value}};
}

inline std::ostream& operator<<(std::ostream& out, const Matrix& m)
{
    for (int i=0;i<3;i++)
    {
        for (int j=0;j<3;j++)
        {
            out << m << " ";
        }
        out << "\n";
    }

    return out;
}

inline Matrix operator*(const Matrix& mat1, Real value)
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

inline Matrix operator+(const Matrix& mat1, Real value)
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

inline Matrix operator+(const Matrix& mat1, const Matrix& mat2)
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