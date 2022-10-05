#pragma once 

#include "CommonTypes.h"
#include "Assert.h"

#include <iostream>
#include <array>
#include <vector>

using Real3  = CommonTypes::Real3;
using Real   = CommonTypes::Real;
using Matrix = CommonTypes::Matrix;

template <typename T, std::size_t dim>
inline std::ostream& operator<<(std::ostream& out, const std::array<T,dim>& arr){
    for (int i=0;i<dim;i++){
        out << arr[i] << " ";
    }

    return out;
}


template <typename T, std::size_t dim>
inline std::array<T,dim> operator+(const std::array<T,dim>& v1, const std::array<T,dim>& v2){
    std::array<T,dim> ret;
    for (int i=0;i<dim;i++){
        ret[i] = v1[i] + v2[i];
    }

    return ret;
}

template <typename T1, typename T2,std::size_t dim>
inline std::array<T1,dim> operator*(const std::array<T1,dim>& v1, const std::array<T2,dim>& v2)
{
    std::array<T1,dim> ret;
    for (int i=0;i<dim;i++){
        ret[i] = v1[i] * v2[i];
    }

    return ret;
}

template <typename T, std::size_t dim>
inline std::array<T,dim> operator-(const std::array<T,dim>& v1, const std::array<T,dim>& v2)
{
    std::array<T,dim> ret;
    for (int i=0;i<dim;i++){
        ret[i] = v1[i] - v2[i];
    }

    return ret;
}

template <typename T, std::size_t dim>
inline std::array<T,dim> operator/(const std::array<T,dim>& v1, const std::array<T,dim>& v2)
{
    std::array<T,dim> ret;
    for (int i=0;i<dim;i++){
        ret[i] = v1[i] / v2[i];
    }

    return ret;
}

template <typename T1, typename T2, std::size_t dim>
inline std::array<T1,dim> operator+(const std::array<T1,dim>& v1, T2 value)
{
    std::array<T1,dim> ret;
    for (int i=0;i<dim;i++){
        ret[i] = v1[i] + value;
    }

    return ret;
}

template <typename T1, typename T2, std::size_t dim>
inline std::array<T1,dim> operator-(const std::array<T1,dim>& v1, T2 value)
{
    std::array<T1,dim> ret;
    for (int i=0;i<dim;i++){
        ret[i] = v1[i] - value;
    }

    return ret;
}

template <typename T1, typename T2, std::size_t dim>
inline std::array<T1,dim> operator*(const std::array<T1,dim>& v1, T2 value)
{
    std::array<T1,dim> ret;
    for (int i=0;i<dim;i++){
        ret[i] = v1[i] * value;
    }

    return ret;
}

template <typename T1, typename T2, std::size_t dim>
inline std::array<T1,dim> operator/(const std::array<T1,dim>& v1, T2 value)
{
    std::array<T1,dim> ret;
    for (int i=0;i<dim;i++){
        ret[i] = v1[i] / value;
    }

    return ret;
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
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
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

/*
    Vectors
*/

template <typename T>
inline std::ostream& operator<<(std::ostream& out, const std::vector<T>& v){
    for (int i=0;i<v.size();i++){
        out << v[i] << " ";
    }

    return out;
}

template <typename T>
inline std::vector<T> operator+(const std::vector<T>& v1, const std::vector<T>& v2){
    ASSERT((v1.size() == v2.size()), \
    "Cannot add the 2 vectors, one is of size " << v1.size() << " while the other is of size " << v2.size());
    int dim = v1.size();
    std::vector<T> ret(dim);
    for (int i=0;i<dim;i++){
        ret[i] = v1[i] + v2[i];
    }

    return ret;
}

template <typename T>
inline std::vector<T> operator-(const std::vector<T>& v1, const std::vector<T>& v2){
    ASSERT((v1.size() == v2.size()), \
    "Cannot add the 2 vectors, one is of size " << v1.size() << " while the other is of size " << v2.size());
    int dim = v1.size();
    std::vector<T> ret(dim);
    for (int i=0;i<dim;i++){
        ret[i] = v1[i] - v2[i];
    }

    return ret;
}

template <typename T>
inline std::vector<T> operator*(const std::vector<T>& v1, const std::vector<T>& v2){
    ASSERT((v1.size() == v2.size()), \
    "Cannot add the 2 vectors, one is of size " << v1.size() << " while the other is of size " << v2.size());
    int dim = v1.size();
    std::vector<T> ret(dim);
    for (int i=0;i<dim;i++){
        ret[i] = v1[i] * v2[i];
    }

    return ret;
}

template <typename T>
inline std::vector<T> operator/(const std::vector<T>& v1, const std::vector<T>& v2){
    ASSERT((v1.size() == v2.size()), \
    "Cannot add the 2 vectors, one is of size " << v1.size() << " while the other is of size " << v2.size());
    int dim = v1.size();
    std::vector<T> ret(dim);
    for (int i=0;i<dim;i++){
        ret[i] = v1[i] / v2[i];
    }

    return ret;
}

template <typename T1, typename T2>
inline std::vector<T1> operator+(const std::vector<T1>& v1, T2 v){
    int dim = v1.size();
    std::vector<T1> ret(dim);
    for (int i=0;i<dim;i++){
        ret[i] = v1[i] + v;
    }

    return ret;
}

template <typename T1, typename T2>
inline std::vector<T1> operator-(const std::vector<T1>& v1, T2 v){
    int dim = v1.size();
    std::vector<T1> ret(dim);
    for (int i=0;i<dim;i++){
        ret[i] = v1[i] - v;
    }

    return ret;
}

template <typename T1, typename T2>
inline std::vector<T1> operator*(const std::vector<T1>& v1, T2 v){
    int dim = v1.size();
    std::vector<T1> ret(dim);
    for (int i=0;i<dim;i++){
        ret[i] = v1[i] * v;
    }

    return ret;
}

template <typename T1, typename T2>
inline std::vector<T1> operator/(const std::vector<T1>& v1, T2 v){
    int dim = v1.size();
    std::vector<T1> ret(dim);
    for (int i=0;i<dim;i++){
        ret[i] = v1[i] / v;
    }

    return ret;
}