#pragma once

#include "tools/Assert.h"

#include <iostream>
#include <vector>
#include <math.h>
#include <type_traits>


template<typename T>
class Matrix;

template<typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& mat);

// a specialized class for 2d Arrays
template<typename T>
class Matrix
{
    public:
        using Iterator = typename std::vector<T>::iterator;
        using cIterator = typename std::vector<T>::const_iterator;

        Matrix(){};
        Matrix(std::size_t NR, std::size_t NC);
        Matrix(std::size_t NR, std::size_t NC, T num);
        ~Matrix(){};

        // copy constructor 
        Matrix(const Matrix& other);

        // indexing operator
        T& operator()(int i, int j);
        // return the address of the first element of the ith row
        T* operator()(int i);
        const T& operator()(int i, int j) const;
        const T* operator()(int i) const;

        // addition or subtraction operators
        Matrix operator+(const Matrix& other) const;
        Matrix operator+(T other) const;

        Matrix operator-(const Matrix& other) const;
        Matrix operator-(T other) const;

        Matrix operator*(const Matrix& other) const;
        Matrix operator*(T other) const;

        Matrix operator/(const Matrix& other) const;
        Matrix operator/(T other) const;

        Matrix& operator=(const Matrix& other);

        Matrix<int> operator==(const Matrix& other) const;
        Matrix<int> operator==(T other) const;

        T sum();

        void resize(std::size_t size);
        void resize(std::size_t NR, std::size_t NC);

        // overload the << operator for std::cout
        friend std::ostream& operator<< <>(std::ostream& os, const Matrix& mat);

        // size
        int size() const{return NR_*NC_;}
        int getNR() const{return NR_;}
        int getNC() const{return NC_;}

        Iterator begin() {return data_.begin();}
        Iterator end() {return data_.end();}
        cIterator cbegin() {return data_.cbegin();}
        cIterator cend() {return data_.cend();}

        // returns the underlying pointer to the data (The contiguous data pointer)
        T* data(){return data_.data();}

        void AssertValidSize(int NR, int NC) const
        {
            ASSERT((NR > 0 && NC>0),"The size is not valid.");
        }

        void AssertSameSize(const Matrix& other) const
        {
            ASSERT((other.getNR() == NR_ && other.getNC() == NC_), "The size is not compatible.");
        }

    private:
        std::vector<T> data_;
        int NR_;
        int NC_;
};

template<typename T>
void Matrix<T>::resize(const std::size_t size)
{
    data_.resize(size);

    for (int i=0;i<size;i++)
    {
        data_[i] = 0;
    }

    // if you resize with one dimension, it will be assumed that NR = dimension
    NR_ = size;
    NC_ = 1;
}

template<typename T>
void Matrix<T>::resize(std::size_t NR, std::size_t NC)
{
    AssertValidSize(NR, NC);

    data_.resize(NR*NC);

    std::fill(data_.begin(), data_.end(), 0);
    NR_ = NR;
    NC_ = NC;
}

template<typename T>
Matrix<T>::Matrix(std::size_t NR, std::size_t NC)
{
    AssertValidSize(NR,NC);
    NR_ = NR;
    NC_ = NC;

    ASSERT((std::is_arithmetic<T>()), "The type is not arithmetic.");
    data_.resize(NR*NC);  
}

template<typename T>
Matrix<T>::Matrix(std::size_t NR, std::size_t NC, T num)
{
    AssertValidSize(NR, NC);
    NR_ = NR;
    NC_ = NC;

    data_.resize(NR*NC); 
}

template<typename T>
Matrix<T>::Matrix(const Matrix& other)
{
    AssertValidSize(other.getNR(), other.getNC());
    NR_ = other.getNR();
    NC_ = other.getNC();

    data_.resize(NR_*NC_);

    for (int i=0;i<NR_*NC_;i++)
    {
        data_[i] = other.data_[i];
    } 
}


template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix& other) const
{
    AssertValidSize(NR_,NC_);
    AssertSameSize(other);

    Matrix<T> output(NR_, NC_);
    for (int i=0;i<NR_;i++)
    {
        for (int j=0;j<NC_;j++)
        {
            output(i,j) = this->operator()(i,j) + other(i,j);
        }
    }

    return output;
}

template<typename T>
Matrix<T> Matrix<T>::operator+(T other) const
{
    AssertValidSize(NR_, NC_);
    AssertSameSize(other);

    Matrix<T> output(NR_, NC_);

    for (int i=0;i<NR_;i++)
    {
        for (int j=0;j<NC_;j++)
        {
            output(i,j) = this->operator()(i,j) + other;
        }
    }

    return output;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix& other) const
{
    AssertValidSize(NR_, NC_);
    AssertSameSize(other);
    Matrix<T> output(NR_, NC_);

    for (int i=0;i<NR_;i++)
    {
        for (int j=0;j<NC_;j++)
        {
            output(i,j) = this->operator()(i,j) - other(i,j);
        }
    }

    return output;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(T other) const
{
    AssertValidSize(NR_, NC_);
    AssertSameSize(other);

    Matrix<T> output(NR_, NC_);

    for (int i=0;i<NR_;i++)
    {
        for (int j=0;j<NC_;j++)
        {
            output(i,j) = this->operator()(i,j) - other;
        }
    }

    return output;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix& other) const
{
    AssertValidSize(NR_, NC_);
    AssertSameSize(other);

    Matrix<T> output(NR_, NC_);

    for (int i=0;i<NR_;i++)
    {
        for (int j=0;j<NC_;j++)
        {
            output(i,j) = this->operator()(i,j) * other(i,j);
        }
    }

    return output;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(T other) const
{
    AssertValidSize(NR_, NC_);
    AssertSameSize(other);

    Matrix<T> output(NR_, NC_);

    for (int i=0;i<NR_;i++)
    {
        for (int j=0;j<NC_;j++)
        {
            output(i,j) = this->operator()(i,j) * other;
        }
    }

    return output;
}

template<typename T>
Matrix<T> Matrix<T>::operator/(const Matrix& other) const
{
    AssertValidSize(NR_, NC_);
    AssertSameSize(other);

    Matrix<T> output(NR_, NC_);

    for (int i=0;i<NR_;i++)
    {
        for (int j=0;j<NC_;j++)
        {
            output(i,j) = this->operator()(i,j)/other(i,j);
        }
    }

    return output;
}

template<typename T>
Matrix<T> Matrix<T>::operator/(T other) const
{
    AssertValidSize(NR_, NC_);
    AssertSameSize(other);

    Matrix<T> output(NR_, NC_);

    for (int i=0;i<NR_;i++)
    {
        for (int j=0;j<NC_;j++)
        {
            output(i,j) = this->operator()(i,j)/other;
        }
    }

    return output;
}

template<typename T>
Matrix<int> Matrix<T>::operator==(const Matrix& other) const
{
    AssertValidSize(NR_, NC_);
    AssertSameSize(other);

    Matrix<int> output(NR_, NC_);
    for (int i=0;i<NR_;i++)
    {
        for (int j=0;j<NC_;j++)
        {
            if (this->operator()(i,j) != other(i,j))
            {
                output(i,j) = 0;
            }
            else
            {
                output(i,j) = 1;
            }
        }
    }

    return output;
}

template<typename T>
Matrix<int> Matrix<T>::operator==(T other) const
{
    AssertValidSize(NR_, NC_);
    AssertSameSize(other);

    Matrix<int> output(NR_, NC_);

    for (int i=0;i<NR_;i++)
    {
        for (int j=0;j<NC_;j++)
        {
            if (this->operator()(i,j) != other ) 
            {
                output(i,j) = 0;
            }
            else
            {
                output(i,j) = 1;
            }
        }
    }

    return output;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& mat)
{
    int NR = mat.getNR();
    int NC = mat.getNC();
    mat.AssertValidSize(NR, NC);

    for (int i=0;i<NR;i++)
    {
        for (int j=0;j<NC;j++)
        {
            os << mat(i,j) << " ";
        }
        os << "\n";
    }

    return os;
}

template<typename T>
T Matrix<T>::sum()
{
    AssertValidSize(NR_, NC_);

    T sum = 0;
    for (int i=0;i<NR_;i++)
    {
        for (int j=0;j<NC_;j++)
        {
            sum += this->operator()(i,j);
        }
    }

    return sum;
}

template<typename T>
T& Matrix<T>::operator()(int i, int j)
{
    AssertValidSize(NR_, NC_);
    ASSERT((i >= 0 && i < NR_), "Rows either exceeded max or less than minimum.");
    ASSERT((j >= 0 && j < NC_), "Columns either exceeded max or less than minimum.");

    int num = i*NC_ + j;

    return data_[num];
}

template<typename T>
const T& Matrix<T>::operator()(int i, int j) const
{
    AssertValidSize(NR_, NC_);
    ASSERT((i >= 0 && i < NR_), "Rows either exceeded max or less than minimum.");
    ASSERT((j >= 0 && j < NC_), "Columns either exceeded max or less than minimum.");

    int num = i*NC_ + j;

    return data_[num];
}

template<typename T>
T* Matrix<T>::operator()(int i)
{
    AssertValidSize(NR_, NC_);
    ASSERT((i >= 0 && i<NR_), "The specified row is larger than the Matrix.");

    int num = i*NC_;

    return &data_[num];
}

template<typename T>
const T* Matrix<T>::operator()(int i) const
{
    AssertValidSize(NR_, NC_);
    ASSERT((i >= 0 && i<NR_), "The specified row is larger than the Matrix.");

    int num = i*NC_;

    return &data_[num];
}

template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix& other)
{
    AssertValidSize(NR_,NC_);
    AssertSameSize(other);

    // copy data 
    for (int i=0;i<NR_*NC_;i++)
    {
        data_[i] = other.data_[i];
    }

    return *this;
}