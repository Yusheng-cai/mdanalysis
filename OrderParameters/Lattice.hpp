#pragma once

#include "tools/CommonTypes.h"

#include <vector>
#include <array>

template <typename T>
class Lattice
{
    public:
        using INT3 = std::array<int,3>;

        Lattice(int nx, int ny, int nz);
        Lattice() {nx_=0;ny_=0;nz_=0;total_size_=0;shape_={{nx_,ny_,nz_}};};

        T& operator()(int a, int b, int c);
        T& operator()(INT3 index3);

        void resize(int nx, int ny, int nz);
        void resize(int nx, int ny, int nz, T value);
        void resize(INT3 shape);
        void resize(INT3 shape, T value);

        INT3 getIndex3(int index);
        int getIndex(int a, int b, int c);
        INT3 getShape() const {return shape_;}
        int getSize() const {return total_size_;}

        Lattice operator* (T x);
        Lattice operator+ (T x);

        T operator[] (int num);

    private:
        std::vector<T> lattice_;
        int nx_,ny_,nz_;
        int total_size_;
        INT3 shape_;
};

template <typename T> 
Lattice<T>::Lattice(int nx, int ny , int nz)
:nx_(nx), ny_(ny), nz_(nz)
{
    total_size_ = nx_ * ny_ * nz_;
    lattice_.resize(total_size_);
    shape_ = {{nx_,ny_,nz_}};
}

template <typename T>
T& Lattice<T>::operator()(int a, int b, int c)
{
    int index = getIndex(a,b,c);

    return lattice_[index];
}

template <typename T>
T& Lattice<T>::operator()(INT3 index)
{
    return operator()(index[0], index[1], index[2]);
}

template <typename T> 
int Lattice<T>::getIndex(int a, int b , int c)
{
    return a + b * nx_ + c * nx_ * ny_;
}

template <typename T>
typename Lattice<T>::INT3 Lattice<T>::getIndex3(int index)
{
    INT3 ind;
    ind[2] = index / (nx_ * ny_);
    ind[1] = (index - ind[2] * nx_ * ny_) / nx_;
    ind[0] = index - ind[1] * nx_ - ind[2] * nx_ * ny_;

    return ind;
}

template <typename T>
Lattice<T> Lattice<T>::operator*(T num)
{
    Lattice<T> newLat(total_size_,0.0);
    for (int i=0;i<total_size_;i++)
    {
        newLat[i] = lattice_[i] * num;
    }
    return newLat;
}

template <typename T>
Lattice<T> Lattice<T>::operator+(T num)
{
    Lattice<T> newLat(total_size_,0.0);
    for (int i=0;i<total_size_;i++)
    {
        newLat[i] = lattice_[i] + num;
    }

    return newLat;
}

template <typename T>
void Lattice<T>::resize(int nx, int ny, int nz, T value)
{
    resize(nx, ny, nz);
    std::fill(lattice_.begin(), lattice_.end(), value);
}

template <typename T>
void Lattice<T>::resize(int nx, int ny, int nz)
{
    nx_ = nx; ny_=ny; nz_=nz;

    total_size_=nx_*ny_*nz_;

    shape_={{nx_,ny_,nz_}};

    lattice_.resize(total_size_);
}

template <typename T>
void Lattice<T>::resize(INT3 shape)
{
    resize(shape[0], shape[1], shape[2]);
}

template <typename T>
void Lattice<T>::resize(INT3 shape, T value)
{
    resize(shape[0], shape[1], shape[2], value);
}

template <typename T>
T Lattice<T>::operator[](int index)
{
    return lattice_[index];
}

