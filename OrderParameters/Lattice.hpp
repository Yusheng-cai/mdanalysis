#pragma once

#include "tools/CommonTypes.h"
#include "tools/CommonOperations.h"

#include <vector>
#include <array>

template <typename T>
class Lattice
{
    public:
        using INT3 = std::array<int,3>;
        using INT2 = std::array<int,2>;

        Lattice(int nx, int ny, int nz);
        Lattice(int nx, int ny, int nz, T num);
        Lattice() {nx_=0;ny_=0;nz_=0;total_size_=0;shape_={{nx_,ny_,nz_}};};
        Lattice(INT3 index);
        Lattice(INT3 index, T num);

        T& operator()(int a, int b, int c);
        T& operator()(INT3 index3);

        void resize(int nx, int ny, int nz);
        void resize(int nx, int ny, int nz, T value);
        void resize(INT3 shape);
        void resize(INT3 shape, T value);

        void fill(T value);

        INT3 getIndex3(int index);
        int getIndex(int a, int b, int c);
        INT3 getShape() const {return shape_;}
        int getSize() const {return total_size_;}

        Lattice operator* (T x);
        Lattice operator+ (T x);
        Lattice& operator+=(const Lattice<T>& other);
        Lattice& operator-=(const Lattice<T>& other);

        T& operator[] (int num);
        T operator[] (int num) const;

        // reduction 
        std::vector<std::vector<T>> reduce(int reduce_dim);

    private:
        std::vector<T> lattice_;
        int nx_,ny_,nz_;
        int total_size_;
        INT3 shape_;
};

template <typename T>
std::vector<std::vector<T>> Lattice<T>::reduce(int reduce_dim){
    int dim1,dim2;
    if (reduce_dim == 0){dim1=ny_;dim2=nz_;}
    else if (reduce_dim == 1){dim1=nx_;dim2=nz_;}
    else {dim1=nx_;dim2=ny_;}
}

template <typename T> 
Lattice<T>::Lattice(int nx, int ny , int nz)
:nx_(nx), ny_(ny), nz_(nz)
{
    total_size_ = nx_ * ny_ * nz_;
    lattice_.clear();
    lattice_.resize(total_size_);
    shape_ = {{nx_,ny_,nz_}};
}

template <typename T>
Lattice<T>::Lattice(int nx, int ny, int nz, T num)
:nx_(nx), ny_(ny), nz_(nz)
{
    total_size_ = nx_ * ny_ * nz_;
    lattice_.clear();
    lattice_.resize(total_size_, num);
    shape_ = {{nx_,ny_,nz_}};
}

template <typename T>
Lattice<T>::Lattice(INT3 index, T num)
:Lattice(index[0], index[1], index[2], num)
{
}

template <typename T>
Lattice<T>::Lattice(INT3 index)
:Lattice(index[0],index[1], index[2])
{
}

template <typename T>
T& Lattice<T>::operator()(int a, int b, int c)
{
    int index = getIndex(a,b,c);
    ASSERT((index < total_size_), "Index " << a << " " << b << " " << c << " out of range.");
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
    if (a < 0){a += nx_;}else{a %= nx_;}
    if (b < 0){b += ny_;}else{b %= ny_;}
    if (c < 0){c += nz_;}else{c %= nz_;}

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
    Lattice<T> newLat(total_size_);
    for (int i=0;i<total_size_;i++){
        newLat[i] = lattice_[i] * num;
    }
    return newLat;
}

template <typename T>
Lattice<T> Lattice<T>::operator+(T num)
{
    Lattice<T> newLat(total_size_);
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
    lattice_.clear();
    resize(shape[0], shape[1], shape[2]);
}

template <typename T>
void Lattice<T>::resize(INT3 shape, T value)
{
    lattice_.clear();
    resize(shape[0], shape[1], shape[2], value);
}

template <typename T>
void Lattice<T>::fill(T value){
    ASSERT((lattice_.size() != 0), "Trying to fill an empty lattice.");
    std::fill(lattice_.begin(), lattice_.end(), value);
}

template <typename T>
T& Lattice<T>::operator[](int index)
{
    return lattice_[index];
}

template <typename T>
T Lattice<T>::operator[](int index) const {
    return lattice_[index];
}

template <typename T>
Lattice<T>& Lattice<T>::operator+=(const Lattice<T>& other){
    for (int i=0;i<total_size_;i++){
        lattice_[i] = lattice_[i] + other[i];
    }

    return *this;
}

template <typename T>
Lattice<T>& Lattice<T>::operator-=(const Lattice<T>& other){
    for (int i=0;i<total_size_;i++){
        lattice_[i] = lattice_[i] - other[i];
    }

    return *this;
}
