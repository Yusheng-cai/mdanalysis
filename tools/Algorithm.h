#pragma once

#include "CommonTypes.h"

#include <algorithm>
#include <map>
#include <numeric>
#include <vector>
#include <random>

namespace Algorithm
{
    using Real = CommonTypes::Real;
    void Permutation(int max, int numSamples, int numTimes, std::vector<std::vector<int>>& samples);
    void Permutation(int max, int numSamples, std::vector<int>& samples);

    template <typename T>
    std::vector<int> argsort(std::vector<T>& vec, bool min=true);

    template <typename T, std::size_t dim>
    int argmin(std::array<T,dim>& arr);

    template <typename T, std::size_t dim>
    int argmax(std::array<T,dim>& arr);

    template <typename T>
    int argmin(std::vector<T>& vec);

    template <typename T>
    int argmax(std::vector<T>& vec);

    template <typename T> 
    std::vector<T> arange(T min, T max, T step);

    template <typename T>
    std::vector<T> linspace(T min, T max, int num);

    template <typename T>
    T max(std::vector<T>& vec);

    template <typename T>
    T min(std::vector<T>& vec);

    template <typename T, std::size_t dim>
    T max(std::array<T,dim>& arr);

    template <typename T>
    bool contain(std::vector<T>& vec, T num);

    template <typename T>
    bool contain(std::vector<T>& vec, T num, int& index);

    template <typename T>
    void unique(std::vector<T>& vec); 

    template <typename T, std::size_t dim>
    bool is_unique(std::array<T,dim>& arr);

    template <typename T, std::size_t dim>
    void sort(std::array<T, dim>& arr);

    template <typename key, typename value>
    bool FindInMap(const std::map<key,value>& map, const key& k, value& v);

    // insert something into map --> this value needs to not exist in the map previously
    template <typename key, typename value>
    bool InsertInMap(const std::map<key,value>& map, const key& k, const value& v);

    template <typename key, typename val>
    bool IsInMap(std::map<key, val>& map, key& k, val& v);
};

template<typename T>
std::vector<T> Algorithm::arange(T start, T stop, T step) {
    std::vector<T> values;
    for (T value = start; value < stop; value += step)
        values.push_back(value);
    return values;
}

template <typename T, std::size_t dim>
int Algorithm::argmin(std::array<T,dim>& arr)
{
    typename std::array<T,dim>::iterator it = std::min_element(arr.begin(), arr.end());

    return it - arr.begin();
}

template <typename T, std::size_t dim>
int Algorithm::argmax(std::array<T,dim>& arr)
{
    typename std::array<T,dim>::iterator it = std::max_element(arr.begin(), arr.end());

    return it - arr.begin();
}

template <typename T>
int Algorithm::argmax(std::vector<T>& vec)
{
    typename std::vector<T>::iterator it = std::max_element(vec.begin(), vec.end());

    return it - vec.begin();
}

template <typename T>
int Algorithm::argmin(std::vector<T>& vec)
{
    typename std::vector<T>::iterator it = std::min_element(vec.begin(), vec.end());

    return it - vec.begin();
}

template <typename T>
std::vector<T> Algorithm::linspace(T min, T max, int num)
{
    T step = (max - min)/num;
    std::vector<T> ret;
    for (int i=0;i<num;i++)
    {
        ret.push_back(min + i*step);
    }

    return ret;
}

template <typename T>
T Algorithm::max(std::vector<T>& vec)
{
    typename std::vector<T>::iterator max_element = std::max_element(vec.begin(), vec.end());

    return *max_element;
}

template <typename T, std::size_t dim>
T Algorithm::max(std::array<T,dim>& arr)
{
    typename std::array<T,dim>::iterator max_element = std::max_element(arr.begin(), arr.end());

    return *max_element;
}

template <typename T>
T Algorithm::min(std::vector<T>& vec)
{
    typename std::vector<T>::iterator min_element = std::min_element(vec.begin(), vec.end());

    return *min_element;
}

template <typename T>
bool Algorithm::contain(std::vector<T>& vec, T num)
{
    return (std::find(vec.begin(), vec.end(), num) != vec.end());
}

template <typename T>
bool Algorithm::contain(std::vector<T>& vec, T num, int& index)
{
    typename std::vector<T>::iterator it = std::find(vec.begin(), vec.end(), num);

    if (it == vec.end())
    {
        index = -1;
        return false;
    }
    else
    {
        index = it - vec.begin();
        return true;
    }
}

template <typename T>
void Algorithm::unique(std::vector<T>& vec)
{
    std::sort(vec.begin(), vec.end());
    typename std::vector<T>::iterator it = std::unique(vec.begin(), vec.end());
    vec.resize(std::distance(vec.begin(), it));
}

template <typename T, std::size_t dim>
bool Algorithm::is_unique(std::array<T,dim>& arr)
{
    std::array<T,dim> temp = arr;
    std::sort(temp.begin(), temp.end());
    typename std::array<T,dim>::iterator pos = std::adjacent_find(std::begin(temp), std::end(temp));
    if (pos != std::end(temp)){return false;}
    else{return true;}
}

template <typename T, std::size_t dim>
void Algorithm::sort(std::array<T,dim>& arr)
{
    std::sort(arr.begin(), arr.end());
}

template <typename key, typename value>
bool Algorithm::FindInMap(const std::map<key,value>& map, const key& k, value& v)
{
    typename std::map<key,value>::const_iterator it = map.find(k);
    if (it != map.end()){v=it->second;return true;}
    else{return false;}
}


template <typename key, typename value>
bool Algorithm::InsertInMap(const std::map<key,value>& map, const key& k, const value& v){
    typename std::map<key,value>::const_iterator it = map.find(k);  
    if (it != map.end()){map.insert(std::make_pair(k,v));}
    else{return false;}

    return true;
}

template <typename T>
std::vector<int> Algorithm::argsort(std::vector<T>& vec, bool min)
{
    std::vector<int> indices(vec.size());
    std::iota(indices.begin(), indices.end(), 0);

    if (min)
    {
        std::sort(indices.begin(), indices.end(),
                [&vec](int left, int right) -> bool {
                    // sort indices according to corresponding array element
                    return vec[left] < vec[right];
                });
    }
    else
    {
        std::sort(indices.begin(), indices.end(), 
                [&vec](int left, int right)-> bool {
                    return vec[left] > vec[right];
                });
    }

    return indices;
}


template <typename key, typename val>
bool Algorithm::IsInMap(std::map<key,val>& map, key& k, val& v)
{
    typename std::map<key,val>::iterator it = map.find(k);

    if (it == map.end())
    {
        return false;
    }

    v = it -> second;

    return true;
}