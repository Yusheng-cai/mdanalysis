#pragma once 

#include <iostream>
#include <vector>
#include <array>

namespace CommonTypes{
    using Real          = double;
    using VectorReal    = std::vector<Real>;
    using Real2         = std::array<Real,2>;
    using Real3         = std::array<Real,3>;
    using Matrix2       = std::array<Real2,2>;
    using Matrix        = std::array<Real3,3>;
    using VecofVecReal2 = std::vector<std::vector<Real2>>;
    using index2        = std::array<int,2>;
    using index3        = std::array<int,3>;
    using VectorReal3   = std::vector<Real3>;

    template<typename T>
    using VectorofVector= std::vector<std::vector<T>>;
};