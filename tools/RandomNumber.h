#pragma once

#include "Assert.h"
#include "CommonTypes.h"

#include <iostream>
#include <ctime>
#include <cmath>
#include <random>
#include <chrono>

class Random 
{
    public:
        using Real          = CommonTypes::Real;
        ~Random() = default;

        // delete copy constructor & assignment operator for singleton design
        Random(const Random& other) = delete;
        Random& operator=(const Random& other) = delete;

        static Random& instance(){return _instance;};

        Real DrawUniform();
        Real DrawExponential(Real lambda);
        Real DrawUniform_minmax(Real min, Real max); 

        static void seed();

        // Set the seed of the object 
        static void set_seed(unsigned int seed){
            instance().set_seedImpl(seed);
            instance().set_seedboolImpl(true);
            }

    private:
        Random(){distribution_ = std::uniform_real_distribution<Real>(0.0,1.0);}

        static Random _instance;
        std::uniform_real_distribution<Real> distribution_;

        unsigned seed_ = 0;
        bool user_defined_seed = false;

        void set_seedImpl(unsigned int seed){seed_ = seed;}
        void set_seedboolImpl(bool b){user_defined_seed = b;};
        void seedImpl();

        std::mt19937 generator_;

        std::random_device rd;
};