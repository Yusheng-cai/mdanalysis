#pragma once

#include "tools/CommonTypes.h"

#include <chrono>
#include <iostream>
#include <string>
#include <map>

class Timer
{
    using Real = CommonTypes::Real;

    public:
        void start();
        void end();
        Real diff();
        
    private:
        std::chrono::high_resolution_clock::time_point s_;
        std::chrono::high_resolution_clock::time_point e_;
};