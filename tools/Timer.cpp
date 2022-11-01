#include "Timer.h"

void Timer::start(){
    s_ = std::chrono::high_resolution_clock::now();
}

void Timer::end(){
    e_ = std::chrono::high_resolution_clock::now();
}

Timer::Real Timer::diff(){
    auto diff = std::chrono::duration_cast<std::chrono::microseconds>(e_-s_);

    return diff.count();
}