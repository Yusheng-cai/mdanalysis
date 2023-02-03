#include "SlabVsr.hpp"

SlabVsr::SlabVsr(const CalculationInput& input)
: Calculation(input)
{
    vsr_ = SRE_ptr(new SRE(input));
}

void SlabVsr::calculate(){
    // first calculate vsr 
    vsr_->calculate();
}