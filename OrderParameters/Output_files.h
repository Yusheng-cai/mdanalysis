#pragma once

#include "tools/InputParser.h"
#include "tools/Assert.h"
#include "tools/CommonTypes.h"

#include <string>
#include <iostream>
#include <iomanip>


class Driver;

class OutputStream
{
    public:
        using Real   = CommonTypes::Real;

        OutputStream(const ParameterPack& input, const Driver& sim);
        ~OutputStream(){};

        // Function that opens the output file
        void Open();
        void close();

        // Helper functions
        void write_header();

        // Check if the file is open
        bool isOpen(){return ofs_.is_open();};

        // Check if the simulation is on step with stride
        bool onStep();
        void printIfOnStep();

    private:
        std::string name_ = "op.out";
        int stride_ = 1;
        std::vector<std::string> value_;

        const Driver& driver_;

        int precision_ = 3;

        // the output stream 
        std::ofstream ofs_;

        bool is_restart_ = false;
};