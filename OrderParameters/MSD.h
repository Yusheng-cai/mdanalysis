#pragma once
#include "FFT.h"
#include "tools/CommonTypes.h"
#include "tools/Assert.h"
#include "Calculation.h"
#include "SimulationState.h"
#include "LinAlgTools.h"

#include <vector>
#include <string>
#include <map>
#include <array>

class MSD : public Calculation 
{
    public:
        using INT2 = std::array<int,2>;
        using Matrix = CommonTypes::Matrix;
        MSD(const CalculationInput& input);

        virtual void calculate() override;
        virtual void finishCalculate() override;

        // function that calculates the MSD for a single atom over time 
        // dimension : (nframes, ndim)
        std::vector<std::vector<Real>> calculateMSD(std::vector<Real3>& data);

        // output function
        void printMSD(std::string name);

    private:
        // (num_frames, num_residue, 3)
        std::vector<std::vector<Real3>> positions_;

        std::string residueName_;

        int numResidues_;

        std::string directionStr_="xyz";

        std::vector<std::string> VectorDirection_;

        std::vector<int> directionIndex_;

        std::vector<std::vector<int>> directionsIndex_;

        int directionSize=3;

        std::vector<std::vector<Real>> MSD_;

        std::map<std::string,std::vector<int>> MapStrToIndex2 = 
        {{"x", {0}}, 
         {"y", {1}}, 
         {"z", {2}}, 
         {"xz", {0,2}}, 
         {"xy", {0,1}},
         {"yz", {1,2}},
         {"xyz",{0,1,2}}};

         std::string RotationMode_;
         Real3 array_;
         bool UseDirector_=false;
         int headindex_=1;
         int tailindex_=2;
         bool rotate_=false;
         Matrix RotationMatrix_;
};