#pragma once 

#include "Calculation.h"
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "parallel/OpenMP.h"
#include "SimulationState.h"
#include "LinAlgTools.h"
#include "Bin.h"

#include <vector>
#include <string>
#include <array>
#include <memory>


// given that a molecule forms a angle costheta with respect to surface and a distance R
// with respect to another molecule, what is the average costheta   

class grcostorientation : public Calculation
{
    public:
        using binptr = std::unique_ptr<Bin>;
        using Range  = CommonTypes::Real2;

        grcostorientation(const CalculationInput& input);

        virtual void calculate();
        virtual void finishCalculate();

        void printHistogram2d(std::string name);
    private:
        int headindex_=1;
        int tailindex_=2;

        std::vector<std::vector<Real>> histogram2d_;
        std::vector<std::vector<Real>> histogram2dPerIter_;

        // The 1d histogram
        std::vector<Real> histogram_;

        // rbin 
        binptr rbin_;
        int numrbins_;

        // tbin
        binptr tbin_;
        int numtbins_=20;
        Range trange_={{-1,1}};

        // surface normal
        Real3 surfaceNormal_={{0,0,1}};

        // residue name 
        std::string resName_;

        // COM vector
        std::vector<Real3> COM_;

        // Angle with surface 
        std::vector<Real> AngleWithSurface_;

        // the uij for each LC molecule
        std::vector<Real3> uij_;

        // max of the r 
        Real rmax_;

        // max of the cosine
        Real maxcos_;
};