#pragma once

#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "SimulationState.h"
#include "tools/CommonOperations.h"

#include <vector>
#include <array>
#include <iostream>
#include <string>
#include <cmath>

class CellGrid{
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;
        using index3= CommonTypes::index3;
        using INT2 = CommonTypes::index2;

        // initialize cell grid 
        CellGrid(SimulationState& simstate, Real dL);

        // update the cell grid dimensions every iteration
        void update();

        // no thread-safe
        std::vector<std::vector<int>> calculateIndices(const std::vector<Real3>& pos);

        index3 getCellGridIndex(const Real3& pos);
        int ConvertGridIndexToIndex(index3& index);
        int getCellGridIntIndex(const Real3& pos);
        std::vector<int> getNeighborIndex(const Real3& pos);
        index3 FixIndex(index3& index);
        int getSize() const {return N_[0] * N_[1] * N_[2];}

    private:
        Real dL_, L_;
        SimulationState& simstate_;
        index3 N_;
        int totalIndices_;

        std::vector<index3> Offsets_;
        std::vector<INT2> searchnum_;
};