#include "CellGrid.h"

CellGrid::CellGrid(SimulationState& simstate, Real dL, int searchnum)
: simstate_(simstate), dL_(dL), searchnum_(searchnum)
{
    for (int i=-searchnum_;i<=searchnum_;i++){
        for (int j=-searchnum_;j<=searchnum_;j++){
            for (int k=-searchnum_;k<=searchnum_;k++){
                index3 ind = {{i,j,k}};
                Offsets_.push_back(ind);
            }
        }
    }
}

void CellGrid::update(){
    Real3 Sides = simstate_.getSimulationBox().getSides();
    // Example like 5 / 3 --> 1 but N should be 2  -- > 0 1
    for (int i=0;i<3;i++){
        N_[i] = (int)std::ceil(Sides[i]/dL_);
    }
    totalIndices_ = N_[0] * N_[1] * N_[2];
}

std::vector<std::vector<int>> CellGrid::calculateIndices(const std::vector<Real3>& pos)
{
    std::vector<std::vector<int>> ret_indices(totalIndices_);

    #pragma omp parallel
    {
        std::vector<std::vector<int>> ret_indices_local(totalIndices_);
        #pragma omp for
        for (int i=0;i<pos.size();i++){
            int index = getCellGridIntIndex(pos[i]);
            ret_indices_local[index].push_back(i);
        }

        #pragma omp critical
        {
            for (int i=0;i<totalIndices_;i++){
                ret_indices[i].insert(ret_indices[i].end(), ret_indices_local[i].begin(), ret_indices_local[i].end());
            }
        }
    }

    return ret_indices;
}

int CellGrid::ConvertGridIndexToIndex(index3& index){
    // convert index into x + y * Nx + z * Nx * Ny
    int sum = index[0] + index[1] * N_[0] + index[2] * N_[1] * N_[0];
    ASSERT((sum >= 0) && (sum < totalIndices_) , "The Cell Grid Index " << index << " is out of bounds.");

    return sum;
}

CellGrid::index3 CellGrid::getCellGridIndex(const Real3& pos)
{
    // first let's shift the position into the simulation box
    Real3 shiftedPos = simstate_.getSimulationBox().shiftIntoBox(pos);

    // Then let's put the shifted Pos into a cell grid
    index3 index;

    for (int i=0;i<3;i++){
        // take care of edge cases when shiftedPos is exactly N[i] * dL
        if (std::abs(shiftedPos[i] - N_[i] * dL_) < 1e-5) {
            index[i] = N_[i] - 1;
        }
        else if (std::abs(shiftedPos[i] - 0) < 1e-5){
            index[i] = 0;
        }
        else{
            index[i] = (int)std::floor(shiftedPos[i] / dL_ );
        }
        ASSERT((index[i] < N_[i]) && (index[i] >= 0), "Index out of range, index is " << index[i] << " max is " << N_[i] << " Position = " << shiftedPos[i] << " dL = " << dL_);
    }

    return index;
}

int CellGrid::getCellGridIntIndex(const Real3& pos)
{
    index3 gridIndex = getCellGridIndex(pos);
    
    return ConvertGridIndexToIndex(gridIndex);
}

CellGrid::index3 CellGrid::FixIndex(index3& index)
{
    index3 ret;

    for (int i=0;i<3;i++){
        if (index[i] < 0){
            index[i] += N_[i];
        }

        ret[i] = index[i] % N_[i];
    }

    return ret;
}

std::vector<int> CellGrid::getNeighborIndex(const Real3& pos)
{
    index3 index = getCellGridIndex(pos);
    std::vector<int> neighborIndex;

    for (auto& off : Offsets_){
        index3 newIndex = index + off;
        newIndex = FixIndex(newIndex);
        int ind  = ConvertGridIndexToIndex(newIndex);
        neighborIndex.push_back(ind);
    }

    return neighborIndex;
}