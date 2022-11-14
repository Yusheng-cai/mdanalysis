#include "CellGrid.h"

CellGrid::CellGrid(SimulationState& simstate, Real dL)
: simstate_(simstate), dL_(dL)
{
    searchnum_.resize(3);
}

void CellGrid::update(){
    Real3 Sides = simstate_.getSimulationBox().getSides();
    // Example like 7 / 3 --> 2 --> last one will be merged with the first 
    for (int i=0;i<3;i++){
        N_[i] = (int)std::floor(Sides[i]/dL_);
        if (N_[i] > 2){
            searchnum_[i] = {{-1,1}};
        }
        else{
            searchnum_[i] = {{0,1}};
        }
    }
    totalIndices_ = N_[0] * N_[1] * N_[2];

    // Now update updates the offsets as well 
    Offsets_.clear();
    for (int i=searchnum_[0][0];i<=searchnum_[0][1];i++){
        for (int j=searchnum_[1][0];j<=searchnum_[1][1];j++){
            for (int k=searchnum_[2][0];k<=searchnum_[2][1];k++){
                Offsets_.push_back({{i,j,k}});
            }
        }
    }
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
        if (std::abs(shiftedPos[i] - 0) < 1e-5){
            index[i] = 0;
        }
        else{
            index[i] = (int)std::floor(shiftedPos[i] / dL_ ) % N_[i];
        }
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