#include "QtensorVoxel.hpp"

QtensorVoxel::QtensorVoxel(const CalculationInput& input)
: Calculation(input)
{
    // first read the number of voxels
    pack_.ReadArrayNumber("VoxelSize", ParameterPack::KeyType::Optional, voxel_size_);

    // then read the residuegroup
    pack_.ReadString("residue", ParameterPack::KeyType::Required, residue_group_);
    initializeResidueGroup(residue_group_);

    // read head and tail index
    pack_.ReadNumber("headindex", ParameterPack::KeyType::Optional, head_index_);
    pack_.ReadNumber("tailindex", ParameterPack::KeyType::Optional, tail_index_);
    head_index_--;
    tail_index_--;

    // reshape qtensor voxels
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            zero_matrix_[i][j] = 0.0;
        }
    }
    Qtensor_Voxels_.resize(voxel_size_, zero_matrix_);
    NumAtoms_Voxels_.resize(voxel_size_, 0.0);

    Qtensor_voxels_buffer_.set_master_object(Qtensor_Voxels_);
    for (auto it = Qtensor_voxels_buffer_.beginworker(); it != Qtensor_voxels_buffer_.endworker(); it++){
        it -> resize(voxel_size_, zero_matrix_);
    }
    NumAtoms_voxels_buffer_.set_master_object(NumAtoms_Voxels_);
    for (auto it = NumAtoms_voxels_buffer_.beginworker(); it != NumAtoms_voxels_buffer_.endworker(); it++){
        it -> resize(voxel_size_, 0.0);
    }
}

QtensorVoxel::INT3 QtensorVoxel::findVoxel(Real3& position){
    Real3 box = simstate_.getSimulationBox().getSides();
    INT3 index;

    for (int i=0;i<3;i++){
        index[i] = std::floor(position[i] / voxel_separation_[i]); 
    }

    return index;
}

void QtensorVoxel::reset_buffers(){
    Qtensor_voxels_buffer_.set_master_object(Qtensor_Voxels_);
    NumAtoms_voxels_buffer_.set_master_object(NumAtoms_Voxels_);

    for (auto it = Qtensor_voxels_buffer_.beginworker(); it != Qtensor_voxels_buffer_.endworker(); it++){
        it -> fill(zero_matrix_);
    }

    for (auto it = NumAtoms_voxels_buffer_.beginworker(); it != NumAtoms_voxels_buffer_.endworker(); it++){
        it -> fill(0);
    }

}

void QtensorVoxel::calculate()
{
    // calculate the COM
    auto& res = getResidueGroup(residue_group_).getResidues();

    // fill the buffers
    reset_buffers();

    // assign the cell grid 
    uij_.resize(res.size());
    COM_.resize(res.size());

    // initialize global qtensor
    Global_Qtensor_.fill({});

    // calculate the center of mass as well as the uij
    #pragma omp parallel
    {
        Matrix Q={};
        // local voxel lattice
        auto& LocalVoxel = Qtensor_voxels_buffer_.access_buffer_by_id();
        auto& LocalAtoms = NumAtoms_voxels_buffer_.access_buffer_by_id();

        #pragma omp for
        for (int i=0;i<COM_.size();i++){
            // calculate the center of mass 
            COM_[i] = calcCOM(res[i]);

            // see which voxel this COM fall
            INT3 voxel_num = findVoxel(COM_[i]);

            Real3 distance;
            Real distsq;
            Real3 headPos = res[i].atoms_[head_index_].positions_;
            Real3 tailPos = res[i].atoms_[tail_index_].positions_;
            simstate_.getSimulationBox().calculateDistance(headPos, tailPos, distance, distsq);
            LinAlg3x3::normalize(distance);

            uij_[i] = distance;

            // calculate the local qtensor = 3/2 ui uj - 1/2 delta ij 
            Matrix localQtensor = LinAlg3x3::LocalQtensor(uij_[i]);
            LocalVoxel(voxel_num) = LocalVoxel(voxel_num) + localQtensor;
            LocalAtoms(voxel_num) += 1;

            // accumate the qtensor
            Q = Q + localQtensor;
        }

        #pragma omp critical
        {
            Global_Qtensor_ = Global_Qtensor_ + Q;
        }
    }

    #pragma omp parallel for
    for (int i=0;i<Qtensor_Voxels_.getSize();i++){
        for (auto it = Qtensor_voxels_buffer_.beginworker(); it != Qtensor_voxels_buffer_.endworker(); it++){
            Qtensor_Voxels_[i] = Qtensor_Voxels_[i] + it->operator[](i);
        }

        for (auto it=NumAtoms_voxels_buffer_.beginworker(); it != NumAtoms_voxels_buffer_.endworker(); it++){
            NumAtoms_Voxels_[i] = NumAtoms_Voxels_[i] + it -> operator[](i);
        }
    }

    Global_Qtensor_ = Global_Qtensor_ / COM_.size();
}

void QtensorVoxel::finishCalculate(){
    Order_Voxels_.resize(voxel_size_, 0.0);
    Director_Voxels_.resize(voxel_size_);

    // first let's divide the qtensor by the number of atoms 
    #pragma omp parallel for
    for (int i=0;i<Qtensor_Voxels_.getSize();i++){
        if (NumAtoms_Voxels_[i] > 0){
            Qtensor_Voxels_[i] = Qtensor_Voxels_[i] / NumAtoms_Voxels_[i];
            auto result = LinAlg3x3::OrderEigenSolver(Qtensor_Voxels_[i]);
            Order_Voxels_[i] = result.first[0];
            Director_Voxels_[i] = result.second[0];
        }
        else{
            Order_Voxels_[i] = 100.0;
            Director_Voxels_[i] = {{0.0,0.0,0.0}};
        }
    }
}

void QtensorVoxel::update(){
    auto box_size = simstate_.getSimulationBox().getSides();
    voxel_separation_ = box_size / voxel_size_;
}