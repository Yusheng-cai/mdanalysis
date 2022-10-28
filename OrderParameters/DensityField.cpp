#include "DensityField.h"

DensityField::DensityField(const DensityFieldInput& input)
:simstate_(input.simstate), nL_(input.nL), n_(input.n), sigma_(input.sigma), isoSurfaceVal_(input.isoval), pbcmesh_(input.pbc)
{
    // calculate the actual cut off
    cutoff_  = n_*sigma_;
    sigmasq_ = sigma_ * sigma_;

    // resize the density field
    density_.resize(nL_,0);
    prefactor_ = std::pow(2.0 * Constants::PI * sigmasq_, -1.5);

    // initialize the buffer object
    for (auto it = density_buffer_.beginworker(); it != density_buffer_.endworker();it++){
        it -> resize(nL_,0);
    }
}

void DensityField::update(){
    L_  = simstate_.getSimulationBox().getSides();
    dL_ = L_ / nL_;

    CalcOffsetIndex();

    // zero out all the buffer obj
    for (auto it = density_buffer_.beginworker(); it != density_buffer_.endworker();it++){
        it -> fill(0);
    }
    density_.fill(0);
}

inline DensityField::Real DensityField::GaussianCoarseGrainFunction(Real distsq){
    return prefactor_*std::exp(-distsq/(2*sigmasq_));
}

void DensityField::CalcOffsetIndex(){
    offsetIndex_.clear();

    // get the differentials in the 3 directions
    Real dx = dL_[0];
    Real dy = dL_[1];
    Real dz = dL_[2]; 

    int Nx_offset = std::round(cutoff_/dx);
    int Ny_offset = std::round(cutoff_/dy);
    int Nz_offset = std::round(cutoff_/dz);

    for( int i=-Nx_offset; i<=Nx_offset;i++){
        for (int j=-Ny_offset; j<=Ny_offset; j++){
            for (int k=-Nz_offset; k<=Nz_offset;k++){
                INT3 id = {{i,j,k}};
                offsetIndex_.push_back(id);
            }
        }
    }
}

DensityField::INT3 DensityField::CalculateLatticePoint(const Real3& pos){
    INT3 index;
    for (int i=0;i<3;i++){
        index[i] = std::round(pos[i] / dL_[i]);
    }

    return index;
}

DensityField::Real3 DensityField::CalculateLatticePosition(const INT3& index){
    Real3 pos;
    for (int i=0;i<3;i++){
        pos[i] = index[i] * dL_[i];
    }

    return pos;
}

void DensityField::CalculateInstantaneousField(const std::vector<Real3>& pos, Mesh& mesh){
    // update
    update();

    if (pbcmesh_){
        mesh.setBoxLength(L_);
    }
    std::cout << "dL = " << dL_ << "\n";
    std::cout << "L = " << L_ << "\n";

    // start calculating InstantaneousInterface
    density_buffer_.set_master_object(density_);
    #pragma omp parallel
    {
        auto& localField = density_buffer_.access_buffer_by_id();

        #pragma omp for 
        for (int i=0;i<pos.size();i++){
            // find the lattice point 
            INT3 Index = CalculateLatticePoint(pos[i]);
            Index = CalculationTools::correctPBCLatticeIndex(Index, nL_);

            // iterate over all the indices and calculate instantaneousinterface
            for(int j=0;j<offsetIndex_.size();j++){
                // offset index 
                INT3 RealIndex = Index + offsetIndex_[j];
                RealIndex = CalculationTools::correctPBCLatticeIndex(RealIndex, nL_);

                // plate the positions on grid
                Real3 latticepos = CalculateLatticePosition(RealIndex);

                // calculate the distance 
                Real3 distance;
                Real distsq;
                simstate_.getSimulationBox().calculateDistance(latticepos, pos[i], distance, distsq);

                Real val = GaussianCoarseGrainFunction(distsq);
                localField(RealIndex) = localField(RealIndex) + val;
            }
        }    
    }

    #pragma omp parallel for
    for (int i=0;i<density_.getSize();i++){
        for (auto it = density_buffer_.beginworker();it != density_buffer_.endworker(); it++){
            density_[i] += it -> operator[](i);
        }
    }

    MarchingCubes_.triangulate_field(density_, mesh, dL_, nL_, isoSurfaceVal_, pbcmesh_);
}