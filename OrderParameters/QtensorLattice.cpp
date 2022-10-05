#include "QtensorLattice.hpp"

namespace CalculationRegistry{
    registry_<QtensorLattice> QtensorLatticeRegister_("QtensorLattice");
}

QtensorLattice::QtensorLattice(const CalculationInput& input)
: Calculation(input)
{
    pack_.ReadArrayNumber("LatticeShape", ParameterPack::KeyType::Required, lattice_shape_);
    pack_.ReadString("residue", ParameterPack::KeyType::Required, resname_);

    // read in head and tail index
    pack_.ReadNumber("headindex", ParameterPack::KeyType::Optional, headIndex_);
    pack_.ReadNumber("tailindex", ParameterPack::KeyType::Optional, tailIndex_);
    headIndex_--;
    tailIndex_--;

    pack_.Readbool("coarse_grain", ParameterPack::KeyType::Optional, coarse_grain_);

    if (coarse_grain_){
        pack_.ReadNumber("sigma", ParameterPack::KeyType::Required, sigma_);
        pack_.ReadNumber("n", ParameterPack::KeyType::Required, n_);
        sigma2_ = sigma_ * sigma_;
        prefactor_ = std::pow(2*Constants::PI*sigma2_, -1.5); 
        inv_factor_ = 1.0 / (2.0 * sigma2_);
    }
    else{
        pack_.ReadNumber("cutoff", ParameterPack::KeyType::Required, cutoff_);
        cutoff_sq_ = cutoff_ * cutoff_;
        pack_.ReadNumber("minDist", ParameterPack::KeyType::Required, min_dist_);
        min_dist_sq_ = min_dist_ * min_dist_;
        cell_ = cellptr(new CellGrid(simstate_, cutoff_, 1));
    }

    pack_.Readbool("performMC", ParameterPack::KeyType::Optional, performMC_);
    if (performMC_){
        pack_.ReadNumber("isoval", ParameterPack::KeyType::Required, isoval_);
        pack_.Readbool("pbc", ParameterPack::KeyType::Optional, pbc_);
    }

    // whether we are reducing the dimensions 
    reduced_ = pack_.ReadArrayNumber("reduced_dimension", ParameterPack::KeyType::Optional, reduced_dimensions_);

    // initialize residue
    initializeResidueGroup(resname_);

    // resize the lattice 
    lattice_.resize(lattice_shape_);
    lattice_Qtensor_.resize(lattice_shape_, {});
    lattice_num_atoms_.resize(lattice_shape_,0.0);
    lattice_biaxiality_.resize(lattice_shape_,{});

    // register outputs
    registerOutputFunction("Director", [this](std::string name) -> void {this -> printDirector(name);});
    registerOutputFunction("Order", [this](std::string name) -> void {this -> printOrder(name);});
    registerOutputFunction("Qtensor",[this](std::string name) -> void {this -> printQtensor(name);});
    registerOutputFunction("velocity", [this](std::string name) -> void {this -> printVelocity(name);});
    registerOutputFunction("reducedOrder", [this](std::string name) -> void {this -> printReducedOrder(name);});
    registerOutputFunction("reducedDirector", [this](std::string name) -> void {this -> printReducedDirector(name);});
    registerOutputFunction("ply", [this](std::string name) -> void {this -> printIsoSurface(name);});

    // per iter outputs
    registerPerIterOutputFunction("Order", [this](std::ofstream& ofs) -> void {this -> printOrderPerIter(ofs);});
    registerPerIterOutputFunction("density", [this](std::ofstream& ofs) -> void {this -> printDensityPerIter(ofs);});
}

void QtensorLattice::calculateCoraseGrain(){
    #pragma omp parallel
    {   
        Lattice<Matrix> mat_local(lattice_shape_, {});
        Lattice<Real> num_atom_local(lattice_shape_,0);

        #pragma omp for
        for (int i=0;i<COM_.size();i++){
            INT3 index3 = CalculationTools::NearestLatticeIndex(COM_[i], dL_, lattice_shape_);
            Matrix Q = LinAlg3x3::LocalQtensor(uij_[i]);
            for (auto off : lattice_offsets_){
                // offset index and position
                INT3 offsetIndex = CalculationTools::correctPBCLatticeIndex(index3 + off, lattice_shape_);
                Real3 pos = dL_ * offsetIndex;

                // find diff between atom and lattice
                Real3 vecdist;
                Real distsq;
                simstate_.getSimulationBox().calculateDistance(COM_[i], pos, vecdist, distsq);

                // find coarse graining factor
                Real coarse_grain_factor = CalculateCoraseGrainFunction(distsq);
                mat_local(offsetIndex) = mat_local(offsetIndex) + Q * coarse_grain_factor;
                num_atom_local(offsetIndex) = num_atom_local(offsetIndex) + coarse_grain_factor;
            }
        }

        #pragma omp critical
        {
            // total
            lattice_Qtensor_ += mat_local;
            lattice_num_atoms_ += num_atom_local;

            // per iter
            lattice_Qtensor_Iter_ += mat_local;
            lattice_num_atoms_Iter_ += num_atom_local;
        }
    }

    #pragma omp parallel for
    for (int i=0;i<lattice_Qtensor_Iter_.getSize();i++){
        INT3 index3 = lattice_Qtensor_Iter_.getIndex3(i);
        if (lattice_num_atoms_Iter_(index3) != 0){
            lattice_Qtensor_Iter_(index3) = lattice_Qtensor_Iter_(index3) * (0.5 / lattice_num_atoms_Iter_(index3));
        }
    }
}

void QtensorLattice::calculateNeighborSearch(){
    // calculate the cell indices 
    std::vector<std::vector<int>> cellIndices = cell_->calculateIndices(COM_);

    // fill the lattice positions 
    #pragma omp parallel for 
    for (int i=0;i<lattice_.getSize();i++){
        INT3 index3 = lattice_.getIndex3(i);
        Real3 pos;
        for (int j=0;j<3;j++){
            pos[j] = index3[j] * dL_[j];
        }
        lattice_(index3) = pos;
    }

    // start the cell grid calculations 
    #pragma omp parallel
    {
        #pragma omp for
        for (int i=0;i<lattice_.getSize();i++){
            INT3 index3 = lattice_.getIndex3(i);
            Real3 latticePos = lattice_(index3);
            int selfIndex = cell_->getCellGridIntIndex(latticePos);
            std::vector<int> neighborIdx = cell_->getNeighborIndex(latticePos);
            neighborIdx.push_back(selfIndex);

            // find the COMs in the neighborhood
            Matrix Qtensor = {};

            // sum of the number of atoms in the neighborhood
            int sum=0;
            // vector that keeps track of all the distances 
            std::vector<Real> vectorDistsq;
            for (int j=0;j<neighborIdx.size();j++){
                // the neighbot index
                int ind = neighborIdx[j];
                for (int k=0;k<cellIndices[ind].size();k++){
                    // the residue index 
                    int resIndex = cellIndices[ind][k];

                    // check the distance between the lattice point and the COM of the LC molecule  
                    Real3 distance;
                    Real distsq;
                    simstate_.getSimulationBox().calculateDistance(latticePos, COM_[resIndex], distance, distsq);

                    // the distance must be less than the cutoff distance
                    if (distsq <= cutoff_sq_){
                        Matrix localQ = LinAlg3x3::LocalQtensor(uij_[resIndex]);
                        Qtensor = Qtensor + localQ;
                        sum += 1;
                        vectorDistsq.push_back(distsq);
                    }
                }
            }

            // update the local per iter
            lattice_Qtensor_Iter_(index3) = Qtensor;
            lattice_num_atoms_Iter_(index3) = sum;

            if (sum != 0){
                Real min = *std::min_element(vectorDistsq.begin(), vectorDistsq.end());

                if (min < min_dist_sq_){
                    lattice_Qtensor_(index3) = lattice_Qtensor_(index3) + Qtensor;
                    lattice_num_atoms_(index3) += sum;
                }
            }
        }
    }

    #pragma omp parallel for
    for (int i=0;i<lattice_Qtensor_Iter_.getSize();i++){
        INT3 index3 = lattice_Qtensor_Iter_.getIndex3(i);
        if (lattice_num_atoms_Iter_(index3) != 0){
            lattice_Qtensor_Iter_(index3) = lattice_Qtensor_Iter_(index3) * (0.5 / lattice_num_atoms_Iter_(index3));
        }
    }
}

void QtensorLattice::calculate()
{
    // assign the cell grid 
    const auto& res = getResidueGroup(resname_).getResidues();

    // resize uij and COM 
    uij_.resize(res.size());
    COM_.resize(res.size());

    // calculate the center of mass as well as the uij
    #pragma omp parallel for
    for (int i=0;i<COM_.size();i++){
        COM_[i] = calcCOM(res[i]);

        Real3 distance;
        Real distsq;

        Real3 headPos = res[i].atoms_[headIndex_].positions_;
        Real3 tailPos = res[i].atoms_[tailIndex_].positions_;

        simstate_.getSimulationBox().calculateDistance(headPos, tailPos, distance, distsq);

        LinAlg3x3::normalize(distance);

        uij_[i] = distance;
    }

    if (coarse_grain_){
        calculateCoraseGrain();
    }
    else{
        calculateNeighborSearch();
    }

}

void QtensorLattice::update(){
    // update the sides 
    Real3 sides = simstate_.getSimulationBox().getSides();
    for (int i=0;i<3;i++){
        dL_[i] = sides[i] / lattice_shape_[i];
    }

    // only calculate lattice offsets if we are doing coarse-grained simulation
    if (coarse_grain_){
        // update the cell grid
        cell_->update();
        CalculateLatticeOffsets();
    }

    lattice_num_atoms_Iter_.resize(lattice_shape_,0);
    lattice_Qtensor_Iter_.resize(lattice_shape_, {});
}

void QtensorLattice::CalculateLatticeOffsets(){
    INT3 num_offsets;
    lattice_offsets_.clear();

    for (int i=0;i<3;i++){
        num_offsets[i] = std::round(n_ * sigma_ / dL_[i]);
    }

    for (int i=-num_offsets[0];i<num_offsets[0];i++){
        for (int j=-num_offsets[1];j<num_offsets[1];j++){
            for (int k=-num_offsets[2];k<num_offsets[2];k++){
                INT3 index = {{i,j,k}};
                lattice_offsets_.push_back(index);
            }
        }
    }
}

void QtensorLattice::finishCalculate(){
    lattice_director_.resize(lattice_shape_);
    lattice_order_.resize(lattice_shape_);

    #pragma omp parallel for
    for (int i=0;i<lattice_Qtensor_.getSize();i++){
        INT3 index3 = lattice_Qtensor_.getIndex3(i);
        if (lattice_num_atoms_(index3) != 0){
            lattice_Qtensor_(index3) = lattice_Qtensor_(index3) * (0.5/lattice_num_atoms_(index3));

            auto res = LinAlg3x3::OrderEigenSolver(lattice_Qtensor_(index3));
            Real biaxi = res.first[1] * 2 + res.first[0];
            Real3 d;
            for (int j=0;j<3;j++){
                d[j] = res.second[j][0];
            }
            lattice_order_(index3) = res.first[0];
            lattice_director_(index3) = d;
            lattice_biaxiality_(index3) = biaxi;
        }
        else{
            lattice_order_(index3) = 100.0;
            lattice_director_(index3) = {{0,0,0}};
            lattice_biaxiality_(index3) = 100.0;
}
    }

    // calculate MC
    if (performMC_){
        mc_.triangulate_field(lattice_order_, m_, dL_, lattice_shape_, isoval_, pbc_);
    }
}

void QtensorLattice::printIsoSurface(std::string name){
    ASSERT((performMC_), "Marching cubes options is not turned on.");
    MeshTools::writePLY(name, m_);
}

void QtensorLattice::printDirector(std::string name){
    std::ofstream ofs;
    ofs.open(name);
    Real3 zero_vec={{0,0,0}};

    for (int i=0;i<lattice_shape_[0];i++){
        for (int j=0;j<lattice_shape_[1];j++){
            for (int k=0;k<lattice_shape_[2];k++){
                if (lattice_director_(i,j,k) != zero_vec){
                    ofs << i << " " << j << " " << k << " " << lattice_director_(i,j,k)[0] << " " << lattice_director_(i,j,k)[1] << " " << lattice_director_(i,j,k)[2] << "\n";
                }
            }
        }
    }

    ofs.close();
}

void QtensorLattice::printOrder(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    for (int i=0;i<lattice_shape_[0];i++){
        for (int j=0;j<lattice_shape_[1];j++){
            for (int k=0;k<lattice_shape_[2];k++){
                ofs << i << " " << j << " " << k << " " << lattice_order_(i,j,k) << " " << \
                lattice_biaxiality_(i,j,k) << "\n";
            }
        }
    }
    
    ofs.close();
}

void QtensorLattice::printQtensor(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);
    for (int i=0;i<lattice_shape_[0];i++){
        for (int j=0;j<lattice_shape_[1];j++){
            for (int k=0;k<lattice_shape_[2];k++){
                ofs << i << " " << j << " " << k << " " << lattice_Qtensor_(i,j,k)[0][0] << " " <<  \
                lattice_Qtensor_(i,j,k)[0][1] << " " << lattice_Qtensor_(i,j,k)[0][2] << " " \
                << lattice_Qtensor_(i,j,k)[1][0] << " " << lattice_Qtensor_(i,j,k)[1][1] << "\n";
            }
        }
    }

    ofs.close();
}

void QtensorLattice::printVelocity(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);
    Real3 zero_vec = {};
    for (int i=0;i<lattice_shape_[0];i++){
        for (int j=0;j<lattice_shape_[1];j++){
            for (int k=0;k<lattice_shape_[2];k++){
                if (lattice_director_(i,j,k) != zero_vec){
                    ofs << i << " " << j << " " << k << " " << lattice_director_(i,j,k)[0] * lattice_order_(i,j,k) \
                     << " " << lattice_director_(i,j,k)[1] * lattice_order_(i,j,k) \
                     << " " << lattice_director_(i,j,k)[2] * lattice_order_(i,j,k) << "\n";
                }
            }
        }
    }

    ofs.close();
}

void QtensorLattice::printReducedDirector(std::string name){
    ASSERT((reduced_), "Did not specify reduced dimensions in the input file!");
    std::ofstream ofs;
    ofs.open(name);
    for (int i=0;i<reduced_director_.size();i++){
        for (int j=0;j<reduced_director_[i].size();j++){
            ofs << i << " " << j << " " << reduced_director_[i][j] << "\n";
        }
    }

    ofs.close();
}

void QtensorLattice::printReducedOrder(std::string name){
    ASSERT((reduced_), "Did not specify reduced dimensions in the input file!");
    std::ofstream ofs;
    ofs.open(name);
    for (int i=0;i<reduced_order_.size();i++){
        for (int j=0;j<reduced_order_[i].size();j++){
            ofs << i << " " << j << " " << reduced_order_[i][j] << "\n";
        }
    }

    ofs.close();
}

void QtensorLattice::printOrderPerIter(std::ofstream& ofs){
    Lattice<Real> order(lattice_shape_,0.0);

    #pragma omp parallel for
    for (int i=0;i<lattice_Qtensor_Iter_.getSize();i++){
        INT3 index3 = lattice_Qtensor_Iter_.getIndex3(i);
        if (lattice_num_atoms_Iter_(index3) != 0){
            auto res = LinAlg3x3::OrderEigenSolver(lattice_Qtensor_Iter_(index3));
            order(index3) = res.first[0];
        }
        else{
            order(index3) = 100;
        }
    }

    for (int i=0;i<lattice_Qtensor_Iter_.getSize();i++){
        INT3 index3 = lattice_Qtensor_Iter_.getIndex3(i);
        ofs << order(index3) << " ";
    }

    ofs << "\n";
}

void QtensorLattice::printDensityPerIter(std::ofstream& ofs){
    for (int i=0;i<lattice_num_atoms_Iter_.getSize();i++){
        ofs << lattice_num_atoms_Iter_[i] << " ";
    }

    ofs << "\n";
}

QtensorLattice::Real QtensorLattice::CalculateCoraseGrainFunction(Real& rsq){
    return prefactor_ * std::exp(-rsq * inv_factor_);
}