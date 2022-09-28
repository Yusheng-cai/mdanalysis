#include "QtensorLattice.hpp"

namespace CalculationRegistry{
    registry_<QtensorLattice> QtensorLatticeRegister_("QtensorLattice");
}

QtensorLattice::QtensorLattice(const CalculationInput& input)
: Calculation(input)
{
    pack_.ReadArrayNumber("LatticeShape", ParameterPack::KeyType::Required, lattice_shape_);
    pack_.ReadString("residue", ParameterPack::KeyType::Required, resname_);

    // read in the cut off (how far do we search near a lattice point for COM)
    pack_.ReadNumber("cutoff", ParameterPack::KeyType::Required, cutoff_);
    cutoff_sq_ = cutoff_ * cutoff_;

    // read in head and tail index
    pack_.ReadNumber("headindex", ParameterPack::KeyType::Optional, headIndex_);
    pack_.ReadNumber("tailindex", ParameterPack::KeyType::Optional, tailIndex_);
    headIndex_--;
    tailIndex_--;

    // the minimum distance in which there needs to be a COM point near the lattice point
    pack_.ReadNumber("minDist", ParameterPack::KeyType::Optional, min_dist_);
    min_dist_sq_ = min_dist_ * min_dist_;

    // whether we are reducing the dimensions 
    reduced_ = pack_.ReadArrayNumber("reduced_dimension", ParameterPack::KeyType::Optional, reduced_dimensions_);

    // initialize residue
    initializeResidueGroup(resname_);

    // resize the lattice 
    lattice_.resize(lattice_shape_);
    lattice_Qtensor_.resize(lattice_shape_, {});
    lattice_num_atoms_.resize(lattice_shape_,0.0);
    lattice_biaxiality_.resize(lattice_shape_,{});

    // cell grid
    cell_  = cellptr(new CellGrid(simstate_, cutoff_, 1));

    // register outputs
    registerOutputFunction("Director", [this](std::string name) -> void {this -> printDirector(name);});
    registerOutputFunction("Order", [this](std::string name) -> void {this -> printOrder(name);});
    registerOutputFunction("Qtensor",[this](std::string name) -> void {this -> printQtensor(name);});
    registerOutputFunction("velocity", [this](std::string name) -> void {this -> printVelocity(name);});
    registerOutputFunction("reducedOrder", [this](std::string name) -> void {this -> printReducedOrder(name);});
    registerOutputFunction("reducedDirector", [this](std::string name) -> void {this -> printReducedDirector(name);});
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
    for (int i=0;i<COM_.size();i++)
    {
        COM_[i] = calcCOM(res[i]);

        Real3 distance;
        Real distsq;

        Real3 headPos = res[i].atoms_[headIndex_].positions_;
        Real3 tailPos = res[i].atoms_[tailIndex_].positions_;

        simstate_.getSimulationBox().calculateDistance(headPos, tailPos, distance, distsq);

        LinAlg3x3::normalize(distance);

        uij_[i] = distance;
    }

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
    #pragma omp parallel for
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

                if (distsq < cutoff_sq_){
                    Matrix localQ = LinAlg3x3::LocalQtensor(uij_[resIndex]);
                    Qtensor = Qtensor + localQ;
                    sum += 1;
                    vectorDistsq.push_back(distsq);
                }
            }
        }

        if (sum != 0){
            Real min = *std::min_element(vectorDistsq.begin(), vectorDistsq.end());

            if (min < min_dist_sq_){
                lattice_Qtensor_(index3) = lattice_Qtensor_(index3) + Qtensor;
                lattice_num_atoms_(index3) += sum;
            }
        }
    }
}

void QtensorLattice::update()
{
    // update the cell grid
    cell_->update();

    // update the sides 
    Real3 sides = simstate_.getSimulationBox().getSides();
    for (int i=0;i<3;i++){
        dL_[i] = sides[i] / lattice_shape_[i];
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
            lattice_order_(index3) = -1.0;
            lattice_director_(index3) = {{0,0,0}};
            lattice_biaxiality_(index3) = -0.0;
        }
    }

    if (reduced_){
        INT3 shape=lattice_Qtensor_.getShape();
        int dim1,dim2;
        dim1 = reduced_dimensions_[0];
        dim2 = reduced_dimensions_[1];
        int size1 = shape[dim1];
        int size2 = shape[dim2];
        reduced_Q_.resize(size1, std::vector<Matrix>(size2,{{}}));
        reduced_num_.resize(size1, std::vector<Real>(size2,0.0));
        reduced_order_.resize(size1, std::vector<Real>(size2));
        reduced_director_.resize(size1, std::vector<Real3>(size2));

        // first let's reduce one of the dimensions 
        #pragma omp parallel
        {
            std::vector<std::vector<Matrix>> localQtensor(size1, std::vector<Matrix>(size2,{{}}));
            std::vector<std::vector<Real>> localNum(size1, std::vector<Real>(size2,0.0));
            #pragma omp for
            for (int i=0;i<lattice_Qtensor_.getSize();i++){
                INT3 index3 = lattice_Qtensor_.getIndex3(i);
                if (lattice_num_atoms_(index3) != 0){
                    localNum[index3[dim1]][index3[dim2]] += lattice_num_atoms_(index3);
                    localQtensor[index3[dim1]][index3[dim2]] = localQtensor[index3[dim1]][index3[dim2]] + lattice_Qtensor_(index3);
                }
            }
            #pragma omp critical
            {
                reduced_Q_ = reduced_Q_ + localQtensor;
                reduced_num_ = reduced_num_ + localNum;
            }
        }

        #pragma omp parallel
        {
            #pragma omp parallel for
            for (int i=0;i<size1;i++){
                for (int j=0;j<size2;j++){
                    if (reduced_num_[i][j] != 0){
                        Matrix thisQ = reduced_Q_[i][j] * (0.5/reduced_num_[i][j]);
                        auto res = LinAlg3x3::OrderEigenSolver(thisQ);
                        reduced_order_[i][j] = res.first[0];
                        Real3 dir;
                        for (int k=0;k<3;k++){
                            dir[k] = res.second[k][0];
                        }
                        reduced_director_[i][j] = dir;
                    }
                }
            }

        }
    }
    
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