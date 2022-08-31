#include "QtensorLattice.hpp"

namespace CalculationRegistry
{
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

    // initialize residue
    initializeResidueGroup(resname_);

    // resize the lattice 
    lattice_.resize(lattice_shape_);
    lattice_Qtensor_.resize(lattice_shape_, {});
    lattice_num_atoms_.resize(lattice_shape_,0.0);

    // cell grid
    cell_  = cellptr(new CellGrid(simstate_, cutoff_, 1));

    // register outputs
    registerOutputFunction("Director", [this](std::string name) -> void {this -> printDirector(name);});
    registerOutputFunction("Order", [this](std::string name) -> void {this -> printOrder(name);});
    registerOutputFunction("Qtensor",[this](std::string name) -> void {this -> printQtensor(name);});
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
    for (int i=0;i<lattice_.getSize();i++)
    {
        INT3 index3 = lattice_.getIndex3(i);
        Real3 pos;
        for (int j=0;j<3;j++)
        {
            pos[j] = index3[j] * dL_[j];
        }
        lattice_(index3) = pos;
    }

    // start the cell grid calculations 
    #pragma omp parallel for
    for (int i=0;i<lattice_.getSize();i++)
    {
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
        for (int j=0;j<neighborIdx.size();j++)
        {
            // the neighbot index
            int ind = neighborIdx[j];

            for (int k=0;k<cellIndices[ind].size();k++)
            {
                // the residue index 
                int resIndex = cellIndices[ind][k];

                // check the distance between the lattice point and the COM of the LC molecule  
                Real3 distance;
                Real distsq;
                simstate_.getSimulationBox().calculateDistance(latticePos, COM_[resIndex], distance, distsq);

                if (distsq < cutoff_sq_)
                {
                    Matrix localQ = LinAlg3x3::LocalQtensor(uij_[resIndex]);
                    Qtensor = Qtensor + localQ;
                    sum += 1;
                    vectorDistsq.push_back(distsq);
                }
            }
        }

        if (sum != 0)
        {
            Real min = *std::min_element(vectorDistsq.begin(), vectorDistsq.end());

            if (min < min_dist_sq_)
            {
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
    for (int i=0;i<3;i++)
    {
        dL_[i] = sides[i] / lattice_shape_[i];
    }
}

void QtensorLattice::finishCalculate()
{
    lattice_director_.resize(lattice_shape_);
    lattice_order_.resize(lattice_shape_);

    #pragma omp parallel for
    for (int i=0;i<lattice_Qtensor_.getSize();i++)
    {
        INT3 index3 = lattice_Qtensor_.getIndex3(i);
        if (lattice_num_atoms_(index3) != 0)
        {
            lattice_Qtensor_(index3) = lattice_Qtensor_(index3) * (0.5/lattice_num_atoms_(index3));

            auto res = LinAlg3x3::OrderEigenSolver(lattice_Qtensor_(index3));
            Real3 d;
            for (int j=0;j<3;j++)
            {
                d[j] = res.second[j][0];
            }
            lattice_order_(index3) = res.first[0];
            lattice_director_(index3) = d;
        }
        else
        {
            lattice_order_(index3) = -1.0;
            lattice_director_(index3) = {{0,0,0}};
        }
    }
}

void QtensorLattice::printDirector(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);
    Real3 zero_vec={{0,0,0}};

    for (int i=0;i<lattice_shape_[0];i++)
    {
        for (int j=0;j<lattice_shape_[1];j++)
        {
            for (int k=0;k<lattice_shape_[2];k++)
            {
                if (lattice_director_(i,j,k) != zero_vec)
                {
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

    for (int i=0;i<lattice_shape_[0];i++)
    {
        for (int j=0;j<lattice_shape_[1];j++)
        {
            for (int k=0;k<lattice_shape_[2];k++)
            {
                ofs << i << " " << j << " " << k << " " << lattice_order_(i,j,k) << "\n";
            }
        }
    }
    
    ofs.close();
}

void QtensorLattice::printQtensor(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);
    for (int i=0;i<lattice_shape_[0];i++)
    {
        for (int j=0;j<lattice_shape_[1];j++)
        {
            for (int k=0;k<lattice_shape_[2];k++)
            {
                ofs << i << " " << j << " " << k << " " << lattice_Qtensor_(i,j,k)[0][0] << " " <<  \
                lattice_Qtensor_(i,j,k)[0][1] << " " << lattice_Qtensor_(i,j,k)[0][2] << " " \
                << lattice_Qtensor_(i,j,k)[1][0] << " " << lattice_Qtensor_(i,j,k)[1][1] << "\n";
            }
        }
    }

    ofs.close();
}