#include "SRE.h"

namespace CalculationRegistry
{
    registry_<SRE> registerSRE("SRE");
}

SRE::SRE(const CalculationInput& input) 
:Calculation(input)
{
    // This could be LC or water
    pack_.ReadString("solvent", ParameterPack::KeyType::Required, SolventName_);
    // This is usually SAM 
    pack_.ReadString("solute", ParameterPack::KeyType::Required, SoluteName_);
    // Read in the alpha parameter
    pack_.ReadNumber("epsilon", ParameterPack::KeyType::Required, epsilon_);
    // read in the cutoff
    pack_.ReadNumber("cutoff", ParameterPack::KeyType::Optional, cutoff_);
    alpha_ = 1.0 / epsilon_;
    // read in the mode (Brute/NS)
    pack_.ReadString("mode", ParameterPack::KeyType::Optional, mode_);
    ASSERT((mode_ == "NS" || mode_ == "Brute"), "The mode can only be Brute or NS.");

    cutoffsq_ = cutoff_ * cutoff_;

    addResidueGroup(SolventName_);
    addResidueGroup(SoluteName_);

    // register the printing function for energy 
    registerPerIterOutputFunction("energy", [this](std::ofstream& ofs) -> void {this -> printEnergyPerIter(ofs);});
    registerOutputFileOutputs("energy", [this](void) -> Real {return this -> getEnergy();});

    cell_ = Cellptr(new CellGrid(simstate_, cutoff_,1));

    getNonZeroCharges();
}

void SRE::getNonZeroCharges()
{
    const auto& solvent = getResidueGroup(SolventName_).getTotalResidue().atoms_;
    const auto& solute  = getResidueGroup(SoluteName_).getTotalResidue().atoms_;

    for (int i=0;i<solvent.size();i++)
    {
        Real charge = solvent[i].charge_;

        if (charge != 0)
        {
            NonZeroSolvent_.push_back(i);
        }
    }
    
    for (int i=0;i<solute.size();i++)
    {
        Real charge = solute[i].charge_;

        if (charge != 0)
        {
            NonZeroSolute_.push_back(i);
        }
    }
}

void SRE::update()
{
    cell_ -> update();
}

void SRE::calculateWithNS()
{
    const auto& solvent = getResidueGroup(SolventName_).getTotalResidue().atoms_;
    const auto& solute  = getResidueGroup(SoluteName_).getTotalResidue().atoms_;

    // make all solute into cellgrid
    int num = cell_->getSize();
    std::vector<std::vector<int>> celllist_(num);
    for (int i=0;i<NonZeroSolute_.size();i++)
    {
        int index = NonZeroSolute_[i];
        Real3 pos = solute[index].positions_;
        index3 id = cell_ -> getCellGridIndex(pos);

        int ind  = cell_ -> getCellGridIntIndex(pos);
        ASSERT((ind >= 0 && ind < num), "Index out of range, max is " << num << " while it is " << ind);

        celllist_[ind].push_back(index);
    }


    // iterate over all the atoms in solvent and solute
    #pragma omp parallel
    {
        Real sum = 0.0;
        #pragma omp for
        for (int i=0;i<NonZeroSolvent_.size();i++)
        {
            int solventInd = NonZeroSolvent_[i];
            Real3 possv = solvent[solventInd].positions_;
            Real qi     = solvent[solventInd].charge_;
            std::vector<int> Neighbors = cell_ -> getNeighborIndex(possv);

            for (int ind : Neighbors)
            {
                ASSERT((ind < num), "Ind is " << ind << " max is " << num);
                for (int index : celllist_[ind])
                {
                    Real3 possl = solute[index].positions_;
                    Real qj     = solute[index].charge_;
                    Real3 dist;
                    Real distsq;

                    simstate_.getSimulationBox().calculateDistance(possv, possl, dist, distsq);

                    if (distsq < cutoffsq_)
                    {
                        Real r = std::sqrt(distsq);

                        sum += factor_ * qi * qj * std::erfc(r * alpha_) / r;
                    }
                }
            }
        }

        #pragma omp critical
        {
            energy_ += sum;
        }
    }
}

void SRE::calculateWithoutNS()
{
    const auto& solvent = getResidueGroup(SolventName_).getTotalResidue().atoms_;
    const auto& solute  = getResidueGroup(SoluteName_).getTotalResidue().atoms_;

    #pragma omp parallel
    {
        Real sum = 0.0;
        #pragma omp for
        for (int i=0;i<NonZeroSolvent_.size();i++)
        {
            int solvIndex = NonZeroSolvent_[i];
            Real3 possv = solvent[solvIndex].positions_;
            Real qi     = solvent[solvIndex].charge_;

            for (int j=0;j<NonZeroSolute_.size();j++)
            {
                int soluIndex = NonZeroSolute_[j];
                Real3 possl = solute[soluIndex].positions_;
                Real qj     = solute[soluIndex].charge_;
                Real3 dist;
                Real distsq;

                simstate_.getSimulationBox().calculateDistance(possv, possl, dist, distsq);

                Real r = std::sqrt(distsq);

                sum += factor_ * qi * qj * std::erfc(r * alpha_)  / r;
            }
        }

        #pragma omp critical
        {
            energy_ += sum;
        }
    }
}

void SRE::calculate()
{
    energy_ = 0.0;

    if (mode_ == "NS")
    {
        calculateWithNS();
    }

    if (mode_ == "Brute")
    {
        calculateWithoutNS();
    }

    std::cout << "Energy = " << energy_ << "\n";
}

void SRE::printEnergyPerIter(std::ofstream& ofs)
{
    ofs << energy_ << "\n";
}