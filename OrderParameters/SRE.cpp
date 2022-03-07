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

    pack_.Readbool("onlyattractive", ParameterPack::KeyType::Optional, onlyattrative_);
    ASSERT((mode_ == "NS" || mode_ == "Brute"), "The mode can only be Brute or NS.");

    cutoffsq_ = cutoff_ * cutoff_;

    addResidueGroup(SolventName_);
    addResidueGroup(SoluteName_);

    // register the printing function for energy 
    registerPerIterOutputFunction("energy", [this](std::ofstream& ofs) -> void {this -> printEnergyPerIter(ofs);});
    registerPerIterOutputFunction("energyperatom", [this](std::ofstream& ofs) -> void {this -> printEnergyPerAtomPerIter(ofs);});
    registerOutputFileOutputs("energy", [this](void) -> Real {return this -> getEnergy();});

    cell_ = Cellptr(new CellGrid(simstate_, cutoff_,1));

    initializeProbeVolumes();
    initializeNotInProbeVolumes();

    // initialize the residue indices for solvent and solute
    initializeSoluteSolventIndices();
}

void SRE::initializeSoluteSolventIndices()
{
    const auto& solvent = getResidueGroup(SolventName_).getResidues();
    const auto& solute  = getResidueGroup(SoluteName_).getResidues();

    int solventResSize = solvent[0].atoms_.size();
    int soluteResSize  = solute[0].atoms_.size();

    // obtain per residue indices --> [atom1, atom2, ...]
    std::vector<int> SoluteIndicesPerResidue;
    std::vector<int> SolventIndicesPerResidue;
    SoluteIndicesPerResidue.resize(soluteResSize);
    SolventIndicesPerResidue.resize(solventResSize);
    std::iota(SoluteIndicesPerResidue.begin(), SoluteIndicesPerResidue.end(), 0);
    std::iota(SolventIndicesPerResidue.begin(), SolventIndicesPerResidue.end(), 0);

    // read in the actual indices 
    std::vector<std::string> soluteInd_;
    pack_.ReadVectorString("soluteIndices", ParameterPack::KeyType::Optional,soluteInd_);
    pack_.ReadVectorNumber("solventIndices", ParameterPack::KeyType::Optional, SolventIndicesPerResidue);
    StringTools::ConvertStringToIndices(soluteInd_, SoluteIndicesPerResidue);

    // make it into a whole vector of local indices 
    for (int i=0;i<solvent.size();i++)
    {
        for (int j=0;j<SolventIndicesPerResidue.size();j++)
        {
            SolventIndices_.push_back(i * solventResSize + SolventIndicesPerResidue[j]);
        }
    }

    for (int i=0;i<solute.size();i++)
    {
        for (int j=0;j<SoluteIndicesPerResidue.size();j++)
        {
            SoluteIndices_.push_back(i* soluteResSize + SoluteIndicesPerResidue[j]);
        }
    }
}

void SRE::getInsidePVIndices()
{
    const auto& solvent  = getResidueGroup(SolventName_).getTotalResidue().atoms_;
    InsidePVIndices_.clear();
    indus_indicators_.clear();

    for (int i=0;i<SolventIndices_.size();i++)
    {
        int index = SolventIndices_[i];
        Real3 pos = solvent[index].positions_;
        Real htildex=1.0;

        if (isInPV(pos, htildex))
        {
            InsidePVIndices_.push_back(index);
            indus_indicators_.push_back(htildex);
        }
    }
}

void SRE::getNonZeroCharges()
{
    const auto& solvent = getResidueGroup(SolventName_).getTotalResidue().atoms_;
    const auto& solute  = getResidueGroup(SoluteName_).getTotalResidue().atoms_;

    NonZeroSolvent_.clear();
    NonZeroSolute_.clear();
    nonzero_indus_indicators_.clear();

    for (int i=0;i<InsidePVIndices_.size();i++)
    {
        int index = InsidePVIndices_[i];
        Real charge = solvent[index].charge_;
        Real htildex= indus_indicators_[i];

        if (charge != 0)
        {
            NonZeroSolvent_.push_back(index);
            nonzero_indus_indicators_.push_back(htildex);
        }
    }
    
    for (int i=0;i<SoluteIndices_.size();i++)
    {
        int index = SoluteIndices_[i];
        Real charge = solute[index].charge_;

        if (charge != 0)
        {
            NonZeroSolute_.push_back(index);
        }
    }
}

void SRE::update()
{
    cell_ -> update();

    getInsidePVIndices();
    getNonZeroCharges();

    PerAtomContribution_.clear();
}

void SRE::calculateWithNS()
{
    const auto& solvent = getResidueGroup(SolventName_).getTotalResidue().atoms_;
    const auto& solute  = getResidueGroup(SoluteName_).getTotalResidue().atoms_;


    PerAtomContribution_.resize(solvent.size());

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

            Real localsum = 0.0;

            for (int ind : Neighbors)
            {
                ASSERT((ind < num), "Ind is " << ind << " max is " << num);
                for (int index : celllist_[ind])
                {
                    Real3 possl = solute[index].positions_;
                    Real qj     = solute[index].charge_;
                    Real qiqj = qi * qj;

                    if (onlyattrative_)
                    {
                        if (qiqj > 0)
                        {
                            continue;
                        }
                    }

                    Real3 dist;
                    Real distsq;

                    simstate_.getSimulationBox().calculateDistance(possv, possl, dist, distsq);

                    if (distsq < cutoffsq_)
                    {
                        Real r = std::sqrt(distsq);

                        Real value = factor_ * qiqj * std::erfc(r * alpha_) / r * nonzero_indus_indicators_[i];

                        sum += value;
                        localsum += value;
                    }
                }
            }

            PerAtomContribution_[solventInd] = localsum;
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

    PerAtomContribution_.clear();
    PerAtomContribution_.resize(solvent.size(),0.0);

    #pragma omp parallel
    {
        Real sum = 0.0;
        #pragma omp for
        for (int i=0;i<NonZeroSolvent_.size();i++)
        {
            int solvIndex = NonZeroSolvent_[i];
            Real3 possv = solvent[solvIndex].positions_;
            Real qi     = solvent[solvIndex].charge_;

            Real localSum = 0.0;

            for (int j=0;j<NonZeroSolute_.size();j++)
            {
                int soluIndex = NonZeroSolute_[j];
                Real3 possl = solute[soluIndex].positions_;
                Real qj     = solute[soluIndex].charge_;
                Real3 dist;
                Real distsq;

                simstate_.getSimulationBox().calculateDistance(possv, possl, dist, distsq);

                Real r = std::sqrt(distsq);
                Real value = factor_ * qi * qj * std::erfc(r * alpha_)  / r; 

                sum += value;
                localSum += value;
            }

            PerAtomContribution_[solvIndex] = localSum;
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

void SRE::printEnergyPerAtomPerIter(std::ofstream& ofs)
{
    int timestep  =simstate_.getTime();
    ofs << timestep << " ";

    for (int i=0;i<PerAtomContribution_.size();i++)
    {
        ofs << PerAtomContribution_[i] << " ";
    }

    ofs << "\n";
}