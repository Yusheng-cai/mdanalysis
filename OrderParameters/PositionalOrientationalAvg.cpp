#include "PositionalOrientationalAvg.h"

namespace CalculationRegistry
{
    registry_<PositionalOrientationalAvg> registerPosOri("PosOriAvg");
}

PositionalOrientationalAvg::PositionalOrientationalAvg(const CalculationInput& input)
: Calculation(input)
{
    pack_.ReadString("residue", ParameterPack::KeyType::Required, ResidueName_);
    pack_.ReadNumber("headindex", ParameterPack::KeyType::Optional, HeadIndex_);
    pack_.ReadNumber("tailindex", ParameterPack::KeyType::Optional, TailIndex_);
    pack_.ReadString("average_type", ParameterPack::KeyType::Required, Avg_type_);

    if (Avg_type_ == "Usr")
    {
        pack_.ReadVectorNumber("Usr_indices", ParameterPack::KeyType::Required, UsrIndices_);
        pack_.ReadNumber("cutoff", ParameterPack::KeyType::Optional, Usr_cuttoff_);
        pack_.ReadNumber("epsilon", ParameterPack::KeyType::Optional, Usr_epsilon_);
        Usr_alpha_ = 1.0 / Usr_epsilon_;

        for (int i=0;i<UsrIndices_.size();i++)
        {
            UsrIndices_[i] = UsrIndices_[i] - 1;
        }
    }
    HeadIndex_--;
    TailIndex_--;

    // read the bins
    auto rbinPack = pack_.findParamPack("rbin", ParameterPack::KeyType::Required);
    RBin_ = Binptr(new Bin(*rbinPack));
    auto CosThetaBinPack = pack_.findParamPack("tbin", ParameterPack::KeyType::Required);
    CosThetaBin_ = Binptr(new Bin(*CosThetaBinPack));
    NumRBins_ = RBin_->getNumbins();
    NumCosThetaBins_ = CosThetaBin_->getNumbins();

    initializeResidueGroup(ResidueName_);

    PosOriHist_.resize(NumRBins_, std::vector<Real>(NumCosThetaBins_, 0.0));
    PosOriAvg_.resize(NumRBins_, std::vector<Real>(NumCosThetaBins_,0.0));

    // register the calculation functions 
    RegisterCalculationFunctions("Usr", [this](const Molecule::residue& res1, const Molecule::residue& res2)-> Real {return this -> PairUsr(res1,res2);});

    // register output functions
    registerOutputFunction("Average", [this](std::string name) -> void {this -> printAverage(name);});
}

void PositionalOrientationalAvg::printAverage(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    ofs << "# i\tj\tRbin\tCosThetaBin\tAverage\tNumberDensity\n";
    for (int i=0;i<NumRBins_;i++)
    {
        for (int j=0;j<NumCosThetaBins_;j++)
        {
            ofs << i << " " << j << " " << RBin_->getCenterLocationOfBin(i) << " " << CosThetaBin_->getCenterLocationOfBin(j) \
            << " " << PosOriAvg_[i][j] << " " << PosOriHist_[i][j] << "\n";
        }
    }

    ofs.close();
}

PositionalOrientationalAvg::Real PositionalOrientationalAvg::PairUsr(const Molecule::residue& res1, const Molecule::residue& res2)
{
    int numAtoms = res1.atoms_.size();

    Real energy  = 0.0;
    for (int i=0;i<UsrIndices_.size();i++)
    {
        int index1 = UsrIndices_[i];
        for (int j=0;j<UsrIndices_.size();j++)
        {
            int index2 = UsrIndices_[j];

            Real3 distance;
            Real dist_sq;

            Real3 Atom1Pos = res1.atoms_[index1].positions_;
            Real3 Atom2Pos = res2.atoms_[index2].positions_;

            simstate_.getSimulationBox().calculateDistance(Atom1Pos, Atom2Pos, distance, dist_sq);

            Real dist = std::sqrt(dist_sq);

            Real qiqj = res1.atoms_[index1].charge_ * res2.atoms_[index2].charge_;

            if (dist <= Usr_cuttoff_)
            {
                energy += qiqj * Usr_factor_ * std::erfc(Usr_alpha_ * dist) / dist; 
            }
        }
    }

    return energy;
}

void PositionalOrientationalAvg::RegisterCalculationFunctions(std::string name, avg_func func)
{
    auto it = MapNameToFunc_.find(name);

    ASSERT((it == MapNameToFunc_.end()), "The function name " << name << " is registered more than once.");

    MapNameToFunc_.insert(std::make_pair(name, func));
}

void PositionalOrientationalAvg::calculate()
{
    auto& res = getResidueGroup(ResidueName_).getResidues();

    uij_.clear();
    uij_.resize(res.size());

    // First find the COM
    for (int i=0;i<res.size();i++)
    {
        COM_[i] = calcCOM(res[i]);

        Real3 HeadPos = res[i].atoms_[HeadIndex_].positions_;
        Real3 TailPos = res[i].atoms_[TailIndex_].positions_;
        Real3 director;
        Real norm;

        simstate_.getSimulationBox().calculateDistance(HeadPos, TailPos, director, norm);
        LinAlg3x3::normalize(director);
        uij_[i] = director;
    }

    #pragma omp parallel
    {
        std::vector<std::vector<Real>> LocalHist(NumRBins_, std::vector<Real>(NumCosThetaBins_,0.0));
        std::vector<std::vector<Real>> LocalAvg(NumRBins_, std::vector<Real>(NumCosThetaBins_, 0.0));
        #pragma omp for
        for (int i=0;i<res.size();i++)
        {
            for (int j=i+1;j<res.size();j++)
            {
                Real3 distance;
                Real dist_sq;
                simstate_.getSimulationBox().calculateDistance(COM_[i], COM_[j], distance, dist_sq);
                Real dist = std::sqrt(dist_sq);
                Real CosTheta = LinAlg3x3::DotProduct(uij_[i], uij_[j]);

                // now we bin
                if (RBin_->isInRange(dist) && CosThetaBin_->isInRange(CosTheta))
                {
                    int Rnum = RBin_->findBin(dist);
                    int Tnum = CosThetaBin_->findBin(CosTheta);

                    LocalHist[Rnum][Tnum] += 1;

                    // calculate the average
                    auto it = MapNameToFunc_.find(Avg_type_);
                    ASSERT((it != MapNameToFunc_.end()), "Averaging type " << Avg_type_ << " is not supported.");
                    Real Value = it -> second(res[i], res[j]);

                    LocalAvg[Rnum][Tnum] += Value;
                }
            }
        }

        #pragma omp critical
        for (int i=0;i<NumRBins_;i++)
        {
            for (int j=0;j<NumCosThetaBins_;j++)
            {
                PosOriHist_[i][j] += LocalHist[i][j];
                PosOriAvg_[i][j]  += LocalAvg[i][j];
            }
        }
    }
}

void PositionalOrientationalAvg::finishCalculate()
{
    int numframes = simstate_.getTotalFrames();

    for (int i=0;i<NumRBins_;i++)
    {
        for (int j=0;j<NumCosThetaBins_;j++)
        {
            if (PosOriHist_[i][j] != 0)
            {
                PosOriAvg_[i][j] = PosOriAvg_[i][j] / PosOriHist_[i][j];
            }
        }
    }

    for (int i=0;i<NumRBins_;i++)
    {
        for (int j=0;j<NumCosThetaBins_;j++)
        {
            PosOriHist_[i][j] = PosOriHist_[i][j] / numframes;
        }
    }
}