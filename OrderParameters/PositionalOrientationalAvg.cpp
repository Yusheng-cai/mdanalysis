#include "PositionalOrientationalAvg.h"

PositionalOrientationalAvg::PositionalOrientationalAvg(const CalculationInput& input)
: Calculation(input)
{
    pack_.ReadString("residue", ParameterPack::KeyType::Required, ResidueName_);
    pack_.ReadNumber("headindex", ParameterPack::KeyType::Optional, HeadIndex_);
    pack_.ReadNumber("tailindex", ParameterPack::KeyType::Optional, TailIndex_);
    pack_.ReadString("average_type", ParameterPack::KeyType::Required, Avg_type_);
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
}

void PositionalOrientationalAvg::RegisterCalculationFunctions()
{

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
            PosOriAvg_[i][j] = PosOriAvg_[i][j] / PosOriHist_[i][j];
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