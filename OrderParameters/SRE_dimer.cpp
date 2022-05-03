#include "SRE_dimer.h"


namespace CalculationRegistry
{
    registry_<SRE_dimer> registerSRE_dimer("SRE_dimer");
}

SRE_dimer::SRE_dimer(const CalculationInput& input)
: Calculation(input)
{
    // initialize the bins 
    const auto binPack = pack_.findParamPack("bin", ParameterPack::KeyType::Required);
    bin_ = Binptr(new Bin(*binPack));
    binnum_ = bin_->getNumbins();
    histogram_attr_.resize(binnum_,0.0);
    histogram_repu_.resize(binnum_,0.0);
    histogram_tota_.resize(binnum_,0.0);
    count_.resize(binnum_,0.0);

    // read residue 
    pack_.ReadString("residue", ParameterPack::KeyType::Required, residueName_);
    initializeResidueGroup(residueName_);

    // read head and tail index
    pack_.ReadNumber("headindex", ParameterPack::KeyType::Optional, head_index_);
    pack_.ReadNumber("tailindex", ParameterPack::KeyType::Optional, tail_index_);
    head_index_--;
    tail_index_--;

    // read Rlim 
    pack_.ReadNumber("Rlim", ParameterPack::KeyType::Required, Rlim_);
    Rlimsq_ = Rlim_ * Rlim_;

    // read sigma
    pack_.ReadNumber("sigma", ParameterPack::KeyType::Optional, sigma_);
    beta_ = 1.0/sigma_;

    // read if only CN 
    pack_.Readbool("onlyCN", ParameterPack::KeyType::Optional, onlyCN_);

    // register outputs
    registerOutputFunction("EnergyHistogram", [this](std::string name) -> void {this -> printHistograms(name);});
}

void SRE_dimer::calculate()
{
    const auto& res = getResidueGroup(residueName_).getResidues();
    COM_.clear();
    COM_.resize(res.size());
    director_.clear();
    director_.resize(res.size());

    // obtain the center of mass of each of the residues
    for (int i=0;i<res.size();i++)
    {
        Real3 COMperAtom_;
        COMperAtom_ = CalculationTools::getCOM(res[i], simstate_, COMIndices_);
        COM_[i] = COMperAtom_;

        Real3 headpos = res[i].atoms_[head_index_].positions_;
        Real3 tailpos = res[i].atoms_[tail_index_].positions_;
        Real3 distance;
        Real distsq;
        simstate_.getSimulationBox().calculateDistance(headpos, tailpos, distance, distsq);
        LinAlg3x3::normalize(distance);

        director_[i] = distance;
    }

    #pragma omp parallel
    {
        std::vector<Real> localAttr(binnum_,0.0);
        std::vector<Real> localRepu(binnum_,0.0);
        std::vector<Real> localtota(binnum_,0.0);
        std::vector<Real> localcount(binnum_,0.0);

        #pragma omp for
        for (int i=0;i<res.size();i++)
        {
            for (int j=i+1;j<res.size();j++)
            {
                Real3 distance;
                Real distsq;

                simstate_.getSimulationBox().calculateDistance(COM_[i], COM_[j], distance, distsq);

                if (distsq <= Rlimsq_)
                {
                    Real attr;
                    Real repul;
                    Real total;

                    calculateSRE(i,j,attr, repul, total);

                    Real dotProduct = LinAlg3x3::DotProduct(director_[i], director_[j]);
                    int binIndex = bin_->findBin(dotProduct);

                    localAttr[binIndex] += attr;
                    localRepu[binIndex] += repul;
                    localtota[binIndex] += total;
                    localcount[binIndex] += 1;
                }
            }
        }

        #pragma omp critical
        {
            for (int i=0;i<binnum_;i++)
            {
                histogram_attr_[i] += localAttr[i];
                histogram_repu_[i] += localRepu[i];
                histogram_tota_[i] += localtota[i];
                count_[i] += localcount[i];
            }
        }
    }
}

void SRE_dimer::finishCalculate()
{
    for (int i=0;i<binnum_;i++)
    {
        if (count_[i] != 0)
        {
            histogram_attr_[i] = histogram_attr_[i] / count_[i];
            histogram_repu_[i] = histogram_repu_[i] / count_[i];
            histogram_tota_[i] = histogram_tota_[i] / count_[i];
        }
    }
}

void SRE_dimer::printHistograms(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);
    int numFrames = simstate_.getTotalFrames();

    ofs << "# BinPos Attr Repul total count\n";

    for (int i=0;i<binnum_;i++)
    {
        ofs << bin_ -> getCenterLocationOfBin(i) << " ";
        ofs << histogram_attr_[i] << " ";
        ofs << histogram_repu_[i] << " ";
        ofs << histogram_tota_[i] << " ";
        ofs << count_[i]/numFrames << "\n";
    }

    ofs.close();
}

void SRE_dimer::calculateSRE(int i, int j, Real& attr, Real& repul, Real& total)
{
    const auto& res = getResidueGroup(residueName_).getResidues();

    auto& res1 = res[i].atoms_;
    auto& res2 = res[j].atoms_;
    attr = 0.0;
    repul = 0.0;

    if (onlyCN_)    
    {
        Int2 indices = {{head_index_, tail_index_}};

        for (int i=0;i<2;i++)
        {
            for (int j=0;j<2;j++)
            {
                Real distsq;
                Real3 distance;

                int index1 = indices[i];
                int index2 = indices[j];
                simstate_.getSimulationBox().calculateDistance(res1[index1].positions_, res2[index2].positions_, distance, distsq);

                Real qiqj = res1[index1].charge_ * res2[index2].charge_;
                Real dist = std::sqrt(distsq);
                Real erfval = std::erfc(dist * beta_);

                if (qiqj < 0)
                {
                    attr += factor_ * qiqj * erfval / dist; 
                }

                if (qiqj > 0)
                {
                    repul += factor_ * qiqj * erfval /dist;
                }
            }
        }
    }
    else
    {
        for (int i=0;i<res1.size();i++)
        {
            for (int j=0;j<res2.size();j++)
            {
                Real distsq;
                Real3 distance;
                simstate_.getSimulationBox().calculateDistance(res1[i].positions_, res2[j].positions_, distance, distsq);

                Real qiqj = res1[i].charge_ * res2[j].charge_;
                Real dist = std::sqrt(distsq);
                Real erfval = std::erfc(dist * beta_);

                if (qiqj < 0)
                {
                    attr += factor_ * qiqj * erfval / dist; 
                }

                if (qiqj > 0)
                {
                    repul += factor_ * qiqj * erfval /dist;
                }
            }
        }
    }

    total = attr + repul;
}