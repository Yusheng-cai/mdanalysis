#include "Pcost2.h"

namespace CalculationRegistry
{
    static registry_<Pcost2> registerPcost2("pcost2");
}

Pcost2::Pcost2(const CalculationInput& input)
:Pcost(input)
{

}


void Pcost2::calculate()
{
    auto& res = getResidueGroup(residueGroupName_).getResidues();
    auto& pv  = simstate_.getProbeVolume(ProbeVolumeName_);

    // clear the histogram per iteration
    histogramPerIter_.clear();
    histogramPerIter_.resize(costBin_->getNumbins(),0.0);

    // find all the COM of the residues in the system
    #pragma omp parallel for
    for (int i=0;i<res.size();i++)
    {
        Real3 com = CalculationTools::getCOM(res[i], simstate_, COMIndices_);
        COM_[i] = com;

        Real3 distance_;
        Real dist_sq;
        Real3 headPos_ = res[i].atoms_[headIndex_].positions_;
        Real3 tailPos_ = res[i].atoms_[tailIndex_].positions_;
        simstate_.getSimulationBox().calculateDistance(headPos_, tailPos_, distance_, dist_sq);

        Real3 normalized_dir = Qtensor::normalize_director(distance_);
        uij_[i] = normalized_dir;
    }

    // check which COM are inside the probevolume
    InsideIndices_.clear();
    for (int i=0;i<COM_.size();i++)
    {
        auto pvOutput = pv.calculate(COM_[i]);
        bool excluded=false;

        for (auto pv2 : NotInprobevolumes_) 
        {
            auto pv0 = pv2 -> calculate(COM_[i]);

            if (pv0.hx_ == 1)
            {
                excluded=true;
            }
        }

        if (! excluded)
        {
            if (pvOutput.hx_ == 1)
            {
                InsideIndices_.push_back(i);
            }
        }
    }

    std::cout << "Insideindices.size = " << InsideIndices_.size() << std::endl;

    std::vector<int> AtomIndicesINPVIter;
    // get the atom indices in the pv per iteration
    for (int i=0;i<InsideIndices_.size();i++)
    {
        for (int j=0;j<res[i].atoms_.size();j++)
        {
            AtomIndicesINPVIter.push_back(res[i].atoms_[j].atomNumber_);
        }
    }
    AtomIndicesInPV_.push_back(AtomIndicesINPVIter);
    
    // starting binning 
    p2tilde_ = 0.0;
    for (int i=0;i<InsideIndices_.size();i++)
    {
        int k = InsideIndices_[i];
        Real cost = std::pow(Qtensor::vec_dot(uij_[k], arr_),2.0);

        p2tilde_ = p2tilde_ + 1.5 * cost - 0.5;

        // ASSERT((cost >= -1 && cost <= 1), "cosine(theta) is not within range of -1 and 1");

        int binNum = costBin_->findBin(cost);
        ASSERT((binNum <= histogram_.size()-1), "Bin number is out of range of histogram.");

        histogram_[binNum] += 1; 
        histogramPerIter_[binNum] += 1;
    }

    p2tilde_ /= InsideIndices_.size();
}