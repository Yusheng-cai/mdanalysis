#include "Pcost.h"

namespace CalculationRegistry
{
    static registry_<Pcost> registerPcost("pcost");
}

Pcost::Pcost(const CalculationInput& input)
:Calculation(input)
{
    auto costBinPack = input.pack_.findParamPack("tbin", ParameterPack::KeyType::Required);
    costBin_ = Binptr(new Bin(*costBinPack));

    input.pack_.ReadString("residue", ParameterPack::KeyType::Required, residueGroupName_);
    input.pack_.ReadNumber("headindex", ParameterPack::KeyType::Required,headIndex_);
    input.pack_.ReadNumber("tailindex", ParameterPack::KeyType::Required,tailIndex_);

    input.pack_.ReadArrayNumber("array", ParameterPack::KeyType::Optional, arr_);
    LinAlg3x3::normalize(arr_);
    input.pack_.ReadNumber("precision", ParameterPack::KeyType::Optional, precision_);

    registerOutputFunction("histogram", [this](std::string name) -> void { this -> printHistogram(name);});
    registerOutputFunction("AtomIndices", [this](std::string name) -> void { this -> printAtomIndices(name);});

    registerPerIterOutputFunction("costheta", [this](std::ofstream& ofs) -> void { this -> printcosthetaHistogramPerIter(ofs);});

    headIndex_--;
    tailIndex_--;

    initializeResidueGroup(residueGroupName_);
    auto res = getResidueGroup(residueGroupName_).getResidues();
    std::vector<Real> Zero(costBin_->getNumbins(),0);

    // resize the molecular director size
    uij_.resize(res.size());

    // resize of histogram to size of bins 
    histogram_.resize(costBin_->getNumbins());
    std::fill(histogram_.begin(), histogram_.end(), 0.0);

    initializeNotInProbeVolumes();
    initializeProbeVolumes();
}

void Pcost::calculate()
{
    auto& res = getResidueGroup(residueGroupName_).getResidues();

    // resize center of mass
    COM_.resize(res.size());

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
    InsideIndices_ = InsidePVIndices(COM_);

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
    for (int i=0;i<InsideIndices_.size();i++)
    {
        int k = InsideIndices_[i];
        Real cost = Qtensor::vec_dot(uij_[k], arr_);

        // ASSERT((cost >= -1 && cost <= 1), "cosine(theta) is not within range of -1 and 1");

        int binNum = costBin_->findBin(cost);
        ASSERT((binNum <= histogram_.size()-1), "Bin number is out of range of histogram.");

        histogram_[binNum] += 1; 
        histogramPerIter_[binNum] += 1;
    }
}

void Pcost::printcosthetaHistogramPerIter(std::ofstream& ofs)
{
    for (int i=0;i<histogramPerIter_.size();i++)
    {
        ofs << histogramPerIter_[i] << "\t";
    }
    ofs << "\n";
}

void Pcost::finishCalculate()
{
    int Numframes = simstate_.getTotalFrames();

    for (int i=0;i<histogram_.size();i++)
    {
        histogram_[i] /= Numframes;
    }
}

void Pcost::printHistogram(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    ofs << std::fixed << std::setprecision(precision_);

    for (int i=0;i<histogram_.size();i++)
    {
        ofs << i << "\t" << histogram_[i] << "\n";
    }
    ofs.close();
}

void Pcost::printAtomIndices(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    ofs << std::fixed << std::setprecision(precision_);

    for (int i=0;i<AtomIndicesInPV_.size();i++)
    {
        int size = AtomIndicesInPV_[i].size();

        for (int j=0;j<size;j++)
        {
            ofs << AtomIndicesInPV_[i][j] << "\t";
        }
        ofs << "\n";
    }

    ofs.close();
}
