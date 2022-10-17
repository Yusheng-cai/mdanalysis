#include "EndToEndDistance.h"

namespace CalculationRegistry
{
    registry_<EndToEndDistance> registerEndToEndDistance("EndToEndDistance");
}

EndToEndDistance::EndToEndDistance(const CalculationInput& input)
: Calculation(input)
{
    pack_.ReadString("residue", ParameterPack::KeyType::Optional, resname_);
    pack_.ReadNumber("headindex", ParameterPack::KeyType::Required, head_index_);
    pack_.ReadNumber("tailindex", ParameterPack::KeyType::Required, tail_index_);
    head_index_--;
    tail_index_--;

    initializeResidueGroup(resname_);

    auto binPack = pack_.findParamPack("bin", ParameterPack::KeyType::Required);
    bin_ = binptr(new Bin(*binPack));
    numbins_ = bin_->getNumbins();
    histogram_.resize(numbins_, 0.0);

    registerOutputFunction("histogram", [this](std::string name)->void {this -> printHistogram(name);});
}

void EndToEndDistance::calculate(){
    auto residues = getResidueGroup(resname_).getResidues();

    for (int i=0;i<residues.size();i++){
        auto& res = residues[i];
        Real3 headpos = res.atoms_[head_index_].positions_;
        Real3 tailpos = res.atoms_[tail_index_].positions_;
        Real r_squared;
        Real3 diff;

        simstate_.getSimulationBox().calculateDistance(headpos, tailpos, diff, r_squared);

        Real r = std::sqrt(r_squared);

        if (bin_->isInRange(r)){
            int index = bin_->findBin(r);
            histogram_[index] += 1;
        }
    }
}

void EndToEndDistance::finishCalculate()
{
    int numframes = simstate_.getTotalFrames();

    for (int i=0;i<numbins_;i++)
    {
        histogram_[i] = histogram_[i] / numframes;
    }
}

void EndToEndDistance::printHistogram(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    ofs << "# R(nm) Count\n";

    for (int i=0;i<numbins_;i++)
    {
        ofs << bin_->getCenterLocationOfBin(i) << " " << histogram_[i] << "\n";
    }
    ofs.close();
}
