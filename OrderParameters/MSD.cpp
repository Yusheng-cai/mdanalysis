#include "MSD.h"

namespace CalculationRegistry
{
    registry_<MSD> registerMSD("MSD");
}

MSD::MSD(const CalculationInput& input)
: Calculation(input)
{
    pack_.ReadString("residue", ParameterPack::KeyType::Required, residueName_);
    pack_.ReadString("direction", ParameterPack::KeyType::Optional, directionStr_);

    auto it = MapStrToIndex2.find(directionStr_);
    ASSERT((it != MapStrToIndex2.end()), "The direction " << directionStr_ << " is not available.");
    directionIndex_ = it -> second;
    directionSize = directionIndex_.size();

    // add the residue group
    addResidueGroup(residueName_);
    auto& resgroup = getResidueGroup(residueName_);
    numAtoms_ = resgroup.getAtomSize();

    // register outputs
    registerOutputFunction("msd", [this](std::string name) -> void {this -> printMSD(name);});
}

void MSD::calculate()
{
    auto& resgroup = getResidueGroup(residueName_).getResidues();

    std::vector<Real3> atomPos;

    for (int i=0;i<resgroup.size();i++)
    {
        int numAtomsInres = resgroup[i].atoms_.size();
        for (int j=0;j<numAtomsInres;j++)
        {
            atomPos.push_back(resgroup[i].atoms_[j].positions_);
        }
    }
    positions_.push_back(atomPos);
}

void MSD::finishCalculate()
{
    int numFrames = positions_.size();
    std::vector<std::vector<Real>> MSDtimeframe(numFrames, std::vector<Real>(numAtoms_,0.0));

    for (int i=0;i<numAtoms_;i++)
    {
        std::vector<Real3> atomPos(numFrames);
        for (int j=0;j<numFrames;j++)
        {
            atomPos[j] = positions_[j][i];
        }

        std::vector<Real> MSD_atom = calculateMSD(atomPos);

        for (int j=0;j<numFrames;j++)
        {
            MSDtimeframe[j][i] = MSD_atom[j];
        }
    }

    // finally average over all the MSD for all the atoms 
    MSD_.resize(numFrames,0.0);
    for (int i=0;i<numFrames;i++)
    {
        for (int j=0;j<numAtoms_;j++)
        {
            MSD_[i] += MSDtimeframe[i][j];
        }
        MSD_[i] /= numAtoms_;
    }
}

void MSD::printMSD(std::string name)
{
    std::ofstream ofs;

    ofs.open(name);

    for (int i=0;i<MSD_.size();i++)
    {
        ofs << MSD_[i] << "\n";
    }

    ofs.close();
}

std::vector<MSD::Real> MSD::calculateMSD(std::vector<Real3>& position)
{
    int size = position[0].size();

    ASSERT((size == directionSize), "The position passed into calculating MSD is " << size << " while the required size is " << directionSize);

    int numFrames = position.size();

    // make a vector of squares 
    std::vector<Real> squares(numFrames, 0.0);
    std::vector<std::vector<Real>> data(numFrames);
    std::vector<std::vector<Real>> AC_vector(size);

    // calculate the squared distances
    for (int i=0;i<numFrames;i++)
    {
        data[i].resize(directionSize);
        for (int j=0;j<directionSize;j++)
        {
            int index = directionIndex_[j];
            squares[i] += std::pow(position[i][index],2.0);
            data[i][j] = position[i][index];
        }
    }

    // calculate the autocorrelation in each dimension
    FFT::autocorrelation(data, AC_vector, true);

    // create the SAB vector
    std::vector<Real> SAB(numFrames);
    for (int i=0;i<numFrames;i++)
    {
        for (int j=0;j<directionSize;j++)
        {
            SAB[i] += AC_vector[j][i];
        }
    }

    Real sumSQ=0.0;
    for (int i=0;i<numFrames-1;i++)
    {
        sumSQ += squares[i];
        std::cout << "squares " << i << " = " << squares[i] << std::endl;
    }
    sumSQ = 2*sumSQ;
    std::cout << "sumSQ = " << sumSQ << std::endl;

    std::vector<Real> MSD(numFrames,0.0);

    for (int i=1;i<numFrames-1;i++)
    {
        sumSQ = sumSQ - squares[i-1] - squares[numFrames-i];

        std::cout << "sumSQ = " << sumSQ << std::endl;

        MSD[i] =  sumSQ/(numFrames-i) - 2*SAB[i];

        std::cout << "MSD " << i << " = " << MSD[i] << std::endl;
    }

    return MSD;
}