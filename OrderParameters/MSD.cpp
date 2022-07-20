#include "MSD.h"

namespace CalculationRegistry
{
    registry_<MSD> registerMSD("MSD");
}

MSD::MSD(const CalculationInput& input)
: Calculation(input)
{
    pack_.ReadString("residue", ParameterPack::KeyType::Required, residueName_);

    // the number of directions 
    pack_.ReadVectorString("directions", ParameterPack::KeyType::Required, VectorDirection_);

    for (int i=0;i<VectorDirection_.size();i++)
    {
        auto it = MapStrToIndex2.find(VectorDirection_[i]);
        ASSERT((it != MapStrToIndex2.end()), "The direction " << VectorDirection_[i] << " is not available.");

        directionsIndex_.push_back(it -> second);
    }

    // add the residue group
    initializeResidueGroup(residueName_);
    const auto& res = getResidueGroup(residueName_);
    numResidues_ = res.getResidues().size();

    // register outputs
    registerOutputFunction("msd", [this](std::string name) -> void {this -> printMSD(name);});
}

void MSD::calculate()
{
    auto& resgroup = getResidueGroup(residueName_).getResidues();

    // numresidue,3
    std::vector<Real3> COM;
    COM.resize(resgroup.size());

    #pragma omp parallel for 
    for (int i=0;i<resgroup.size();i++)
    {
        COM[i] = calcCOM(resgroup[i], COMIndices_);
    }

    positions_.push_back(COM);
}

void MSD::finishCalculate()
{
    int numFrames = positions_.size();
    int numdir    = directionsIndex_.size();

    // numdir, numFrames, numResidues
    std::vector<std::vector<std::vector<Real>>> MSDtimeframe(numdir, std::vector<std::vector<Real>>(numFrames, std::vector<Real>(numResidues_,0.0)));

    for (int i=0;i<numResidues_;i++)
    {
        // store the residue positions for one particular residue (timeframe, 3)
        std::vector<Real3> ResiduePos(numFrames);
        for (int j=0;j<numFrames;j++)
        {
            ResiduePos[j] = positions_[j][i];
        }

        // (numdir, numFrames) calculate the mean squared displacement
        std::vector<std::vector<Real>> MSD_residue = calculateMSD(ResiduePos);

        // copy over the mean squared displacement per residue to a larger array 
        for (int j=0;j<numdir;j++)
        {
            for (int k=0;k<numFrames;k++)
            {
                MSDtimeframe[j][k][i] = MSD_residue[j][k];
            }
        }
    }

    // finally average over all the MSD for all the atoms 
    MSD_.resize(numdir, std::vector<Real>(numFrames,0.0));

    for (int i=0;i<numdir;i++)
    {
        for (int j=0;j<numFrames;j++)
        {
            for (int k=0;k<numResidues_;k++)
            {
                MSD_[i][j] += MSDtimeframe[i][j][k];
            }
            MSD_[i][j] /= numResidues_;
        }
    }
}

void MSD::printMSD(std::string name)
{
    std::vector<Real> timestamps = simstate_.gettimestamps();

    std::ofstream ofs;

    ofs.open(name);

    int numdir = directionsIndex_.size();
    int numframes = MSD_[0].size();

    // write the header
    ofs << "# time[ps]\t";
    for (auto dir : VectorDirection_)
    {
        ofs << dir << "\t";
    }
    ofs << "\n";

    for (int i=0;i<numframes;i++)
    {
        ofs << timestamps[i] << "\t";
        for (int j=0;j<numdir;j++)
        {
            ofs << MSD_[j][i] << "\t";
        }
        ofs << "\n";
    }

    ofs.close();
}

std::vector<std::vector<MSD::Real>> MSD::calculateMSD(std::vector<Real3>& position)
{
    int numFrames = position.size();
    int numdir = directionsIndex_.size();

    // make a vector of squares 
    // which direction, numFrames
    std::vector<std::vector<Real>> squares(numdir, std::vector<Real>(numFrames, 0.0));

    // set up the MSD 
    // which direction, numFrames
    std::vector<std::vector<Real>> MSD(numdir, std::vector<Real>(numFrames,0.0));

    // which direction, numFrames,  directionIndex
    std::vector<std::vector<std::vector<Real>>> data(numdir);

    // calculate the squared distances
    for (int i=0;i<numdir;i++)
    {
        int directionSize = directionsIndex_[i].size();
        data[i].resize(numFrames, std::vector<Real>(directionSize,0.0));
        for (int j=0;j<numFrames;j++)
        {
            for (int k=0;k<directionSize;k++)
            {
                int index = directionsIndex_[i][k];
                squares[i][j] += std::pow(position[j][index],2.0);
                data[i][j][k] = position[j][index];
            }
        }
    }

    for (int i=0;i<directionsIndex_.size();i++)
    {
        int directionSize = directionsIndex_[i].size();
        std::vector<std::vector<Real>> AC_vector;

        ASSERT((data[i][0].size() == directionSize), "The size of the data does not match the size of the directions.");

        // calculate the autocorrelation in each dimension
        FFT::autocorrelation(data[i], AC_vector);

        // create the SAB vector
        std::vector<Real> SAB(numFrames);
        for (int j=0;j<numFrames;j++)
        {
            for (int k=0;k<directionSize;k++)
            {
                SAB[j] += AC_vector[j][k];
            }
        }

        // sum up the squares
        Real sumSQ=0.0;
        for (int j=0;j<numFrames;j++)
        {
            sumSQ += squares[i][j];
        }
        sumSQ = 2*sumSQ;

        // start calculating the mean squared diplacement
        MSD[i][0] = (sumSQ - 2 *SAB[0])/numFrames;

        for (int j=1;j<numFrames;j++)
        {
            sumSQ = sumSQ - squares[i][j-1] - squares[i][numFrames-j];
            MSD[i][j] =  (sumSQ - 2 * SAB[j])/(numFrames-j);
        }
    }

    return MSD;
}