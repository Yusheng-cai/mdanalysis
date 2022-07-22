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

    // added a functionality to be able to rotate to some predefined direction
    // we always rotate such that the director is the z axis 
    rotate_ = pack_.ReadString("RotateToDefinedDirector", ParameterPack::KeyType::Optional, RotationMode_);
    if (rotate_)
    {
        if (RotationMode_ == "usedirector")
        {
            UseDirector_=true;
            pack_.ReadNumber("headindex", ParameterPack::KeyType::Optional, headindex_);
            pack_.ReadNumber("tailindex", ParameterPack::KeyType::Optional, tailindex_);
            headindex_--;
            tailindex_--;
        }
        else if (RotationMode_ == "array")
        {
            pack_.ReadArrayNumber("array", ParameterPack::KeyType::Required, array_);
            UseDirector_=false;
        }
        else
        {
            ASSERT((true==false), "Rotation mode " << RotationMode_ << " is not supported, only usedirection/array is supported.");
        }
    }
}

void MSD::calculate()
{
    const auto& resgroup = getResidueGroup(residueName_).getResidues();

    // are we rotating to some other frame of reference
    if (rotate_)
    {
        if (UseDirector_)
        {
            Matrix Qtensor = {};
            #pragma omp parallel
            {
                Matrix Qlocal = {};
                #pragma omp for
                for (int i=0;i<resgroup.size();i++)
                {
                    Real3 headposition = resgroup[i].atoms_[headindex_].positions_;
                    Real3 tailposition = resgroup[i].atoms_[tailindex_].positions_;
                    Real3 r;
                    Real rsq;
                    simstate_.getSimulationBox().calculateDistance(headposition, tailposition, r, rsq);

                    LinAlg3x3::normalize(r);
                    Matrix singleQ = LinAlg3x3::LocalQtensor(r);

                    LinAlg3x3::matrix_accum_inplace(Qlocal, singleQ);
                }

                #pragma omp critical
                LinAlg3x3::matrix_accum_inplace(Qtensor, Qlocal);
            }

            LinAlg3x3::matrix_mult_inplace(Qtensor, 1.0/(2.0 * resgroup.size()));
            auto res = LinAlg3x3::OrderEigenSolver(Qtensor);
            for (int i=0;i<3;i++)
            {
                array_[i] = res.second[i][0];
            }
        }

        // matrix that rotates array to [0,0,1] which is the z axis 
        RotationMatrix_ = LinAlg3x3::GetRotationMatrix(array_, {{0,0,1}});
    }

    // numresidue,3
    std::vector<Real3> COM;
    COM.resize(resgroup.size());

    #pragma omp parallel for 
    for (int i=0;i<resgroup.size();i++)
    {
        Real3 c = calcCOM(resgroup[i], COMIndices_);

        if (rotate_)
        {
            c = LinAlg3x3::MatrixDotVector(RotationMatrix_, c);
        }

        COM[i] = c;
    }

    Real3 sides = simstate_.getSimulationBox().getSides();

    positions_.push_back(COM);
    boxSides_.push_back(sides);
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

        // unwrap the positions
        for (int j=1;j<numFrames;j++)
        {
            Real3 shift = simstate_.getSimulationBox().calculateShift(ResiduePos[j], ResiduePos[j-1], boxSides_[j]);
            for (int k=0;k<3;k++)
            {
                ResiduePos[j][k] += shift[k];
            }
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
    Real dt = timestamps[1] - timestamps[0];

    std::ofstream ofs;

    ofs.open(name);

    int numdir = directionsIndex_.size();
    int numframes = MSD_[0].size();

    // write the header
    ofs << "# lagtime[ps]\t";
    for (auto dir : VectorDirection_)
    {
        ofs << dir << "\t";
    }
    ofs << "\n";

    for (int i=0;i<numframes;i++)
    {
        ofs << dt * i << "\t";
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
                data[i][j][k] = position[j][index];
                squares[i][j] += std::pow(data[i][j][k],2.0);
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