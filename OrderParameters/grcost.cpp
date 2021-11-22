#include "grcost.h"

grcost::grcost(const CalculationInput& input)
: Calculation(input)
{
    initializeProbeVolumes();
    initializeNotInProbeVolumes();

    auto tbinPack = input.pack_.findParamPack("tbin", ParameterPack::KeyType::Required);
    auto rbinPack = input.pack_.findParamPack("rbin", ParameterPack::KeyType::Required);

    // initialize tbin and rbin
    tbin_ = binptr(new Bin(*tbinPack));
    rbin_ = binptr(new Bin(*rbinPack));

    input.pack_.ReadString("residue", ParameterPack::KeyType::Required, residueName_);
    initializeResidueGroup(residueName_);

    input.pack_.ReadNumber("headindex", ParameterPack::KeyType::Optional, headindex_);
    input.pack_.ReadNumber("tailindex", ParameterPack::KeyType::Optional, tailindex_);
    headindex_--;
    tailindex_--;
}

void grcost::calculate()
{
    // first let's calculate all the COM
    auto& res = getResidueGroup(residueName_).getResidues();
    COM_.clear();
    COM_.resize(res.size());
    uij_.clear();
    uij_.resize(res.size());

    Qtensor_.fill({});

    #pragma omp parallel for
    for (int i=0;i<res.size();i++)
    {
        COM_[i] = calcCOM(res[i]);

        Real3 headpos = res[i].atoms_[headindex_].positions_;
        Real3 tailpos = res[i].atoms_[tailindex_].positions_;

        Real3 distance;
        Real distance_sq;

        simstate_.getSimulationBox().calculateDistance(headpos, tailpos, distance, distance_sq);

        LinAlg3x3::normalize(distance);

        uij_[i] = distance;
    }

    InsideIndices_ = InsidePVIndices(COM_);
    int N = InsideIndices_.size();

    // first let's calculate the Qtensor inside of the probe volume 
    #pragma omp parallel
    {
        Matrix localQ={};
        #pragma omp for
        for (int i=0;i<InsideIndices_.size();i++)
        {
            int index = InsideIndices_[i];

            Matrix singleQ = LinAlg3x3::LocalQtensor(uij_[index]);

            Qtensor::matrix_accum_inplace(localQ, singleQ);
        }

        #pragma omp critical
        {
            Qtensor::matrix_accum_inplace(Qtensor_, localQ);
        }
    }

    Qtensor::matrix_mult_inplace(Qtensor_, 1.0/(2.0*N));
    auto result = Qtensor::orderedeig_Qtensor(Qtensor_);
    // get the director out of the result
    for (int i=0;i<3;i++)
    {
        director_[i] = result.first[i][0];
    }

    neighborDistance_.clear();
    neighborDistance_.resize(N, std::vector<Real>(N,0));
    rdotdirector_.clear();
    rdotdirector_.resize(N, std::vector<Real>(N,0.0));
    histogramPerIter_.clear();
    histogramPerIter_.resize(N, std::vector<int>(N,0));

    for (int i=0;i<N;i++)
    {
        for (int j=i+1;j<N;j++)
        {
            Real3 distance;
            Real distance_sq;

            int index1 = InsideIndices_[i];
            int index2 = InsideIndices_[j];

            simstate_.getSimulationBox().calculateDistance(COM_[index1], COM_[index2], distance, distance_sq);

            LinAlg3x3::normalize(distance);
            Real rdotd = LinAlg3x3::DotProduct(distance, director_);
            Real dist = std::sqrt(distance_sq);

            if (tbin_->isInRange(rdotd) && rbin_->isInRange(dist))
            {
                int i1 = tbin_->findBin(rdotd);
                int i2 = rbin_->findBin(dist);

                histogramPerIter_[i1][i2] += 1;
            }
        }
    }
}