#include "grcost.h"

namespace CalculationRegistry
{
    registry_<grcost> registergrcost("grcost");
}

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

    // read information about residues 
    input.pack_.ReadString("residue", ParameterPack::KeyType::Required, residueName_);
    initializeResidueGroup(residueName_);
    auto& res = getResidueGroup(residueName_).getResidues();
    numresidues_ = res.size();
    numatoms_ = res[0].atoms_.size();
    distanceCOMIndices_.resize(numatoms_,0);
    std::iota(distanceCOMIndices_.begin(), distanceCOMIndices_.end(),1);
    input.pack_.ReadVectorNumber("distanceCOM", ParameterPack::KeyType::Optional, distanceCOMIndices_);
    distanceCOMIndices_ = distanceCOMIndices_ - 1;

    input.pack_.ReadNumber("headindex", ParameterPack::KeyType::Optional, headindex_);
    input.pack_.ReadNumber("tailindex", ParameterPack::KeyType::Optional, tailindex_);
    headindex_--;
    tailindex_--;

    numtbins_ = tbin_ -> getNumbins();
    numrbins_ = rbin_ -> getNumbins();

    histogramDotProduct_.resize(numrbins_, std::vector<Real>(numtbins_,0.0));
    histogram_.resize(numrbins_, std::vector<Real>(numtbins_, 0.0));

    // read in the array 
    arrayRead_ = input.pack_.ReadArrayNumber("array", ParameterPack::KeyType::Optional, arr_);

    registerOutputFunction("grcost", [this](std::string name) -> void {this -> printgrcost(name);});
    registerOutputFunction("histogram", [this](std::string name) -> void {this -> printhistogram(name);});
}


void grcost::calculate()
{
    // first let's calculate all the COM
    auto& res = getResidueGroup(residueName_).getResidues();
    COM_.clear();
    COM_.resize(res.size());
    uij_.clear();
    uij_.resize(res.size());
    distanceCOM_.clear();
    distanceCOM_.resize(res.size());

    Qtensor_.fill({});

    #pragma omp parallel for
    for (int i=0;i<res.size();i++)
    {
        COM_[i] = calcCOM(res[i]);
        distanceCOM_[i] = CalculationTools::getCOM(res[i], simstate_, distanceCOMIndices_);

        Real3 headpos = res[i].atoms_[headindex_].positions_;
        Real3 tailpos = res[i].atoms_[tailindex_].positions_;

        Real3 distance;
        Real distance_sq;

        simstate_.getSimulationBox().calculateDistance(headpos, tailpos, distance, distance_sq);

        LinAlg3x3::normalize(distance);

        uij_[i] = distance;
    }

    OutsideIndices_.clear();
    InsideIndices_ = InsidePVIndices(COM_, OutsideIndices_);
    int N = InsideIndices_.size();

    // first let's calculate the Qtensor inside of the probe volume 
    #pragma omp parallel
    {
        Matrix localQ={};
        #pragma omp for
        for (int i=0;i<InsideIndices_.size();i++){
            int index = InsideIndices_[i];

            Matrix singleQ = LinAlg3x3::LocalQtensor(uij_[index]);

            LinAlg3x3::matrix_accum_inplace(localQ, singleQ);
        }

        #pragma omp critical
        {
            LinAlg3x3::matrix_accum_inplace(Qtensor_, localQ);
        }
    }

    LinAlg3x3::matrix_mult_inplace(Qtensor_, 1.0/(2.0*N));
    auto result = LinAlg3x3::OrderEigenSolver(Qtensor_);

    // get the director out of the result
    if (! arrayRead_){
        for (int i=0;i<3;i++){
            director_[i] = result.second[i][0];
        }
    }
    else{director_ = arr_;}

    histogramDPPerIterBuffer_.set_master_object(histogramDotProductPerIter_);
    histogramPerIterBuffer_.set_master_object(histogramPerIter_);

    #pragma omp parallel
    {
        auto& DPPerIter   = histogramDPPerIterBuffer_.access_buffer_by_id();
        auto& histPerIter = histogramPerIterBuffer_.access_buffer_by_id();

        DPPerIter.clear();
        DPPerIter.resize(numrbins_, std::vector<Real>(numtbins_,0.0));
        histPerIter.clear();
        histPerIter.resize(numrbins_, std::vector<int>(numtbins_,0.0));

        #pragma omp for
        for (int i=0;i<N;i++){
            for (int j=i+1;j<N;j++){
                Real3 distance;
                Real distance_sq;

                int index1 = InsideIndices_[i];
                int index2 = InsideIndices_[j];

                simstate_.getSimulationBox().calculateDistance(distanceCOM_[index1], distanceCOM_[index2], distance, distance_sq);

                LinAlg3x3::normalize(distance);
                Real rdotd = LinAlg3x3::DotProduct(distance, director_);
                Real dist = std::sqrt(distance_sq);

                Real dotproduct = LinAlg3x3::DotProduct(uij_[index1], uij_[index2]);

                if (tbin_->isInRange(rdotd) && rbin_->isInRange(dist)){
                    int i1 = rbin_->findBin(dist);
                    int i2 = tbin_->findBin(rdotd);

                    histPerIter[i1][i2] += 1;
                    DPPerIter[i1][i2] += dotproduct;
                }
            }
        }

        #pragma omp for
        for (int i=0;i<N;i++){
            for (int j=0;j<OutsideIndices_.size();j++){
                Real3 distance;
                Real distance_sq;

                int index1 = InsideIndices_[i];
                int index2 = OutsideIndices_[j];

                simstate_.getSimulationBox().calculateDistance(distanceCOM_[index1], distanceCOM_[index2], distance, distance_sq);

                LinAlg3x3::normalize(distance);
                Real rdotd = LinAlg3x3::DotProduct(distance, director_);
                Real dist = std::sqrt(distance_sq);

                Real dotproduct = LinAlg3x3::DotProduct(uij_[index1], uij_[index2]);

                if (tbin_->isInRange(rdotd) && rbin_->isInRange(dist)){
                    int i1 = rbin_->findBin(dist);
                    int i2 = tbin_->findBin(rdotd);

                    histPerIter[i1][i2] += 1;
                    DPPerIter[i1][i2] += dotproduct;
                }
            }
        }

        #pragma omp critical
        {
            for (int i=0;i<numrbins_;i++){
                for (int j=0;j<numtbins_;j++){
                    histogram_[i][j] += histPerIter[i][j];
                    histogramDotProduct_[i][j] += DPPerIter[i][j];
                }
            }
        }
    }
}

void grcost::finishCalculate()
{
    int numframes = simstate_.getTotalFrames();

    for (int i=0;i<numrbins_;i++){
        for (int j=0;j<numtbins_;j++){
            if (histogram_[i][j] != 0){
                histogramDotProduct_[i][j] /= histogram_[i][j];
            }
            histogram_[i][j] /= numframes;
        }
    }
}

void grcost::printgrcost(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    ofs << "# index1 index2 value\n";

    for (int i=0;i<numrbins_;i++){
        for (int j=0;j<numtbins_;j++){
            ofs << i << " " << j << " " << histogramDotProduct_[i][j] << "\n";
        }
    }

    ofs.close();
}

void grcost::printhistogram(std::string name){
    std::ofstream ofs;
    ofs.open(name);

    ofs << "# index1 index2 value\n";
    for (int i=0;i<numrbins_;i++){
        for (int j=0;j<numtbins_;j++){
            ofs << i << " " << j << " " << histogram_[i][j] << "\n";
        }
    }
    ofs.close();
}