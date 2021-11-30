#include "gcost.h"

namespace CalculationRegistry
{
    registry_<gcost> registergcost("gcost");
}

gcost::gcost(const CalculationInput& input)
: Calculation(input)
{
    initializeProbeVolumes();
    initializeNotInProbeVolumes();

    auto binPack = input.pack_.findParamPack("bin", ParameterPack::KeyType::Required);
    bin_ = binptr(new Bin(*binPack));
    numbins_ = bin_ ->getNumbins();

    input.pack_.ReadNumber("headindex", ParameterPack::KeyType::Optional,headindex_);
    input.pack_.ReadNumber("tailindex", ParameterPack::KeyType::Optional, tailindex_);
    headindex_--;
    tailindex_--;

    input.pack_.ReadString("residue", ParameterPack::KeyType::Required, residueName_);
    initializeResidueGroup(residueName_);

    // costheta bins is always assumed to go from -1 to 1
    input.pack_.ReadNumber("numtbins", ParameterPack::KeyType::Optional, numtbins_);
    tbin_ = binptr(new Bin(trange_, numtbins_));

    registerOutputFunction("histogram", [this](std::string name) -> void {this->printHistogram(name);});

    input.pack_.ReadNumber("index", ParameterPack::KeyType::Optional, index_);

    histogramDotProduct_.resize(numbins_,0.0);

    registerCalcFunc(1, [this](Real3& ui, Real3& uj) -> Real {return this -> calcg1(ui,uj);});
    registerCalcFunc(2, [this](Real3& ui, Real3& uj) -> Real {return this -> calcg2(ui,uj);});

    // initialize histogram 2d which contains (R, t)
    histogramDotProduct2d_.resize(numbins_, std::vector<Real>(numtbins_,0.0));
}

gcost::Real gcost::calcFactor(Real3& ui, Real3& uj)
{
    auto it = MapIndexToFcn_.find(index_);

    ASSERT((it != MapIndexToFcn_.end()), "The index " << index_ << " is not registered.");

    Real val = it->second.operator()(ui,uj);

    return val;
}

gcost::Real gcost::calcg1(Real3& ui, Real3& uj)
{
    return LinAlg3x3::DotProduct(ui,uj);
}

gcost::Real gcost::calcg2(Real3& ui, Real3& uj)
{
    Real dotProduct = LinAlg3x3::DotProduct(ui, uj);
    return 1.5 * dotProduct * dotProduct  - 0.5;
}

void gcost::registerCalcFunc(int index, fcn function)
{
    auto it = MapIndexToFcn_.find(index);

    ASSERT((it == MapIndexToFcn_.end()), "The index " << index << " is already registered.");

    MapIndexToFcn_.insert(std::make_pair(index, function));
}

void gcost::calculate()
{
    const auto& res = getResidueGroup(residueName_).getResidues(); 
    COM_.clear();
    uij_.clear();
    COM_.resize(res.size());
    uij_.resize(res.size());
    histogramDotProductPerIter_.clear();
    histogramDotProductPerIter_.resize(numbins_,0.0);
    histogramPerIter_.clear();
    histogramPerIter_.resize(numbins_,0);

    // find the COM 
    #pragma omp parallel for
    for (int i=0;i<COM_.size();i++)
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

    InsideIndicesBuffer_.set_master_object(InsideIndices_);

    // first find which set of COM are inside of PV
    #pragma omp parallel
    {
        auto& buffer_ = InsideIndicesBuffer_.access_buffer_by_id();
        buffer_.clear();

        #pragma omp for
        for (int i=0;i<COM_.size();i++)
        {
            bool inPV = false;
            for (auto pv : NotInprobevolumes_)
            {
                auto out = pv -> calculate(COM_[i]);

                if (out.hx_ == 1)
                {
                    // break breaks out of the closest enclosing for loop
                    inPV = true;
                    break;
                }
            }

            if (! inPV)
            {
                for (auto pv : probevolumes_)
                {
                    auto out = pv -> calculate(COM_[i]);
                    if (out.hx_ == 1)
                    {
                        buffer_.push_back(i);
                        break;
                    }
                }
            }
        }
    }


    int size = InsideIndices_.size();
    for (auto it = InsideIndicesBuffer_.beginworker(); it != InsideIndicesBuffer_.endworker(); it ++)
    {
        size += it -> size();
    }

    InsideIndices_.reserve(size);

    for (auto it = InsideIndicesBuffer_.beginworker(); it != InsideIndicesBuffer_.endworker(); it ++)
    {
        InsideIndices_.insert(InsideIndices_.end(), it -> begin(), it -> end());
    }


    neighborDistance_.clear();
    dotProduct_.clear();
    neighborDistance_.resize(InsideIndices_.size(), std::vector<Real>(InsideIndices_.size(),0.0));
    dotProduct_.resize(InsideIndices_.size(), std::vector<Real>(InsideIndices_.size(),0.0));

    // find the distance between pair of COM 
    #pragma omp parallel for
    for (int i=0;i<InsideIndices_.size();i++)
    {
        for (int j=0;j<InsideIndices_.size();j++)
        {
            // find the distance between ith and jth COM 
            Real3 distance;
            Real val;

            int index1 = InsideIndices_[i];
            int index2 = InsideIndices_[j];
            simstate_.getSimulationBox().calculateDistance(COM_[index1], COM_[index2], distance, val);

            Real3 u1 = uij_[index1];
            Real3 u2 = uij_[index2];
            Real dotproduct = calcFactor(u1, u2);

            neighborDistance_[i][j] = std::sqrt(val);
            dotProduct_[i][j] = dotproduct;
        }
    }

    // for a single residue, let's bin it 
    histogramDotProductPerIterbuffer_.set_master_object(histogramDotProductPerIter_);
    histogramPerIterbuffer_.set_master_object(histogramPerIter_);
    histogramDotProduct2dbuffer_.set_master_object(histogramDotProduct2d_);

    #pragma omp parallel
    {
        auto& dotbuffer = histogramDotProductPerIterbuffer_.access_buffer_by_id();
        auto& histbuffer = histogramPerIterbuffer_.access_buffer_by_id();
        auto& hist2dbuffer = histogramDotProduct2dbuffer_.access_buffer_by_id();

        dotbuffer.clear();
        dotbuffer.resize(numbins_,0);

        histbuffer.clear();
        histbuffer.resize(numbins_,0);

        hist2dbuffer.clear();
        hist2dbuffer.resize(numbins_, std::vector<Real>(numtbins_,0.0));

        #pragma omp for
        for (int i=0;i<InsideIndices_.size();i++)
        {
            for (int j=i+1;j<InsideIndices_.size();j++)
            {
                Real dist = neighborDistance_[i][j];
                Real dotp = dotProduct_[i][j];

                if (bin_ ->isInRange(dist))
                {
                    // find the bins
                    int binnum = bin_ -> findBin(dist);
                    int tbinum = tbin_ -> findBin(dotp);

                    // add it to the histograms
                    histbuffer[binnum] += 1;
                    dotbuffer[binnum] += dotProduct_[i][j];
                    hist2dbuffer[binnum][tbinum] += 1;
                }
            }
        }

        #pragma omp critical 
        {
            for (int i=0;i<numbins_;i++)
            {
                histogramDotProductPerIter_[i] += dotbuffer[i];
                histogramPerIter_[i] += histbuffer[i];

                for (int j=0;j<numtbins_;j++)
                {
                    histogramDotProduct2d_[i][j] += hist2dbuffer[i][j];
                }
            }
        }
    }

    // divide histogram dot product by histogram
    for (int i=0;i<histogramDotProductPerIter_.size();i++)
    {
        if (histogramPerIter_[i] != 0)
        {
            histogramDotProductPerIter_[i] /= histogramPerIter_[i];
        }
    }

    // add histogramdotproductperiter to histogramdotproduct
    for (int i=0;i<histogramDotProduct_.size();i++)
    {
        histogramDotProduct_[i] += histogramDotProductPerIter_[i];
    }
}

void gcost::printHistogram(std::string name)
{
    std::ofstream ofs_;
    ofs_.open(name);

    for (int i=0;i<histogramDotProduct_.size();i++)
    {
        ofs_ << bin_->getLeftLocationOfBin(i) << " " <<  histogramDotProduct_[i] << "\n";
    }

    ofs_.close();
}

void gcost::finishCalculate()
{
    int numframes = simstate_.getTotalFrames();
    std::cout << "num frames = " << numframes << std::endl;

    for (int i=0;i<histogramDotProduct_.size();i++)
    {
        histogramDotProduct_[i] /= numframes;
    }

    // take care of histogram2d --> we can normalize over columns (costheta)
    for (int i=0;i<numbins_;i++)
    {
        Real sum=0.0;
        for (int j=0;j<numtbins_;j++)
        {
            sum += histogramDotProduct2d_[i][j];
        }

        for (int j=0;j<numtbins_;j++)
        {
            histogramDotProduct2d_[i][j] /= sum;
        }
    }
}