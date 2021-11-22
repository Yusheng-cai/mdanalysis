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

    registerOutputFunction("histogram", [this](std::string name) -> void {this->printHistogram(name);});

    histogramDotProduct_.resize(numbins_,0.0);
}

void gcost::calculate()
{
    const auto& res = getResidueGroup(residueName_).getResidues(); 
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
        COM_[i] = CalculationTools::getCOM(res[i], simstate_, COMIndices_);

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

    for (auto it = InsideIndicesBuffer_.beginworker(); it != InsideIndicesBuffer_.endworker(); it ++)
    {
        InsideIndices_.insert(InsideIndices_.end(), it -> begin(), it -> end());
    }

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
            Real dotproduct = LinAlg3x3::DotProduct(u1,u2);

            neighborDistance_[i][j] = val;
            dotProduct_[i][j] = dotproduct;
        }
    }

    // for a single residue, let's bin it 
    for (int i=0;i<InsideIndices_.size();i++)
    {
        for (int j=i+1;j<InsideIndices_.size();j++)
        {
            Real dist = neighborDistance_[i][j];

            if (bin_ ->isInRange(dist))
            {
                int binnum = bin_ -> findBin(dist);
                histogramDotProductPerIter_[binnum] += dotProduct_[i][j];
                histogramPerIter_[binnum] += 1;
            }
        }
    }

    // divide histogram dot product by histogram
    for (int i=0;i<histogramDotProductPerIter_.size();i++)
    {
        histogramDotProductPerIter_[i] /= histogramPerIter_[i];
    }

    for (int i=0;i<histogramDotProductPerIter_.size();i++)
    {
        std::cout << histogramDotProductPerIter_[i] << std::endl;
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
        ofs_ << histogramDotProduct_[i] << "\n";
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
}