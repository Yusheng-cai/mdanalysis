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

    input.pack_.ReadNumber("headindex1", ParameterPack::KeyType::Optional,headindex1_);
    input.pack_.ReadNumber("tailindex1", ParameterPack::KeyType::Optional, tailindex1_);
    input.pack_.ReadNumber("headindex2", ParameterPack::KeyType::Optional, headindex2_);
    input.pack_.ReadNumber("tailindex2", ParameterPack::KeyType::Optional, tailindex2_);
    headindex1_--;
    tailindex1_--;

    // read 2 residues 
    input.pack_.ReadString("residue1", ParameterPack::KeyType::Required, residueName_);
    initializeResidueGroup(residueName_, COMIndices1_, "COMIndices1", COM1_);
    input.pack_.ReadString("residue2", ParameterPack::KeyType::Required, residueName2_);
    initializeResidueGroup(residueName2_, COMIndices2_, "COMIndices2", COM2_);
    auto& res1 = getResidueGroup(residueName_);
    numresidues1_ = res1.getResidues().size();
    numatoms1_ = res1[0].atoms_.size();
    auto& res2 = getResidueGroup(residueName2_);
    numresidues2_ = res2.getResidues().size();
    numatoms2_ = res2[0].atoms_.size();

    // initialize the distance COM
    initializeDistanceCOM();

    // costheta bins is always assumed to go from -1 to 1
    input.pack_.ReadNumber("numtbins", ParameterPack::KeyType::Optional, numtbins_);
    tbin_ = binptr(new Bin(trange_, numtbins_));

    // read in the index for calculation
    input.pack_.ReadNumber("index", ParameterPack::KeyType::Optional, index_);

    histogramDotProduct_.resize(numbins_,0.0);

    registerCalcFunc(1, [this](Real3& ui, Real3& uj) -> Real {return this -> calcg1(ui,uj);});
    registerCalcFunc(2, [this](Real3& ui, Real3& uj) -> Real {return this -> calcg2(ui,uj);});

    // initialize histogram 2d which contains (R, t)
    histogramDotProduct2d_.resize(numbins_, std::vector<Real>(numtbins_,0.0));

    registerOutputFunction("histogram", [this](std::string name) -> void {this->printHistogram(name);});
    registerOutputFunction("histogram2d", [this](std::string name) -> void {this ->printHistogram2d(name);});
}

void gcost::initializeDistanceCOM()
{
    distanceCOMIndices1_.resize(numatoms1_,0);
    std::iota(distanceCOMIndices1_.begin(), distanceCOMIndices1_.end(), 1);
    pack_.ReadVectorNumber("distanceCOM1", ParameterPack::KeyType::Optional, distanceCOMIndices1_);
    for (int i=0;i<distanceCOMIndices1_.size();i++)
    {
        distanceCOMIndices1_[i] -= 1;
    }

    distanceCOMIndices2_.resize(numatoms2_,0);
    std::iota(distanceCOMIndices2_.begin(), distanceCOMIndices2_.end(), 1);
    pack_.ReadVectorNumber("distanceCOM2", ParameterPack::KeyType::Optional, distanceCOMIndices1_);
    for (int i=0;i<distanceCOMIndices2_.size();i++)
    {
        distanceCOMIndices2_[i] -= 1;
    }
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
    const auto& res1 = getResidueGroup(residueName_).getResidues(); 
    const auto& res2 = getResidueGroup(residueName2_).getResidues();
    COM1_.clear();
    uij1_.clear();
    COM2_.clear();
    uij2_.clear();
    distanceCOM1_.clear();
    distanceCOM1_.resize(res1.size());
    COM2_.resize(res2.size());
    uij2_.resize(res2.size());

    histogramDotProductPerIter_.clear();
    histogramDotProductPerIter_.resize(numbins_,0.0);
    histogramPerIter_.clear();
    histogramPerIter_.resize(numbins_,0);
    histogramDotProduct2dPerIter_.clear();
    histogramDotProduct2dPerIter_.resize(numbins_, std::vector<Real>(numtbins_,0.0));

    // find the COM of 1 & 2
    #pragma omp parallel
    {
        #pragma omp for
        for (int i=0;i<numresidues1_;i++)
        {
            COM1_[i] = calcCOM(res1[i]); 
            distanceCOM1_[i] = CalculationTools::getCOM(res1[i], simstate_, distanceCOMIndices1_);

            Real3 headpos = res1[i].atoms_[headindex1_].positions_;
            Real3 tailpos = res1[i].atoms_[tailindex1_].positions_;
            Real3 distance;
            Real distance_sq;
            simstate_.getSimulationBox().calculateDistance(headpos, tailpos, distance, distance_sq);

            LinAlg3x3::normalize(distance);

            uij1_[i] = distance;
        }

        #pragma omp for
        for (int i=0;<numresidues2_;i++)
        {
            COM2_[i] = calcCOM(res2[i]);
            distanceCOM2_[i] = CalculationTools::getCOM(res2[i], simstate_, distanceCOMIndices2_);

            Real3 headpos = res2[i].atoms_[headindex2_].positions_;
            Real3 tailpos = res2[i].atoms_[tailindex2_].positions_;
            Real3 distance;
            Real distance_sq;
            simstate_.getSimulationBox().calculateDistance(headpos, tailpos, distance, distance_sq);

            LinAlg3x3::normalize(distance);

            uij2_[i] = distance;
        }
    }

    InsideIndices1_ = InsidePVIndices(COM1_, OutsideIndices1_);
    InsideIndices2_ = InsidePVIndices(COM2_, OutsideIndices2_);
    int size1 = InsideIndices1_.size();
    int size2 = InsideIndices2_.size();

    neighborDistance_.clear();
    dotProduct_.clear();
    neighborDistance_.resize(size1, std::vector<Real>(size2,0.0));
    dotProduct_.resize(size1, std::vector<Real>(size2,0.0));

    // find the distance between pair of COM 
    #pragma omp parallel for
    for (int i=0;i<size;i++)
    {
        for (int j=i+1;j<size;j++)
        {
            // find the distance between ith and jth COM 
            Real3 distance;
            Real val;

            int index1 = InsideIndices_[i];
            int index2 = InsideIndices_[j];
            simstate_.getSimulationBox().calculateDistance(distanceCOM_[index1], distanceCOM_[index2], distance, val);

            Real3 u1 = uij_[index1];
            Real3 u2 = uij_[index2];
            Real dotproduct = calcFactor(u1, u2);

            neighborDistance_[i][j] = std::sqrt(val);
            dotProduct_[i][j] = dotproduct;
        }
    }

    // calculate the distribution b/t inside the pv and outside pv
    neighborOutsideDistance_.clear();
    dotProductOutside_.clear();
    neighborOutsideDistance_.resize(size, std::vector<Real>(OutsideIndices_.size(), 0.0));
    dotProductOutside_.resize(size, std::vector<Real>(OutsideIndices_.size(),0.0));

    ASSERT((OutsideIndices_.size() + InsideIndices_.size() == res.size()), "The number of residues inside pv + number of residues \
    outside pv does not add up to " << res.size() << " and is equal to " << OutsideIndices_.size() + InsideIndices_.size());

    #pragma omp parallel for
    for (int i=0;i<size;i++)
    {
        for (int j=0;j<OutsideIndices_.size();j++)
        {
            Real3 distance;
            Real val;

            int index1 = InsideIndices_[i];
            int index2 = OutsideIndices_[j];

            simstate_.getSimulationBox().calculateDistance(distanceCOM_[index1], distanceCOM_[index2], distance, val);

            neighborOutsideDistance_[i][j] = std::sqrt(val);

            Real3 u1 = uij_[index1];
            Real3 u2 = uij_[index2];

            dotProductOutside_[i][j] = calcFactor(u1,u2);
        }
    }

    // for a single residue, let's bin it 
    histogramDotProductPerIterbuffer_.set_master_object(histogramDotProductPerIter_);
    histogramPerIterbuffer_.set_master_object(histogramPerIter_);
    histogramDotProduct2dbuffer_.set_master_object(histogramDotProduct2dPerIter_);

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
        for (int i=0;i<size;i++)
        {
            for (int j=i+1;j<size;j++)
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

        #pragma omp for 
        for (int i=0;i<size;i++)
        {
            for (int j=0;j<OutsideIndices_.size();j++)
            {
                Real dist = neighborOutsideDistance_[i][j];
                Real cos = dotProductOutside_[i][j];

                if (bin_ -> isInRange(dist))
                {
                    int binnum = bin_ -> findBin(dist);
                    int tbinnum = tbin_ -> findBin(cos);

                    histbuffer[binnum] += 1;
                    dotbuffer[binnum] += dotProductOutside_[i][j];

                    hist2dbuffer[binnum][tbinnum] += 1;
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

    // divide histogram dot product by histogram --> finding average costheta  --> <cos(theta)>
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

void gcost::printHistogram2d(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    for (int i=0;i<numbins_;i++)
    {
        for (int j=0;j<numtbins_;j++)
        {
            ofs << i << "\t" << j << "\t" << bin_ -> getCenterLocationOfBin(i) << "\t" << tbin_ -> getCenterLocationOfBin(j) \
            << "\t" << histogramDotProduct2d_[i][j] << "\n";
        }
    }

    ofs.close();
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

    // take care of histogram2d --> we can normalize over columns (costheta) , no need to time average 
    for (int i=0;i<numbins_;i++)
    {
        for (int j=0;j<numtbins_;j++)
        {
            histogramDotProduct2d_[i][j] /= numframes;
        }
    }
}