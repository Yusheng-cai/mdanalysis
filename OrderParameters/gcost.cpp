#include "gcost.h"

namespace CalculationRegistry
{
    registry_<gcost> registergcost("gcost");
}

gcost::gcost(const CalculationInput& input)
: Calculation(input)
{
    // initialize the probe volumes
    initializeProbeVolumes();
    initializeNotInProbeVolumes();

    // check for probe volume initialization
    ASSERT((probevolumes_.size() > 0 || NotInprobevolumes_.size() > 0), "No probe volume passed in");

    // intialize the bins 
    auto binPack = input.pack_.findParamPack("bin", ParameterPack::KeyType::Required);
    bin_ = binptr(new Bin(*binPack));
    numbins_ = bin_ ->getNumbins();

    // read inputs
    pack_.ReadNumber("headindex1", ParameterPack::KeyType::Optional,headindex1_);
    pack_.ReadNumber("tailindex1", ParameterPack::KeyType::Optional, tailindex1_);
    pack_.ReadNumber("headindex2", ParameterPack::KeyType::Optional, headindex2_);
    pack_.ReadNumber("tailindex2", ParameterPack::KeyType::Optional, tailindex2_);
    pack_.Readbool("selfinteraction", ParameterPack::KeyType::Optional, selfinteraction_);
    headindex1_--;
    tailindex1_--;
    headindex2_--;
    tailindex2_--;

    // read 2 residues 
    input.pack_.ReadString("residue1", ParameterPack::KeyType::Required, residueName1_);
    initializeResidueGroup(residueName1_, "COMIndices1", COMIndices1_, COM1_);
    input.pack_.ReadString("residue2", ParameterPack::KeyType::Required, residueName2_);
    initializeResidueGroup(residueName2_, "COMIndices2", COMIndices2_, COM2_);
    auto& res1 = getResidueGroup(residueName1_);
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
    input.pack_.ReadString("name", ParameterPack::KeyType::Optional, CalcFuncName_);

    // register the calculation function --> either g1 or g2 
    registerCalcFunc("g1", [this](Real3& ui, Real3& uj) -> Real {return this -> Calculateg1(ui,uj);});
    registerCalcFunc("g2", [this](Real3& ui, Real3& uj) -> Real {return this -> Calculateg2(ui,uj);});

    // initialize histogram 2d which contains (R, t)
    histogramDotProduct_.resize(numbins_,0.0);
    histogramDotProduct2d_.resize(numbins_, std::vector<Real>(numtbins_,0.0));
    histogram_.resize(numbins_,0.0);
    histogramDotProduct2drdf_.resize(numbins_, std::vector<Real>(numtbins_,0.0));

    registerOutputFunction("histogram", [this](std::string name) -> void {this->printHistogram(name);});
    registerOutputFunction("jointdistribution", [this](std::string name) -> void {this ->printHistogram2d(name);});
    registerOutputFunction("rdfhist", [this](std::string name) -> void {this -> printrdfhist2d(name);});
    registerOutputFunction("numneighbors", [this](std::string name) -> void {this->printnumneighbors(name);});

    volume_.resize(bin_->getNumbins(),0.0);
    for (int i=0;i<bin_->getNumbins();i++)
    {
        volume_[i] = 4 * Constants::PI * std::pow(bin_ ->getLeftLocationOfBin(i),2) * bin_ ->getStep();
    }
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
    pack_.ReadVectorNumber("distanceCOM2", ParameterPack::KeyType::Optional, distanceCOMIndices2_);
    for (int i=0;i<distanceCOMIndices2_.size();i++)
    {
        distanceCOMIndices2_[i] -= 1;
    }
}

void gcost::printnumneighbors(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    for (int i=0;i<numbins_;i++)
    {
        for (int j=0;j<numtbins_;j++)
        {
            ofs << i << "\t" << j << "\t" << bin_->getLeftLocationOfBin(i) << "\t" << bin_->getLeftLocationOfBin(j) << "\t" << \
            histogramDotProduct2drdf_[i][j] << "\n";
        }
    }

    ofs.close();
}

void gcost::printrdfhist2d(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);
    std::vector<std::vector<Real>> rdf2d(numbins_, std::vector<Real>(numtbins_,0.0));

    for (int i=1;i<numbins_;i++)
    {
        for (int j=0;j<numtbins_;j++)
        {
            rdf2d[i][j] = histogramDotProduct2drdf_[i][j]/volume_[i];
        }
    }

    for (int i=0;i<numbins_;i++)
    {
        for (int j=0;j<numtbins_;j++)
        {
            ofs << i << "\t" << j << "\t" << bin_ -> getCenterLocationOfBin(i) << "\t" << tbin_ -> getCenterLocationOfBin(j) << "\t" \
            << rdf2d[i][j] << "\n";
        }
    }
    ofs.close();
}

gcost::Real gcost::CalculateFactor(Real3& ui, Real3& uj)
{
    auto it = MapNameToCalcFunc_.find(CalcFuncName_);

    ASSERT((it != MapNameToCalcFunc_.end()), "The name " << CalcFuncName_ << " is not registered.");

    Real val = it->second.operator()(ui,uj);

    return val;
}

gcost::Real gcost::Calculateg1(Real3& ui, Real3& uj)
{
    return LinAlg3x3::DotProduct(ui,uj);
}

gcost::Real gcost::Calculateg2(Real3& ui, Real3& uj)
{
    Real dotProduct = LinAlg3x3::DotProduct(ui, uj);
    return 1.5 * dotProduct * dotProduct  - 0.5;
}

void gcost::registerCalcFunc(std::string name, fcn function)
{
    auto it = MapNameToCalcFunc_.find(name);

    ASSERT((it == MapNameToCalcFunc_.end()), "The function name " << name << " is already registered.");

    MapNameToCalcFunc_.insert(std::make_pair(name, function));
}

void gcost::calculate()
{
    const auto& res1 = getResidueGroup(residueName1_).getResidues(); 
    const auto& res2 = getResidueGroup(residueName2_).getResidues();
    uij1_.clear();
    uij2_.clear();
    distanceCOM1_.clear();
    distanceCOM1_.resize(res1.size());
    distanceCOM2_.clear();
    distanceCOM2_.resize(res2.size());
    uij1_.resize(res1.size());
    uij2_.resize(res2.size());

    histogramDotProduct2dPerIter_.clear();
    histogramDotProduct2dPerIter_.resize(numbins_, std::vector<Real>(numtbins_,0.0));

    // find the COM of 1 & 2
    #pragma omp parallel
    {
        #pragma omp for
        for (int i=0;i<numresidues1_;i++)
        {
            COM1_[i] = CalculationTools::getCOM(res1[i], simstate_, COMIndices1_); 
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
        for (int i=0;i<numresidues2_;i++)
        {
            COM2_[i] = CalculationTools::getCOM(res2[i], simstate_, COMIndices2_);
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

    int insideNum = 0;
    // find the distance between pair of COM 
    #pragma omp parallel
    {
        // local histogram for theta bins
        std::vector<Real> histogramdotproductLocal(numbins_,0.0);
        std::vector<Real> histogramLocal(numbins_,0.0);
        std::vector<std::vector<Real>> hist2dLocal(numbins_, std::vector<Real>(numtbins_,0.0));
        int insideNumLocal=0;

        // find distances inside PV
        #pragma omp for
        for (int i=0;i<numresidues1_;i++)
        {
            if (isInPV(COM1_[i]))
            {
                insideNumLocal += 1;
                for (int j=0;j<numresidues2_;j++)
                {
                    // find the distance between ith and jth COM 
                    Real3 distance;
                    Real val;

                    if (selfinteraction_)
                    {
                        if (i ==j)
                        {
                            continue;
                        }
                    }

                    simstate_.getSimulationBox().calculateDistance(distanceCOM1_[i], distanceCOM2_[j], distance, val);

                    Real3 u1 = uij1_[i];
                    Real3 u2 = uij2_[j];
                    Real dotproduct = CalculateFactor(u1, u2);

                    val = std::sqrt(val);

                    // 0 distance doesn't make sense --> probably means self interaction
                    if (bin_ -> isInRange(val))
                    {
                        int binnum = bin_ -> findBin(val);
                        histogramdotproductLocal[binnum] += dotproduct;
                        histogramLocal[binnum] += 1;

                        if (tbin_ -> isInRange(dotproduct))
                        {
                            int tbinnum = tbin_ -> findBin(dotproduct);
                            hist2dLocal[binnum][tbinnum] += 1;
                        }
                    }
                }
            }
        }

        #pragma omp critical
        for (int i=0;i<numbins_;i++)
        {
            histogram_[i] += histogramLocal[i];
            histogramDotProduct_[i] += histogramdotproductLocal[i];
            for (int j=0;j<numtbins_;j++)
            {
                histogramDotProduct2dPerIter_[i][j] += hist2dLocal[i][j];
                histogramDotProduct2d_[i][j] += hist2dLocal[i][j];
            }
        }

        #pragma omp critical
        insideNum += insideNumLocal;
    }

    // add to the hist2drdf
    for (int i=0;i<numbins_;i++)
    {
        for (int j=0;j<numtbins_;j++)
        {
            histogramDotProduct2drdf_[i][j] += histogramDotProduct2dPerIter_[i][j] / insideNum;
        }
    }
}

void gcost::printHistogram2d(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    Real dcos = tbin_ -> getStep();
    Real dr   = bin_ -> getStep();
    Real sum  = 0.0;

    for (int i=0;i<numbins_;i++)
    {
        for (int j=0;j<numtbins_;j++)
        {
            sum += histogramDotProduct2d_[i][j] * dcos * dr;
        }
    }

    std::vector<std::vector<Real>> temp(numbins_, std::vector<Real>(numtbins_,0));
    for (int i=0;i<numbins_;i++)
    {
        for (int j=0;j<numtbins_;j++)
        {
            temp[i][j] = histogramDotProduct2d_[i][j] / sum;
        }
    }

    for (int i=0;i<numbins_;i++)
    {
        for (int j=0;j<numtbins_;j++)
        {
            ofs << i << "\t" << j << "\t" << bin_ -> getCenterLocationOfBin(i) << "\t" << tbin_ -> getCenterLocationOfBin(j) \
            << "\t" << temp[i][j] << "\n";
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

    // Total cosine theta divided by total histogram
    for (int i=0;i<histogramDotProduct_.size();i++)
    {
        if (histogram_[i] != 0)
        {
            histogramDotProduct_[i] /= histogram_[i];
        }
    }

    // take care of histogram2d --> we can normalize over columns (costheta) , no need to time average 
    for (int i=0;i<numbins_;i++)
    {
        for (int j=0;j<numtbins_;j++)
        {
            histogramDotProduct2d_[i][j] /= numframes;

            // time averaged number of neighbors 
            histogramDotProduct2drdf_[i][j] /= numframes;
        }
    }
}