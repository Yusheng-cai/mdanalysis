#include "CylindricalRDF.h"

namespace CalculationRegistry
{
    registry_<CylindricalRDF> registerCylindricalRDF("CylindricalRDF");
}

CylindricalRDF::CylindricalRDF(const CalculationInput& input)
: Calculation(input)
{
    // read in the 2 residues 
    pack_.ReadString("residue1", ParameterPack::KeyType::Required, residue1_);
    pack_.ReadString("residue2", ParameterPack::KeyType::Required, residue2_);

    // initialize the residue groups 
    initializeResidueGroup(residue1_, "COMIndices1", COMIndices1_, COM1_);
    initializeResidueGroup(residue2_, "COMIndices2", COMIndices2_, COM2_);
    const auto& res1 = getResidueGroup(residue1_);
    numres_ = res1.getResidues().size();

    // initialize the density
    rho_ = 0.0;

    // read in the 2 bins 
    auto rbinPack = pack_.findParamPack("rbin", ParameterPack::KeyType::Required);
    rbin_ = binptr(new Bin(*rbinPack));
    Real dr = rbin_->getStep();
    auto zbinPack = pack_.findParamPack("zbin", ParameterPack::KeyType::Required);
    zbin_ = binptr(new Bin(*zbinPack));
    Real dz = zbin_->getStep();

    // read in the array indicator z
    pack_.Readbool("usedirector", ParameterPack::KeyType::Optional, usedirector_);

    if (! usedirector_)
    {
        pack_.ReadArrayNumber("array", ParameterPack::KeyType::Optional, direction_);
        LinAlg3x3::normalize(direction_);
    }
    else
    {
        pack_.ReadNumber("headindex", ParameterPack::KeyType::Optional, headindex_);
        pack_.ReadNumber("tailindex", ParameterPack::KeyType::Optional, tailindex_);
        headindex_--;
        tailindex_--;
    }

    // resize cylindrical rdf
    numrbin_ = rbin_->getNumbins();
    numzbin_ = zbin_->getNumbins();
    cylindrical_rdf_.resize(numzbin_, std::vector<Real>(numrbin_,0.0));

    // calculate cylindrical shell volume
    cylindrical_volume_.resize(numzbin_, std::vector<Real>(numrbin_,0.0));
    for (int i=0;i<numzbin_;i++)
    {
        for (int j=0;j<numrbin_;j++)
        {
            cylindrical_volume_[i][j] = 2 * Constants::PI * rbin_->getRightLocationOfBin(j) * dr * dz;
        }
    }

    // read in the properties if we are doing any
    pack_.Readbool("Usr", ParameterPack::KeyType::Optional, calcUsr_);
    cylindrical_Usr_.resize(numzbin_, std::vector<Real>(numrbin_,0.0));
    cylindrical_Usr_attr_.resize(numzbin_, std::vector<Real>(numrbin_,0.0));
    cylindrical_Usr_repul_.resize(numzbin_, std::vector<Real>(numrbin_,0.0));
    if (calcUsr_)
    {
        ReadResidueIndices(residue1_, "UsrIndices1", UsrIndices1_);
        ReadResidueIndices(residue2_, "UsrIndices2", UsrIndices2_);
        Real sigma;
        pack_.ReadNumber("sigma", ParameterPack::KeyType::Optional, sigma);
        pack_.ReadNumber("Usrcutoff", ParameterPack::KeyType::Optional, Usr_cutoff_);
        beta_ = 1.0/sigma;
    }

    // read if the calculation is for interface, this will affect the density calculation
    pack_.Readbool("interface", ParameterPack::KeyType::Optional, interface_);

    // register the printing function
    registerOutputFunction("RDF", [this](std::string name)-> void {this -> printCylindricalRDF(name);});
    registerOutputFunction("Usr", [this](std::string name)-> void {this -> printCylindricalUsr(name);});
    registerOutputFunction("UsrAttr", [this](std::string name)-> void {this -> printCylindricalUsrAttr(name);});
}

void CylindricalRDF::calculate()
{
    // get the residues 
    const auto& res1 = getResidueGroup(residue1_).getResidues();
    const auto& res2 = getResidueGroup(residue2_).getResidues();

    // obtain the volume 
    if (interface_)
    {
        Real area=1.0;
        Real3 sides = simstate_.getSimulationBox().getSides();
        for (int i=0;i<2;i++)
        {
            area *= sides[i];
        }
        rho_ += res1.size()/area;
    }
    else
    {
        Real volume = simstate_.getSimulationBox().getVolume();
        rho_ += res1.size()/volume;
    }

    // calculate director using residue1
    if (usedirector_)
    {
        Matrix Qtensor = {};
        #pragma omp parallel
        {
            Matrix QtensorLocal = {};
            #pragma omp for
            for (int i=0;i<res1.size();i++)
            {
                Real3 headpos = res1[i].atoms_[headindex_].positions_;
                Real3 tailpos = res1[i].atoms_[tailindex_].positions_;

                Real3 uij;
                Real distsq;

                simstate_.getSimulationBox().calculateDistance(headpos, tailpos, uij, distsq);
                LinAlg3x3::normalize(uij);

                Matrix Qtensor_single = LinAlg3x3::LocalQtensor(uij);
                LinAlg3x3::matrix_accum_inplace(QtensorLocal, Qtensor_single);
            }

            #pragma omp critical
            LinAlg3x3::matrix_accum_inplace(Qtensor, QtensorLocal);
        }

        LinAlg3x3::matrix_mult_inplace(Qtensor, 1.0/(2.0*res1.size()));
        auto res = LinAlg3x3::OrderEigenSolver(Qtensor);
        for (int i=0;i<3;i++)
        {
            direction_[i] = res.second[i][0];
        }
    }

    // first get the COMs
    #pragma omp parallel
    {
        #pragma omp for
        for (int i=0;i<COM1_.size();i++)
        {
            COM1_[i] = CalculationTools::getCOM(res1[i], simstate_, COMIndices1_);
        }

        #pragma omp for
        for (int i=0;i<COM2_.size();i++)
        {
            COM2_[i] = CalculationTools::getCOM(res2[i], simstate_, COMIndices2_);
        }
    }

    // calculate the distances 
    #pragma omp parallel
    {
        std::vector<std::vector<Real>> local_cylindrical_rdf(numzbin_, std::vector<Real>(numrbin_,0.0));
        std::vector<std::vector<Real>> local_cylindrical_Usr(numzbin_, std::vector<Real>(numrbin_,0.0));
        std::vector<std::vector<Real>> local_cylindrical_Usr_attr(numzbin_, std::vector<Real>(numrbin_, 0.0));
        std::vector<std::vector<Real>> local_cylindrical_usr_repul(numzbin_, std::vector<Real>(numrbin_,0.0));
        #pragma omp for
        for (int i=0;i<COM1_.size();i++)
        {
            for (int j=0;j<COM2_.size();j++)
            {
                Real3 dist;
                Real distsq;

                simstate_.getSimulationBox().calculateDistance(COM1_[i], COM2_[j], dist, distsq);

                if (distsq != 0)
                {
                    // dot the dist with the array 
                    Real h = LinAlg3x3::DotProduct(dist, direction_);
                    Real rsq = distsq - h*h;
                    Real r = std::sqrt(rsq);

                    // find the 2 bins 
                    if (zbin_->isInRange(h) && rbin_->isInRange(r))
                    {
                        int hindex = zbin_->findBin(h);
                        int rindex = rbin_->findBin(r);

                        local_cylindrical_rdf[hindex][rindex] += 1;

                        if (calcUsr_)
                        {
                            Real dist = std::sqrt(distsq);
                            if (dist < Usr_cutoff_)
                            {
                                Real attr, repul, total;
                                CalculationTools::CalculateUsrBetweenPair(res1[i], res2[j], simstate_, dist, UsrIndices1_, UsrIndices2_, beta_, attr, repul, total);
                                local_cylindrical_Usr[hindex][rindex] += total;
                                local_cylindrical_Usr_attr[hindex][rindex] += attr;
                                local_cylindrical_usr_repul[hindex][rindex] += repul;
                            }
                        }
                    }
                }
            }
        }

        #pragma omp critical
        for (int i=0;i<numzbin_;i++)
        {
            for (int j=0;j<numrbin_;j++)
            {
                cylindrical_rdf_[i][j] += local_cylindrical_rdf[i][j];
                cylindrical_Usr_[i][j] += local_cylindrical_Usr[i][j];
                cylindrical_Usr_attr_[i][j] += local_cylindrical_Usr_attr[i][j];
                cylindrical_Usr_repul_[i][j] += local_cylindrical_usr_repul[i][j];
            }
        }
    }
}

void CylindricalRDF::finishCalculate()
{
    int totalFrames = simstate_.getTotalFrames();
    rho_ = rho_/totalFrames;

    if (calcUsr_)
    {
        for (int i=0;i<numzbin_;i++)
        {
            for (int j=0;j<numrbin_;j++)
            {
                cylindrical_Usr_[i][j] /= cylindrical_rdf_[i][j];
                cylindrical_Usr_attr_[i][j] /= cylindrical_rdf_[i][j];
                cylindrical_Usr_repul_[i][j] /= cylindrical_rdf_[i][j];
            }
        }
    }

    for (int i=0;i<numzbin_;i++)
    {
        for (int j=0;j<numrbin_;j++)
        {
            cylindrical_rdf_[i][j] = cylindrical_rdf_[i][j] / (totalFrames * numres_* rho_ * cylindrical_volume_[i][j]);
        }
    }
}

void CylindricalRDF::printCylindricalUsr(std::string name)
{
    ASSERT(calcUsr_, "Usr is not calculated.");

    std::ofstream ofs;
    ofs.open(name);
    for (int i=0;i<numzbin_;i++)
    {
        for (int j=0;j<numrbin_;j++)
        {
            ofs << cylindrical_Usr_[i][j] << " ";
        }
        ofs << "\n";
    }

    ofs.close();
}

void CylindricalRDF::printCylindricalUsrAttr(std::string name)
{
    ASSERT(calcUsr_, "Usr is not calculated.");

    std::ofstream ofs;
    ofs.open(name);
    for (int i=0;i<numzbin_;i++)
    {
        for (int j=0;j<numrbin_;j++)
        {
            ofs << cylindrical_Usr_attr_[i][j] << " ";
        }
        ofs << "\n";
    }

    ofs.close();
}

void CylindricalRDF::printCylindricalRDF(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    for (int i=0;i<numzbin_;i++)
    {
        for (int j=0;j<numrbin_;j++)
        {
            ofs << cylindrical_rdf_[i][j] << " ";
        }
        ofs << "\n";
    }

    ofs.close();
}