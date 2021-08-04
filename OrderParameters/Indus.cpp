#include "Indus.h"

namespace OrderParametersRegistry
{
    registry_<Indus> register_Indus("indus");
}

Indus::Indus(const OrderParametersInput& input)
:OrderParameters(input), pv_(const_cast<ProbeVolumeRegistry&>(input.pv_registry_))
{
    input.pack_.ReadString("probevolume", ParameterPack::KeyType::Required, pvName_);
    input.pack_.ReadString("atomgroup", ParameterPack::KeyType::Required, atomGroupName_);
    input.pack_.ReadString("name", ParameterPack::KeyType::Required, name_);

    registerOutput("n", [this](void)->Real {return this->getN();});
    registerOutput("ntilde", [this](void)->Real {return this->getNtilde();});

    addAtomGroup(atomGroupName_);
}

void Indus::calculate()
{
    auto& ag  = getAtomGroup(atomGroupName_);
    auto& atoms = ag.getAtoms();
    auto& pv  = pv_.getProbeVolume(pvName_);
    IndusIndicesBuffer_.set_master_object(indusIndices_);

    for (int i=0;i<derivativeOutputs_.size();i++)
    {
        // clear the current derivatives
        derivativeOutputs_[i].clear();

        // set the master object
        derivativeOutputs_[i].setMasterObject();
    }

    auto& derivativesSet = accessDerivatives(atomGroupName_);

    #pragma omp parallel
    {
        Real local_N = 0;
        Real local_Ntilde = 0.0;

        auto& indicesbuffer_ = IndusIndicesBuffer_.access_buffer_by_id();
        #pragma omp for
        for (int i=0;i<atoms.size();i++)
        { 
            ProbeVolumeOutput output = pv.calculate(atoms[i].position);
            
            ASSERT((output.htilde_x_ <= 1 && output.htilde_x_ >= 0), "htilde_x_ is not less or equal to 1 and larger or equal to 0, it is equal to " << output.htilde_x_);

            // if htilde_x larger than 0
            if(output.htilde_x_ > 0)
            {
                 indicesbuffer_.push_back(atoms[i].index);

                 derivativesSet.insertOMP(atoms[i].index, output.dhtilde_dx_);
            }

            local_Ntilde += output.htilde_x_;
            local_N      += output.hx_;
        }

        #pragma omp critical
        {
            Ntilde_ += local_Ntilde;
            N_ += local_N;
        }
    }

    derivativesSet.CombineAndClearOMPBuffer();

    // for (int i=0;i<derivativesSet.size();i++)
    // {
    //     auto& d = derivativesSet.getAtomDerivativeByIndex(i);
    //     auto deriv = d.derivatives;
    //     auto index = d.index;
    //     int local_index = ag.GlobalIndices2AtomGroupIndices(d.index);
        
    //     std::cout << "index = " << d.index << std::endl;
    //     auto atom = ag.getAtomByIndex(local_index);
    //     ASSERT((atom.index == index), "Bad.");

    //     std::cout << "Positions is : ";
    //     for (int i=0;i<3;i++)
    //     {
    //         std::cout << atom.position[i] << " ";
    //     }
    //     std::cout << "\n";

    //     std::cout << "Derivative is : ";
    //     for (int j=0;j<3;j++)
    //     {
    //         std::cout << deriv[j] << " ";
    //     }
    //     std::cout << "\n";
    // }


    int size = 0; 
    for (auto it = IndusIndicesBuffer_.beginworker(); it != IndusIndicesBuffer_.endworker();it++)
    {
        size += it->size();
    }

    indusIndices_.reserve(size);

    for (auto it = IndusIndicesBuffer_.beginworker();it != IndusIndicesBuffer_.endworker();it++)
    {
        indusIndices_.insert(indusIndices_.end(), it->begin(), it->end());
    }
}

void Indus::update()
{
    // set orderParameters to zero
    N_ = 0;
    Ntilde_ = 0;
    indusIndices_.clear();
    IndusIndicesBuffer_.clearBuffer();

    // update the probeVolume as needed 
    auto& ProbeV = pv_.accessProbeVolume(pvName_);
    ProbeV.update();
}