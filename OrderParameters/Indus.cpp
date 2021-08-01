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
}

void Indus::calculate()
{
    auto& ag  = simstate_.getAtomGroup(atomGroupName_);
    auto& pos = ag.getAtomPositions();
    auto& pv  = pv_.getProbeVolume(pvName_);
    IndusIndicesBuffer_.set_master_object(indusIndices_);

    #pragma omp parallel
    {
        Real local_N = 0;
        Real local_Ntilde = 0.0;

        auto& indicesbuffer_ = IndusIndicesBuffer_.access_buffer_by_id();
        #pragma omp for
        for (int i=0;i<pos.size();i++)
        { 
            int globalIndices = ag.AtomGroupIndices2GlobalIndices(i); 

            ProbeVolumeOutput output = pv.calculate(pos[i]);
            
            ASSERT((output.htilde_x_ <= 1 && output.htilde_x_ >= 0), "htilde_x_ is not less or equal to 1 and larger or equal to 0, it is equal to " << output.htilde_x_);

            // if htilde_x larger than 0
            if(output.htilde_x_ > 0)
            {
                 indicesbuffer_.push_back(globalIndices);
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