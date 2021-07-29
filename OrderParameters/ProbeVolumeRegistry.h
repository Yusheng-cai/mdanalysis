#pragma once 
#include "ProbeVolume.h"
#include "tools/Assert.h"

#include <map>
#include <string>
#include <memory>

class ProbeVolumeRegistry
{
    public:
        using ProbeVolumePtr = std::unique_ptr<ProbeVolume>; 

        ProbeVolumeRegistry(){};
        ~ProbeVolumeRegistry(){};

        void registerProbeVolume(const std::string name, const ProbeVolume* pv_ptr);
        ProbeVolume& getProbeVolume(const std::string name) const;

    private:
        std::map<std::string, ProbeVolumePtr> MapName2PV;  
};