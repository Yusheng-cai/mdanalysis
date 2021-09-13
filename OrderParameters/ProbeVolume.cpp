#include "SimulationState.h"

ProbeVolume::ProbeVolume(ProbeVolumeInput& input):simstate_(input.simstate), simbox_(simstate_.getSimulationBox()){};