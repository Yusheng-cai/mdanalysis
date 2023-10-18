#include "ClathrateCageOccupancy.h"
#include "SimulationState.h"
#include "Mesh.h"

namespace CalculationRegistry
{
    registry_<ClathrateCageOccupancy> registerCCO("ClathrateCageOccupancy");
}


ClathrateCageOccupancy::ClathrateCageOccupancy(const CalculationInput& input)
 : Calculation(input)
{
    pack_.ReadString("Substrate", ParameterPack::KeyType::Required, substrate_name_);
    pack_.ReadString("index_file", ParameterPack::KeyType::Required, index_file_);
    pack_.ReadNumber("radius", ParameterPack::KeyType::Required, radius_);
    pack_.ReadNumber("threshold", ParameterPack::KeyType::Required, threshold_);
    pack_.ReadArrayNumber("Ray", ParameterPack::KeyType::Optional, Ray_);
    radius_sq_ = radius_ * radius_;
    pack_.ReadNumber("sigma", ParameterPack::KeyType::Optional, sigma_);
    pack_.ReadNumber("n", ParameterPack::KeyType::Optional, n_);
    pack_.ReadArrayNumber("nL", ParameterPack::KeyType::Optional, nL_);
    pack_.ReadNumber("isoval", ParameterPack::KeyType::Optional, isoval_);
    pack_.Readbool("pbcMesh", ParameterPack::KeyType::Optional, pbcMesh_);
    LinAlg3x3::normalize(Ray_);
    DensityFieldInput inp= {simstate_, nL_, sigma_, n_, isoval_, pbcMesh_};
    density_ = densityptr(new DensityField(inp));


    // add the atom groups
    addAtomgroup(substrate_name_);

    readFile(index_file_);

    // initialize cell grid
    cell_ = cellptr(new CellGrid(simstate_,radius_));

    // resize neighbor waters
    const auto& substrate_pos = getAtomGroup(substrate_name_);
    neighbor_water_.resize(substrate_pos.getNumAtoms(), 0);
    neighbor_water_iter_.resize(substrate_pos.getNumAtoms(),0);

    registerPerIterOutputFunction("Substrate_nearby", [this](std::ofstream& ofs) -> void {this -> printAtomsNearby(ofs);});
    registerPerIterOutputFunction("Substrate_within", [this](std::ofstream& ofs) -> void {this -> printAtomsWithin(ofs);});
}

void ClathrateCageOccupancy::calculate(){
    // correct for atom group
    std::vector<Real3> total_pos = simstate_.getTotalAtomPos();
    const auto& substrate_atom = getAtomGroup(substrate_name_).getAtoms();
    const auto& substrate_pos = getAtomGroup(substrate_name_).getAtomPositions();
    int time_frame = simstate_.getFrameNumber();
    std::vector<Real3> ClathratePos;

    // get the ice positions
    for (int i=0;i<indices_[time_frame].size();i++){
        int index = indices_[time_frame][i];
        ASSERT((index < total_pos.size()), "Index out of range.");
        Real3 shiftedpos = simstate_.getSimulationBox().shiftIntoBox(total_pos[index]);
        ClathratePos.push_back(shiftedpos);
    }

    // obtain the mesh from the instantaneous interface
    Mesh mesh;
    density_->CalculateInstantaneousField(ClathratePos, mesh);

    // get the triangles, vertices of the mesh
    const auto& tri = mesh.gettriangles();
    const auto& v   = mesh.getvertices();
    Real3 boxLength = simstate_.getSimulationBox().getSides();

    #pragma omp parallel
    {
        int local_num=0;
        std::vector<int> local_index;
        #pragma omp for
        for (int i=0;i<substrate_pos.size();i++){
            int num_intersect_pos=0;
            int num_intersect_neg=0;
            // define the ice position
            Real3 O = substrate_pos[i];
            for (auto ti : tri){
                // declare the A B C of the triangle
                Real3 A,B,C;
                A = v[ti[0]].position_, B=v[ti[1]].position_, C=v[ti[2]].position_;

                // shift periodic triangle into whole
                MeshTools::ShiftPeriodicTriangle(v, ti.triangleindices_, boxLength, A, B, C);

                Real t,u,v,t2;
                Real3 neg_Ray = Ray_ * (-1.0);
                if (MeshTools::MTRayTriangleIntersection(A,B,C,O,Ray_, t,u,v)){
                    if (t > 0){
                        num_intersect_pos += 1;
                        continue;
                    }
                }

                if (MeshTools::MTRayTriangleIntersection(A,B,C,O,neg_Ray,t,u,v)){
                    if (t > 0){
                        num_intersect_neg += 1;
                    }
                }

                if ((num_intersect_neg == 1) && (num_intersect_pos==1)){
                    break;
                }
            }
            if ((num_intersect_neg == 1) && (num_intersect_pos==1)){
                local_num += 1;
                local_index.push_back(substrate_atom[i].index + 1);
            }
        }

        #pragma omp critical
        {
            num_within_volume_ += local_num;
            methane_within_.insert(methane_within_.end(), 
                                    local_index.begin(), 
                                    local_index.end());
        }
    }

    const auto& substrate = getAtomGroup(substrate_name_);
    std::vector<std::vector<int>> substrate_cell_indices = cell_->calculateIndices(substrate_pos);
    neighbor_water_iter_.clear();
    neighbor_water_iter_.resize(substrate_pos.size(),0.0);

    #pragma omp parallel
    {
        std::vector<Real> local_vec(substrate_pos.size(), 0);

        #pragma omp for
        for (int i=0;i<indices_[time_frame].size();i++){
            int water_ind = indices_[time_frame][i];
            int water_intind = cell_->getCellGridIntIndex(total_pos[water_ind]);
            std::vector<int> neighbor_index = cell_->getNeighborIndex(total_pos[water_ind]);

            for (auto& vec_neigh : substrate_cell_indices){
                for (int ind : vec_neigh){
                    Real3 dist;
                    Real dist_sq;
                    simstate_.getSimulationBox().calculateDistance(total_pos[water_ind], substrate_pos[ind], dist, dist_sq);

                    if (dist_sq <= radius_sq_){
                        local_vec[ind] += 1;
                    }
                }
            }
        }

        #pragma omp critical
        {
            neighbor_water_ = neighbor_water_ + local_vec;
            neighbor_water_iter_ = neighbor_water_iter_ + local_vec;
        }
    }

    const auto& substrate_atoms = substrate.getAtoms();
    for (int i=0;i<substrate_atoms.size();i++){
        if (neighbor_water_iter_[i] >= threshold_){
            num_with_neighbors_ += 1;
            methane_nearby_.push_back(substrate_atom[i].index + 1);
        }
    }

    std::cout << "Num with neighbors = " << num_with_neighbors_ << "\n";
    std::cout << "Num within = " << num_within_volume_ << "\n";
}


void ClathrateCageOccupancy::printAtomsNearby(std::ofstream& ofs){
    int time = simstate_.getTime();

    ofs << time << " ";
    for (int i=0;i<methane_nearby_.size();i++){
        ofs << methane_nearby_[i] << " ";
    }
    ofs << "\n";
}

void ClathrateCageOccupancy::printAtomsWithin(std::ofstream& ofs){
    int time = simstate_.getTime();

    ofs << time << " ";
    for (int i=0;i<methane_within_.size();i++){
        ofs << methane_within_[i] << " ";
    }
    ofs << "\n";
}

void ClathrateCageOccupancy::update(){
    num_within_volume_ = 0;
    num_with_neighbors_ = 0;
    cell_->update();

    methane_within_.clear();
}

void ClathrateCageOccupancy::readFile(const std::string& filename){
    std::ifstream ifs;
    ifs.open(filename);
    ASSERT((ifs.is_open()), "File " << filename << " cannot be opened.");

    std::stringstream ss;
    std::string sentence;
    while(std::getline(ifs, sentence)){
        ss.str(sentence);
        std::vector<int> local_index;
        int index;
        // get rid of the first one
        ss >> index;
        while (ss >> index){
            local_index.push_back(index-1);
        }
        indices_.push_back(local_index);
        ss.clear();
    }
    ifs.close(); 

    ASSERT((simstate_.getTotalFrames() <= indices_.size()), "The index file has less frames than the simulation.");
}

