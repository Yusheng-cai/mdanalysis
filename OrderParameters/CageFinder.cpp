//
//                  ***GRADE_v1.00***
//  MyFunctions.hpp
//
//  Created by Farbod Mahmoudinobar and Cristiano L. Dias on 12/6/18.
//  Copyright © 2018 Farbod Mahmoudinobar and Cristiano L. Dias.
//
//  This file is part of GRADE.
//
//  GRADE is a free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  GRADE is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with GRADE.  If not, see <https://www.gnu.org/licenses/>.


#include "CageFinderAlgorithm.hpp"
#include "CageFinder.h"
#include "tools/Algorithm.h"
#include "SimulationState.h"
#include "LinAlgTools.h"
#include <map>

namespace CalculationRegistry
{
    registry_<CageFinder> registerCF("CageFinder");
}


CageFinder::CageFinder(const CalculationInput& input)
:Calculation(input)
{
    // read in the residue
    input.pack_.ReadString("solvent", ParameterPack::KeyType::Required, atomgroup_);
    input.pack_.ReadString("solute", ParameterPack::KeyType::Required, solute_group_);
    input.pack_.ReadNumber("cutoff", ParameterPack::KeyType::Optional, cut_off_);
    input.pack_.ReadNumber("plane_cutoff", ParameterPack::KeyType::Optional, plane_cutoff_);
    input.pack_.ReadNumber("internal_angle_cutoff", ParameterPack::KeyType::Optional, internal_angle_cutoff_);
    input.pack_.ReadNumber("Solute_512_cutoff", ParameterPack::KeyType::Optional,solute_512_cutoff_);
    input.pack_.ReadNumber("Solute_62512_cutoff", ParameterPack::KeyType::Optional, solute_62512_cutoff_);
    solute_512_cutoffsq_ = solute_512_cutoff_ * solute_512_cutoff_;
    solute_62512_cutoffsq_ = solute_62512_cutoff_ * solute_62512_cutoff_;
    addAtomgroup(atomgroup_);
    addAtomgroup(solute_group_);
    cut_off_sq_ = cut_off_ * cut_off_;

    cell_ = cellptr(new CellGrid(simstate_, cut_off_));

    registerPerIterOutputFunction("Cage512", [this](std::ofstream& ofs)-> void {this -> PrintCage512(ofs);});
    registerPerIterOutputFunction("Cage62512", [this](std::ofstream& ofs)-> void {this -> PrintCage62512(ofs);});
    registerPerIterOutputFunction("Occupation", [this](std::ofstream& ofs) -> void {this -> PrintPerIterOccupation(ofs);});
    registerPerIterOutputFunction("NonOccupied62512", [this](std::ofstream& ofs) -> void {this -> PrintNonOccupied62512(ofs);});
    registerPerIterOutputFunction("NonOccupied512", [this](std::ofstream& ofs) -> void {this -> PrintNonOccupied512(ofs);});
    registerPerIterOutputFunction("CageWater", [this](std::ofstream& ofs)-> void {this ->PrintCageWater(ofs);});
    registerPerIterOutputFunction("GuestCages512", [this](std::ofstream& ofs)-> void {this -> PrintGuestCages512(ofs);});
    registerPerIterOutputFunction("GuestCages62512", [this](std::ofstream& ofs) -> void {this -> PrintGuestCages62512(ofs);});
    registerPerIterOutputFunction("Guest", [this](std::ofstream& ofs) -> void {this -> PrintGuest(ofs);});
}

void CageFinder::update(){
    cell_->update();
    ring5_.clear();
    ring6_.clear();
    Cage_512_rings_.clear();
    Cage_62512_rings_.clear();
    solvent_neighbor_indices_.clear();
    occupied_512_ = 0;
    occupied_62512_ = 0;
    Non_occupied_cage62512_.clear();
    Non_occupied_cage512_.clear();
    Cage_512_atoms_.clear();
    Cage_62512_atoms_.clear();
    guest_indices_.clear();
}

void CageFinder::PrintPerIterOccupation(std::ofstream& ofs){
    int step = simstate_.getStep();

    if (step == 0){
        ofs << "# Cage512\Cage62512\tSoluteIn512\tSoluteIn62512\n";
    }

    ofs << Cage_512_rings_.size() << " " << Cage_62512_rings_.size() << " " << 
    occupied_512_ << " " << occupied_62512_ << "\n";
}

void CageFinder::PrintNonOccupied512(std::ofstream& ofs){
    const auto& ag = getAtomGroup(atomgroup_).getAtoms();
    int step = simstate_.getStep();

    ofs << step << " ";

    std::map<int, bool> map;
    for (int cage_ind : Non_occupied_cage512_){
        bool p;
        auto cage = Cage_512_rings_[cage_ind];
        for (int i=0;i<cage.size();i++){
            auto atom_indices = ring5_[cage[i]];

            for (int a : atom_indices){
                if (! Algorithm::IsInMap(map, a, p)){
                    ofs << ag[a].index+1 << " ";
                    Algorithm::InsertInMap(map, a ,true);
                }
            }
        }
    }

    ofs << "\n";

}

void CageFinder::FindUniqueCageWater(){
    // reset the vectors
    UniqueCage62512Water_.clear();
    UniqueCage512Water_.clear();
    UniqueCageWater_.clear();

    // start the calculation
    const auto& ag = getAtomGroup(atomgroup_).getAtoms();

    // iterate over the 62512 rings 
    for (int i=0;i<Cage_62512_rings_.size();i++){
        for (int j=0;j<Cage_62512_rings_[i].size();j++){
            // 6-ring
            if (j == 0 || j == 7){
                ASSERT((Cage_62512_rings_[i][j] < ring6_.size()), "Index out of range.");
                std::vector<int> atom_ids = ring6_[Cage_62512_rings_[i][j]];
                ASSERT((atom_ids.size() == 7), "Ring 6 Wrong");
                for (int k=0;k<atom_ids.size()-1;k++){
                    UniqueCageWater_.push_back(ag[atom_ids[k]].index+1);
                    UniqueCage62512Water_.push_back(ag[atom_ids[k]].index+1);
                }
            }
            // else it is a 5-ring
            else{
                std::vector<int> atom_ids = ring5_[Cage_62512_rings_[i][j]];
                ASSERT((atom_ids.size() == 6), "Ring 5 Wrong");
                for (int k=0;k<atom_ids.size()-1;k++){
                    UniqueCageWater_.push_back(ag[atom_ids[k]].index+1);
                    UniqueCage62512Water_.push_back(ag[atom_ids[k]].index+1);
                }
            }
        }
    }

    // find the unique atom indices for 62512
    Algorithm::unique(UniqueCage62512Water_);

    // do the 5-ring 
    for (int i=0;i<Cage_512_rings_.size();i++){
        for (int j=0;j<Cage_512_rings_[i].size();j++){
            std::vector<int> atom_ids = ring5_[Cage_512_rings_[i][j]];
            ASSERT((atom_ids.size() == 6), "Wrong");
            for (int k=0;k<atom_ids.size()-1;k++){
                UniqueCage512Water_.push_back(ag[atom_ids[k]].index+1);
                UniqueCageWater_.push_back(ag[atom_ids[k]].index + 1);
            }
        }
    }

    // find the unique atom indices for 512
    Algorithm::unique(UniqueCage512Water_);

    // find the unique atom indices total
    Algorithm::unique(UniqueCageWater_);
}

void CageFinder::calculate(){
    const auto& ag_pos = getAtomGroup(atomgroup_).getAtomPositions();
    neighbor_list_.clear();
    neighbor_list_.resize(ag_pos.size());

    solvent_neighbor_indices_ = cell_->calculateIndices(ag_pos);

    #pragma omp parallel for
    for (int i=0;i<ag_pos.size();i++){
        Real3 pos_curr =ag_pos[i];
        std::vector<int> neighbor_index = cell_->getNeighborIndex(pos_curr);
        std::vector<int> neighbors;
        std::vector<Real> distances;
        for (auto& vec_neigh : solvent_neighbor_indices_){
            for (int ind : vec_neigh){
                if (ind != i){
                    Real3 pos_next = ag_pos[ind];
                    Real dist_sq;
                    Real3 dist_vec;
                    simstate_.getSimulationBox().calculateDistance(pos_curr, pos_next,
                                                                dist_vec, dist_sq);
                    if (dist_sq <= cut_off_sq_){
                        neighbors.push_back(ind);
                        distances.push_back(dist_sq);
                    }
                }
            }
        }

        std::vector<int> argsort_index = Algorithm::argsort(distances);
        ASSERT((distances.size() >= 4), "Insufficient neighbors");

        for (int j=0;j<4;j++){
            neighbor_list_[i].push_back(neighbors[argsort_index[j]]);
        }
    }


    std::vector<std::vector<int>> Cup_512, Cup62512, Cage_512,  Cage_62512;
    Algorithm::find_all_cycles(neighbor_list_, 5, ring5_);
    Algorithm::find_all_cycles(neighbor_list_, 6, ring6_);

    CheckCoplanar(ring6_);
    CheckCoplanar(ring5_);
    CheckConvex(ring5_);
    CheckConvex(ring6_);
    RemoveDuplicateFaces(ring5_);
    RemoveDuplicateFaces(ring6_);

    std::vector<std::vector<int>> My_Neigh_Ring5, My_Neigh_ring6_ring5;
    std::vector<unsigned long> N_Neigh_Ring5, N_ring6_ring5_neigh;;

    CageFinderAlgorithm::find_shared_edges_ring5(ring5_.size(), 
                                                ring5_, My_Neigh_Ring5, 
                                                N_Neigh_Ring5);
    CageFinderAlgorithm::find_shared_edges_ring6_ring5(ring5_, ring6_, 
                             My_Neigh_ring6_ring5, N_ring6_ring5_neigh);

    CageFinderAlgorithm::cup_512_Finder(ring5_, ring5_.size(), 
                                        N_Neigh_Ring5, My_Neigh_Ring5, Cup_512);

    for (int i=0;i<ring6_.size();i++){
        std::vector<std::vector<int>> temp_nn_list(My_Neigh_ring6_ring5[i].size());
        for (int j=0;j<My_Neigh_ring6_ring5[i].size();j++){
            for (int k=j;k<My_Neigh_ring6_ring5[i].size();k++)
            {
                int first_ring5_index = My_Neigh_ring6_ring5[i][j];
                int second_ring5_index = My_Neigh_ring6_ring5[i][k];
                if (CageFinderAlgorithm::compare_adjacant(ring5_[first_ring5_index], 
                            ring5_[second_ring5_index], 5,5,ring6_[i]))
                            {
                                temp_nn_list[j].push_back(k);
                                temp_nn_list[k].push_back(j);
                            }
            }
        }

        std::vector<std::vector<int>> all_c;
        Algorithm::find_all_cycles(temp_nn_list, 6, all_c);
        RemoveDuplicate(all_c);

        for (auto& c : all_c){
            std::vector<int> temp_ind = {i};
            for (int j=0;j<c.size();j++){
                temp_ind.push_back(My_Neigh_ring6_ring5[i][c[j]]);
            }
            Cup62512.push_back(temp_ind);
        }
    }


    int count512 = CageFinderAlgorithm::remove_duplicates_map(Cup_512);
    int count62512 = Cup62512.size();

    CageFinderAlgorithm::cage_Finder(Cup_512, count512, My_Neigh_Ring5, Cage_512, 
                                    Cage_512_rings_, "1");
    CageFinderAlgorithm::cage_Finder(Cup62512, count62512, My_Neigh_Ring5, Cage_62512, 
                                      Cage_62512_rings_, "1");
    RemoveDuplicate(Cage_512_rings_);
    RemoveDuplicate(Cage_62512_rings_);

    // find the unique atoms in a single cage
    getUniqueAtomsInCages();

    FindOccupancy512();
    FindOccupancy62512();

    // find the unique waters in all cages
    FindUniqueCageWater();

}

void CageFinder::RemoveDuplicate(std::vector<std::vector<int>>& faces){
    std::map<std::vector<int>, int> map;
    std::vector<std::vector<int>> unique_faces;

    for (auto& f : faces){
        std::vector<int> new_f(f.size());
        std::partial_sort_copy(f.begin(), f.end(), new_f.begin(), new_f.end());
        int val;

        if (! Algorithm::IsInMap(map, new_f, val)){
            unique_faces.push_back(f);
            Algorithm::InsertInMap(map, new_f,1);
        }
    }

    faces = unique_faces;
}

void CageFinder::PrintNonOccupied62512(std::ofstream& ofs){
    const auto& ag = getAtomGroup(atomgroup_).getAtoms();
    int step = simstate_.getStep();

    ofs << step << " ";

    std::map<int, bool> map;
    for (int cage_ind : Non_occupied_cage62512_){
        bool p;
        auto big_cage = Cage_62512_rings_[cage_ind];
        for (int i=0;i<big_cage.size();i++){
            if (i == 0 || i == 7){
                auto atom_indices = ring6_[big_cage[i]];

                for (int a : atom_indices){
                    if (! Algorithm::IsInMap(map, a, p)){
                        ofs << ag[a].index + 1 << " ";
                        Algorithm::InsertInMap(map, a ,true);
                    }
                }
            }
            else{
                auto atom_indices = ring5_[big_cage[i]];

                for (int a : atom_indices){
                    if (! Algorithm::IsInMap(map, a, p)){
                        ofs << ag[a].index+1 << " ";
                        Algorithm::InsertInMap(map, a ,true);
                    }
                }
            }
        }
    }

    ofs << "\n";
}


void CageFinder::getUniqueAtomsInCages(){
    // get the unique atoms for 512 and 62512 cages 
    const auto& solvent = getAtomGroup(atomgroup_).getAtoms();
    Cage_512_atoms_.resize(Cage_512_rings_.size());
    Cage_62512_atoms_.resize(Cage_62512_rings_.size());

    // first obtain the unique atoms in a cage 512
    for (int i=0;i<Cage_512_rings_.size();i++){
        auto& cage = Cage_512_rings_[i];
        std::map<int,bool> Seen_atoms;
        bool a;
        for (auto& ring : cage){
            for (int atom_ind : ring5_[ring]){
                if (! Algorithm::IsInMap(Seen_atoms, atom_ind, a)){
                    Cage_512_atoms_[i].push_back(atom_ind);
                }
            }
        }
    }

    // obtain the unique atoms in a cage 62512
    for (int i=0;i<Cage_62512_rings_.size();i++){
        auto& cage = Cage_62512_rings_[i];
        std::map<int,bool> Seen_atoms;
        bool a;
        for (int j=0;j<cage.size();j++){
            std::vector<int> ring;
            if (j == 0 || j == 7){ring = ring6_[cage[j]];}
            else{ring = ring5_[cage[j]];}

            for (int atom_ind : ring){
                if (! Algorithm::IsInMap(Seen_atoms, atom_ind, a)){
                    Cage_62512_atoms_[i].push_back(atom_ind);
                }
            }
        }
    }
}


void CageFinder::FindOccupancy512(){
    const auto& guest  = getAtomGroup(solute_group_).getAtoms();
    const auto& solute = getAtomGroup(solute_group_).getAtomPositions();
    const auto& solv   = getAtomGroup(atomgroup_).getAtoms();
    const auto& solvent= getAtomGroup(atomgroup_).getAtomPositions();

    #pragma omp parallel
    {
        int occupied_512_local = 0;
        std::map<int,std::vector<int>> localmap;
        std::vector<int> local_nonooccupied;
        std::vector<int> local_guest;
        #pragma omp for
        for (int i=0;i<Cage_512_atoms_.size();i++){
            bool found_one_inside=false;
            for (int j=0;j<solute.size();j++){
                bool IsInCage = true;
                Real3 solute_p = solute[j];
                for (int atom_ind : Cage_512_atoms_[i]){
                    Real3 dist;
                    Real dist_sq;
                    simstate_.getSimulationBox().calculateDistance(
                        solute_p, 
                        solvent[atom_ind], dist, dist_sq
                    );

                    if (dist_sq > solute_512_cutoffsq_){
                        IsInCage = false;
                        goto endloop;
                    }
                }

                endloop:
                if (IsInCage){
                    std::vector<int> global_vec;
                    for (int ind : Cage_512_atoms_[i]){
                        global_vec.push_back(solv[ind].index);
                    }
                    local_guest.push_back(guest[j].index);
                    occupied_512_local += 1;
                    localmap.insert(std::make_pair(guest[j].index, global_vec));
                    found_one_inside = true;
                    break;
                }
            }

            if (! found_one_inside){
                local_nonooccupied.push_back(i);
            }
        }

        #pragma omp critical
        {
            occupied_512_ += occupied_512_local;
            mapGuestToCage512_.insert(localmap.begin(), localmap.end());
            Non_occupied_cage512_.insert(Non_occupied_cage512_.end(), 
                                        local_nonooccupied.begin(), 
                                        local_nonooccupied.end());
            guest_indices_.insert(guest_indices_.end(), 
                                  local_guest.begin(), 
                                  local_guest.end());
        }
    }    
}

void CageFinder::PrintGuest(std::ofstream& ofs){
    int num_step = simstate_.getStep();

    ofs << num_step << " ";

    for (int i=0;i<guest_indices_.size();i++){
        ofs << guest_indices_[i] << " ";
    }
    ofs << "\n";
}

void CageFinder::FindOccupancy62512(){
    const auto& guest  = getAtomGroup(solute_group_).getAtoms();
    const auto& solute = getAtomGroup(solute_group_).getAtomPositions();
    const auto& solv   = getAtomGroup(atomgroup_).getAtoms();
    const auto& solvent= getAtomGroup(atomgroup_).getAtomPositions();

    #pragma omp parallel
    {
        int occupied_62512_local = 0;
        std::map<int,std::vector<int>> localmap;
        std::vector<int> local_nonoccupied;
        std::vector<int> local_guest;
        #pragma omp for
        for (int i=0;i<Cage_62512_atoms_.size();i++){
            bool found_one_inside = false;
            for (int j=0;j<solute.size();j++){
                Real3 solute_p = solute[j];
                bool IsInCage = true;
                for (int atom_ind : Cage_62512_atoms_[i]){
                    Real3 dist;
                    Real dist_sq;

                    simstate_.getSimulationBox().calculateDistance(
                        solute_p, 
                        solvent[atom_ind], dist, dist_sq
                    );

                    if (dist_sq > solute_62512_cutoffsq_){
                        IsInCage = false;
                        goto endloop;
                    }
                }

                endloop:
                if(IsInCage){
                    std::vector<int> global_vec;
                    for (int ind : Cage_62512_atoms_[i]){
                        global_vec.push_back(solv[ind].index);
                    }
                    occupied_62512_local += 1;
                    local_guest.push_back(guest[j].index);
                    localmap.insert(std::make_pair(guest[j].index, global_vec));
                    found_one_inside = true;
                    break;
                }
            }

            if ( ! found_one_inside){
                local_nonoccupied.push_back(i);
            }
        }

        #pragma omp critical
        {
            occupied_62512_ += occupied_62512_local;
            mapGuestToCage62512_.insert(localmap.begin(), localmap.end());
            Non_occupied_cage62512_.insert(Non_occupied_cage62512_.end(),
                                        local_nonoccupied.begin(), 
                                        local_nonoccupied.end());
            guest_indices_.insert(guest_indices_.end(),
                                   local_guest.begin(),
                                    local_guest.end());
        }
    }
}


void CageFinder::CheckConvex(std::vector<std::vector<int>>& Faces){
    const auto& ag_pos = getAtomGroup(atomgroup_).getAtomPositions();
    std::vector<std::vector<int>> newF;
    #pragma omp parallel
    {
        std::vector<std::vector<int>> newF_local;
        #pragma omp for
        for (int i=0;i<Faces.size();i++){
            auto f = Faces[i];
            Real correct_angle = 180.0 * (f.size() - 2);

            Real angle = 0;
            for (int j=0;j<f.size();j++){
                int before = j-1;
                int after = (j+1) % f.size();
                if (before < 0){before += f.size();}
                Real3 vec3_before, vec3_after;
                Real dist_sq;
                simstate_.getSimulationBox().calculateDistance(
                    ag_pos[f[before]],
                    ag_pos[f[j]],
                    vec3_before, dist_sq
                );
                simstate_.getSimulationBox().calculateDistance(
                    ag_pos[f[after]],
                    ag_pos[f[j]],
                    vec3_after, dist_sq
                );

                LinAlg3x3::normalize(vec3_before);
                LinAlg3x3::normalize(vec3_after);

                Real dot = LinAlg3x3::DotProduct(vec3_after, vec3_before);
                angle += std::acos(dot) * 180 / Constants::PI;
            }

            if (std::abs(angle - correct_angle) < internal_angle_cutoff_){
                newF_local.push_back(f);
            }
        }

        #pragma omp critical
        {
            newF.insert(newF.end(), newF_local.begin(), newF_local.end());
        }

    }
    
    Faces = newF;
}

void CageFinder::PrintGuestCages512(std::ofstream& ofs){
    int timef = simstate_.getFrameNumber();
    ofs << "# Timeframe " << timef << "\n";
    for (auto it = mapGuestToCage512_.begin(); it!=mapGuestToCage512_.end(); it++){
        ofs << it->first << " ";
        for (int atom_ind : it->second){
            ofs << atom_ind << " ";
        }
        ofs << "\n";
    }
}

void CageFinder::PrintGuestCages62512(std::ofstream& ofs){
    int timef = simstate_.getFrameNumber();
    ofs << "# Timeframe " << timef << "\n";
    for (auto it = mapGuestToCage62512_.begin(); it !=mapGuestToCage62512_.end(); it++){
        ofs << it->first << " "; 
        for (int atom_ind : it ->second){
            ofs << atom_ind << " ";
        }
        ofs << "\n";
    }
}

void CageFinder::PrintCageWater(std::ofstream& ofs){
    ofs << 0 << " ";
    for (int i=0;i<UniqueCageWater_.size();i++){
        ofs << UniqueCageWater_[i] << " ";
    }
    ofs << "\n";
}

void CageFinder::PrintCage62512(std::ofstream& ofs){
    ofs << 0 << " ";
    for (int i=0;i<UniqueCage62512Water_.size();i++){
        ofs << UniqueCage62512Water_[i] << " ";
    }
    ofs << "\n";
}


void CageFinder::PrintCage512(std::ofstream& ofs){
    ofs << 0 << " ";
    for (int i=0;i<UniqueCage512Water_.size();i++){
        ofs << UniqueCage512Water_[i] << " ";
    }
    ofs << "\n";
}

void CageFinder::RemoveDuplicateFaces(std::vector<std::vector<int>>& faces){
    std::map<std::vector<int>, int> map;
    std::vector<std::vector<int>> unique_faces;

    for (auto& f : faces){
        std::vector<int> new_f(f.size());
        std::partial_sort_copy(f.begin(), f.end(), new_f.begin(), new_f.end());
        int val;

        if (! Algorithm::IsInMap(map, new_f, val)){
            f.push_back(f[0]);
            unique_faces.push_back(f);
            Algorithm::InsertInMap(map, new_f, 1);
        }
    }

    faces.clear();
    faces.insert(faces.end(), unique_faces.begin(), unique_faces.end());
}

void CageFinder::CheckCoplanar(std::vector<std::vector<int>>& Faces){
    const auto& ag_pos = getAtomGroup(atomgroup_).getAtomPositions();
    std::vector<std::vector<int>> newF;
    #pragma omp parallel
    {
        std::vector<std::vector<int>> newF_local;
        #pragma omp for
        for (int i=0;i<Faces.size();i++){
            std::vector<int> f = Faces[i];
            int num_sides = f.size();
            Eigen::MatrixXf pos(num_sides,3);

            for (int j=0;j<num_sides;j++){
                int ind = f[j];
                Real3 e;
                Real d;
                simstate_.getSimulationBox().calculateDistance(
                                ag_pos[ind], 
                                ag_pos[f[0]], 
                                e,d
                );
                for (int k=0;k<3;k++){
                    pos(j,k) = e[k];
                }
            }

            auto mean_pos = pos.colwise().mean();
            for (int j=0;j<num_sides;j++){
                for (int k=0;k<3;k++){
                    pos(j,k) -= mean_pos[k];
                }
            }

            auto covmat = pos.transpose() * pos;

            Eigen::EigenSolver<Eigen::MatrixXf> eigensolver;
            eigensolver.compute(covmat);
            Eigen::VectorXf eigen_values = eigensolver.eigenvalues().real();
            Eigen::MatrixXf eigen_vectors = eigensolver.eigenvectors().real();

            Real3 ev;
            for (int j=0;j<3;j++){ev[j] = eigen_values[j];}
            int min_index = Algorithm::argmin(ev);
            ASSERT((min_index < 3), "Error");
            auto evec = eigen_vectors.col(min_index);

            bool isPlanar = true;
            for (int j=0;j<num_sides;j++){
                Real dot_val = pos.row(j) * evec;
                if (std::abs(dot_val) > plane_cutoff_){
                    isPlanar = false;
                    break;
                }
            }

            if (isPlanar){newF_local.push_back(f);}
        }

        #pragma omp critical
        {
            newF.insert(newF.end(), newF_local.begin(), newF_local.end());
        }
    }
    

    // #pragma omp parallel
    // {
    //     std::vector<std::vector<int>> newF_local;
    //     #pragma omp for
    //    for (int i=0;i<Faces.size();i++){
    //         std::vector<int> f = Faces[i];
    //         int num_sides = f.size() - 1;
    //         std::vector<bool> isPlanar_vec;
    //         std::vector<Real3> p;
    //         for (int j=0;j<num_sides;j++){
    //             std::vector<Real3> normals;
    //             std::vector<Real3> edges;
    //             Real total_angle = 0;
    //             p.push_back(ag_pos[f[j]]);
    //             // f are [1,2,3,4,5,6,1]
    //             int count = 0;
    //             int next_index=(j + 1) % num_sides;
    //             while (count < num_sides-1){
    //                 Real3 e;
    //                 Real d;
    //                 simstate_.getSimulationBox().calculateDistance(
    //                     ag_pos[f[next_index]],
    //                     ag_pos[f[j]],
    //                     e,d
    //                 );

    //                 LinAlg3x3::normalize(e);
    //                 edges.push_back(e);

    //                 next_index = (next_index + 1) % num_sides;
    //                 count ++;
    //             }
    //             std::cout << "edges = " << edges.size() << "\n";

    //             for (int j=0;j<edges.size()-1;j++){
    //                 Real3 n = LinAlg3x3::CrossProduct(edges[j], edges[j+1]);
    //                 LinAlg3x3::normalize(n);
    //                 normals.push_back(n);
    //             }
    //             std::cout << "normals size = " << normals.size() << "\n";
    //             std::cout << "total_angle = " << total_angle << "\n";

    //             bool isPlanar = true;
    //             for (int j=0;j<normals.size()-1;j++){
    //                 Real dot = LinAlg3x3::DotProduct(normals[j], normals[j+1]);
    //                 Real angle = std::acos(dot) * 180 / Constants::PI;
    //                 std::cout << angle << " ";

    //                 if (! (angle < 45 || (180 - angle) < 45)){
    //                     isPlanar = false;
    //                 }
    //             }
    //             std::cout << "\n";
    //             isPlanar_vec.push_back(isPlanar);
    //         }
    //         bool total_isPlanar = isPlanar_vec[0];

    //         bool all_equal=true; 
    //         for (int j=1;j<isPlanar_vec.size();j++)
    //         {
    //             if (isPlanar_vec[j] != total_isPlanar){
    //                 all_equal = false;
    //             }
    //         }

    //         if (all_equal){
    //             for (auto a : p){
    //                 std::cout << a << "\n";
    //             }
    //             ASSERT((false), "Wrong");
    //         }

    //             // for (int k=j;j<f.size()-1;j++){
    //             //     Real3 e;
    //             //     Real d;
    //             //     simstate_.getSimulationBox().calculateDistance(ag_pos[f[j]],
    //             //                                                 ag_pos[f[0]], 
    //             //                                                 e, d);
    //             //     LinAlg3x3::normalize(e);
    //             //     edges.push_back(e);
    //             // }

    //             // for (int j=0;j<edges.size()-1;j++){
    //             //     Real3 n = LinAlg3x3::CrossProduct(edges[j], edges[j+1]);
    //             //     LinAlg3x3::normalize(n);
    //             //     normals.push_back(n);
    //             // }

    //             // bool isPlanar = true;
    //             // for (int j=0;j<normals.size()-1;j++){
    //             //     Real dot = LinAlg3x3::DotProduct(normals[j], normals[j+1]);
    //             //     Real angle = std::acos(dot) * 180 / Constants::PI;

    //             //     if (! (angle < 45 || (180 - angle) < 45)){
    //             //         isPlanar = false;
    //             //     }
    //             // }

    //         if (total_isPlanar){
    //             newF.push_back(f);
    //         }

    //         }
            

    //     #pragma omp critical
    //     {
    //         newF.insert(newF.end(), newF_local.begin(), newF_local.end());
    //     }
    // }
    Faces.clear();
    Faces.insert(Faces.end(), newF.begin(), newF.end());
}
