#include "ChillPlus.h"

namespace CalculationRegistry
{
    registry_<ChillPlus> registerChillPlus("ChillPlus");
}

ChillPlus::ChillPlus(const CalculationInput& input)
: Calculation(input)
{
    // read the solvation shell radius
    pack_.ReadNumber("solvation_shell_r", ParameterPack::KeyType::Optional, solvation_shell_r_);
    solvation_shell_r_squared_ = solvation_shell_r_ * solvation_shell_r_;

    // create the cell grid object
    cell_ = cellptr(new CellGrid(simstate_, solvation_shell_r_));

    pack_.ReadNumber("harmonics_degree", ParameterPack::KeyType::Optional, harmonics_degree_);
    num_m_ = harmonics_degree_ * 2 + 1;
    for (int i=-harmonics_degree_;i<=harmonics_degree_;i++){
        m_vec_.push_back(i);
    }
    ASSERT((m_vec_.size() == num_m_), "The size of m vector is incorrect.");

    // read atom group
    pack_.ReadString("atomgroup", ParameterPack::KeyType::Required, atomgroup_name_);
    addAtomgroup(atomgroup_name_);

    // read surface atom group
    pack_.ReadVectorString("surfaceAtomGroups", ParameterPack::KeyType::Optional, surface_atomgroups_);
    pack_.ReadNumber("surface_shell_r", ParameterPack::KeyType::Optional, surface_shell_r_);
    cell_surface_ = cellptr(new CellGrid(simstate_, surface_shell_r_));
    surface_shell_r_squared_ = surface_shell_r_ * surface_shell_r_;
    for (auto s : surface_atomgroups_){
        addAtomgroup(s);
    }

    // initialize the probevolumes
    initializeProbeVolumes();
    initializeNotInProbeVolumes();

    // register per iter output functions
    registerPerIterOutputFunction("ice_like_atoms", [this](std::ofstream& ofs) -> void {this -> printIceLikeAtoms(ofs);});
    registerPerIterOutputFunction("ice_types", [this](std::ofstream& ofs) -> void {this -> printIceTypeNum(ofs);});

    // register output file output
    registerOutputFileOutputs("num_ice_like_atoms", [this](void) -> Real {return this -> getNumIceLikeAtoms();});
    registerOutputFileOutputs("num_chillplus_atoms", [this](void) -> Real {return this -> getNumIceLikeChillPlus();});

    // whether we are finding true ice 
    pack_.Readbool("FindTrueIce", ParameterPack::KeyType::Optional, findtrueice_);
    if (findtrueice_){
        pack_.ReadNumber("sigma", ParameterPack::KeyType::Required, sigma_);
        pack_.ReadNumber("n", ParameterPack::KeyType::Required, n_);
        pack_.ReadArrayNumber("nL", ParameterPack::KeyType::Required, nL_);
        pack_.ReadArrayNumber("Ray", ParameterPack::KeyType::Required, Ray_);
        pack_.ReadNumber("isoval", ParameterPack::KeyType::Required, isoval_);
        pack_.Readbool("pbcMesh", ParameterPack::KeyType::Optional, pbcMesh_);
        cutMesh_ = pack_.ReadArrayNumber("volume", ParameterPack::KeyType::Optional, volume_);
        LinAlg3x3::normalize(Ray_);
        DensityFieldInput input = {simstate_, nL_, sigma_, n_, isoval_, pbcMesh_};
        density_ = densityptr(new DensityField(input));
    }

    // whether we are performing surface correction for number of ice like atoms 
    pack_.Readbool("SurfaceCorrection", ParameterPack::KeyType::Optional, surface_correction_);
    if (surface_correction_){
        pack_.ReadNumber("surface_cutoff", ParameterPack::KeyType::Optional, surface_cutoff_);
        pack_.ReadNumber("ice_cutoff", ParameterPack::KeyType::Optional, ice_cutoff_);
        pack_.ReadNumber("surface_threshold", ParameterPack::KeyType::Optional, surface_threshold_);
        pack_.ReadNumber("ice_threshold", ParameterPack::KeyType::Optional, ice_threshold_);
        ice_cutoff_sq_ = ice_cutoff_ * ice_cutoff_;
        surface_cutoff_sq_ = surface_cutoff_ * surface_cutoff_;

        cell_ice_corr_ = cellptr(new CellGrid(simstate_, ice_cutoff_));
        cell_surface_corr_ = cellptr(new CellGrid(simstate_, surface_cutoff_));
    }
}

void ChillPlus::calculate(){
    // obtain the atom groups
    const auto& pos = getAtomGroup(atomgroup_name_).getAtomPositions();
    const auto& ag  = getAtomGroup(atomgroup_name_);

    // for the positions --> calculate whether or not they are in the PV
    IsInsideProbeVolume_.clear(); 
    IsInsideProbeVolume_.resize(pos.size(), 0);
    for (int i=0;i<pos.size();i++){
        if (isInPV(pos[i])){
            IsInsideProbeVolume_[i] = 1;
        }
    }

    // make cell list for all the water atoms 
    water_cell_list_ = cell_->calculateIndices(pos);

    // define all the necessary arrays
    std::vector<std::vector<int>> neighbor_indices(pos.size());
    std::vector<std::vector<Real>> neighbor_distance(pos.size());
    std::vector<std::vector<Real3>> neighbor_vector_distance(pos.size());
    std::vector<std::vector<std::complex<Real>>> qlm(pos.size());
    std::vector<std::vector<Real>> cij(pos.size());

    std::fill(type_nums_.begin(), type_nums_.end(), 0);

    // define the cell indices for surface atom groups
    surface_cell_list_.clear();
    for (auto s : surface_atomgroups_){
        const auto& posA = getAtomGroup(s).getAtomPositions();
        surface_cell_list_.push_back(cell_surface_->calculateIndices(posA));
    }


    // first calculate within water itself
    #pragma omp parallel for
    for (int i=0;i<pos.size();i++){
        std::vector<int> neighbor_cell_indices = cell_->getNeighborIndex(pos[i]);
        for (int neighbor_cell_ind : neighbor_cell_indices){
            for (int neighbor_ind : water_cell_list_[neighbor_cell_ind]){
                if (neighbor_ind != i){
                    Real3 distance;
                    Real distsq;
                    simstate_.getSimulationBox().calculateDistance(pos[neighbor_ind], pos[i], distance, distsq);

                    // if the other water is within the first hydration shell
                    if (distsq <= solvation_shell_r_squared_){
                        neighbor_indices[i].push_back(neighbor_ind);
                        neighbor_distance[i].push_back(std::sqrt(distsq));
                        neighbor_vector_distance[i].push_back(distance);
                    }
                }
            }
        }
    }

    #pragma omp parallel for
    for (int i=0;i<neighbor_indices.size();i++){
        // if larger than 4, then we sort and reassign the neighbor indices with the 4 of the closest distances
        if (neighbor_distance[i].size() > 4){
            std::vector<int> argsortIndices = Algorithm::argsort(neighbor_distance[i], true);
            std::vector<int> ind(4);
            std::vector<Real> dist(4);
            std::vector<Real3> vec_dist(4); 

            // reorder the indices, vector distances and distances 
            for (int j=0;j<4;j++){
                ind[j] = neighbor_indices[i][argsortIndices[j]];
                vec_dist[j] = neighbor_vector_distance[i][argsortIndices[j]];
                dist[j] = neighbor_distance[i][argsortIndices[j]];
            }

            // update the neighbor indices , distance and the vectors 
            neighbor_indices[i] = ind;
            neighbor_distance[i] = dist;
            neighbor_vector_distance[i] = vec_dist;
        }
    }

    // calculate the qlm = 1/4 * \sum_{j=1}^{4} Ylm (rij)
    #pragma omp parallel for
    for (int i=0;i<neighbor_indices.size();i++){
        std::vector<std::complex<Real>> ql_i;
        for (int m : m_vec_){
            std::complex<Real> qlm_i = {{0,0}};
            Real factor = 1.0 / neighbor_indices[i].size();
            for (int j=0;j<neighbor_indices[i].size();j++){
                Real3 dist = neighbor_vector_distance[i][j];

                // normalized distance
                dist = dist / neighbor_distance[i][j];

                // azimuthal angle phi (0, 2 pi)
                // first quadrant
                Real phi = std::atan2(dist[1],dist[0]);
                if (phi < 0){
                    phi += 2 * Constants::PI;
                }

                // zenithal angle theta (0, pi)
                Real theta = std::acos(dist[2]);

                // calculate the spherical harmonics
                auto res = boost::math::spherical_harmonic(harmonics_degree_, m, theta, phi);
                qlm_i += res;
            }
            ql_i.push_back(qlm_i * factor);
        }
        qlm[i] = ql_i;
    }

    // once we have all the qlm, then we can calculate cij = qi * qj / |qi| |qj|
    #pragma omp parallel for 
    for (int i=0;i<pos.size();i++){
        cij[i].resize(neighbor_indices[i].size());
        for (int j=0;j<neighbor_indices[i].size();j++){
            int neighborIdx = neighbor_indices[i][j];
            std::complex<Real> top = {{0,0}};
            std::complex<Real> qli = {{0,0}};
            std::complex<Real> qlneighbor = {{0,0}};
            for (int m=0;m<m_vec_.size();m++){
                top += qlm[i][m] * std::conj(qlm[neighborIdx][m]);
                qli += qlm[i][m] * std::conj(qlm[i][m]);
                qlneighbor += qlm[neighborIdx][m] * std::conj(qlm[neighborIdx][m]);
            }
            qli = std::sqrt(qli);
            qlneighbor = std::sqrt(qlneighbor);
            cij[i][j] = top.real()/(qli.real() * qlneighbor.real());
        }
    }



    std::vector<int> ice_t(6,0);
    std::vector<INT2> bonds(cij.size(), {0,0});
    // calculate whether or not an atom is ice-like or not 
    for (int i=0;i<cij.size();i++){
        // only check if it's ice if it's inside the probe volume
        if (IsInsideProbeVolume_[i]){
            // initialize the number of staggered or eclipse bonds
            INT2 bond = {{0,0}};
            for (int j=0;j<cij[i].size();j++){
                if (isStaggered(cij[i][j])){
                    bond[BondTypes::staggered] += 1;
                }
                else if (isEclipse(cij[i][j])){
                    bond[BondTypes::eclipse] += 1;
                }
            }
            // assign the bond counts
            bonds[i] = bond;
        }
    }

    //  now identify the ice 
    //  type            E    S    neighbor
    //  hexagonal       1    3    4
    //  cubic ice       0    4    4
    //  interfacial ice any  2    4    must have at least one first neighbor with more than 2 staggered bond 
    //                  0    3    4    must have at least one first neighbor with more than 1 staggered bond
    //  clathrate       4    0    4 
    //  interfacial c   3    any  4 
    //  liquid          N/A  N/A  any
    // vector of bool keeping track of whether an atom is ice like
    is_ice_like_.clear();
    is_ice_like_.resize(pos.size(), 0);
    types_.clear();
    types_.resize(pos.size(), ChillPlusTypes::Liquid);
    ice_like_indices_.clear();
    // now we perform the chill plus algorithm for identifying ice
    for (int i=0;i<cij.size();i++){
        auto bond = bonds[i];
        
        // define type based on the number of bonds
        int type  = ChillPlusTypes::Liquid;

        // chill plus identification
        if ((neighbor_indices[i].size() == 4) && IsInsideProbeVolume_[i]){
            // if it's in the map but does not belong to interfacial, hexagonal or cubic --> water
            if (Algorithm::IsInMap(mapBondToIceType_, bond, type)){
                // check if it is interfacial
                if ((type == ChillPlusTypes::Interfacial)){
                    // count number of neighbor with 1. more than 2 staggered bonds, 2. more than 1 staggered bonds 
                    int num_staggered_2=0;
                    int num_staggered_1=0;
                    for (int j=0;j<4;j++){
                        if (bonds[neighbor_indices[i][j]][BondTypes::staggered] > 2){
                            num_staggered_2++;
                        }

                        if (bonds[neighbor_indices[i][j]][BondTypes::staggered] > 1){
                            num_staggered_1++;
                        }
                    }


                    // if itself has 2 staggered bonds, then at least a neighbor need to have more than 2 staggered bonds 
                    if (bonds[i][BondTypes::staggered]==2){
                        if (num_staggered_2>=1){
                            ice_like_indices_.push_back(i);
                            is_ice_like_[i]=1;
                        }
                        else{
                            type = ChillPlusTypes::Liquid;
                        }
                    }
                    // if itself has 3 staggered bonds, then at least a neighbor need to have more than 1 staggered bonds
                    else if (bonds[i][BondTypes::staggered]==3){
                        if (num_staggered_1>=1){
                            ice_like_indices_.push_back(i);
                            is_ice_like_[i]=1;
                        }
                        else{
                            type = ChillPlusTypes::Liquid;
                        }
                    }
                }
                // check if it is either hexagonal or cubic 
                else if ((type == ChillPlusTypes::Hexagonal) || (type == ChillPlusTypes::Cubic)){
                    // identify i to ice like indices --> local index
                    ice_like_indices_.push_back(i);
                    is_ice_like_[i] = 1;
                }
                // else it is liquid (clathrate doesn't count)
            }
        }

        types_[i]=type;
    }


    //then calculate wrt surface atom group
    #pragma omp parallel
    {
        std::vector<int> local_ind;
        #pragma omp for
        for (int i=0;i<pos.size();i++){
            // at this point, we only check for waters which have 3 neighbors and is identified as water
            if ((neighbor_indices[i].size() == 3) && (! is_ice_like_[i]) && (IsInsideProbeVolume_[i])){
                // find the neighbor cell indices 
                std::vector<int> neighbor_cell_indices = cell_surface_->getNeighborIndex(pos[i]);

                // iterate over the surface atomgroups --> check if there is an water atom within surface_shell_r 
                bool is_surface=false;
                for (int j=0;j<surface_atomgroups_.size();j++){
                    const auto& PosSurface = getAtomGroup(surface_atomgroups_[j]).getAtomPositions();
                    for (int neighbor_cell_ind : neighbor_cell_indices){
                        for (int neighbor_ind : surface_cell_list_[j][neighbor_cell_ind]){
                            Real3 distance;
                            Real distsq;
                            simstate_.getSimulationBox().calculateDistance(PosSurface[neighbor_ind], pos[i], distance, distsq);

                            // if the distance to surface is less than the threshold 
                            if (distsq <= surface_shell_r_squared_){
                                is_surface=true;

                                // go to the label EXE
                                goto EXE;
                            } 
                        }
                    }
                }


                EXE:
                if (is_surface){
                    int num_nn_3S=0, num_nn_2S_1E=0, num_nn_3E=0;
                    // iterate over all the neighbors of the current atom of interest 
                    for (int k=0;k<neighbor_indices[i].size();k++){
                        int neighbor_idx = neighbor_indices[i][k];

                        // if this neighbor does not have 4 neighbors, then we skip
                        if (neighbor_indices[neighbor_idx].size() != 4){
                            continue;
                        }

                        if (bonds[neighbor_idx][BondTypes::staggered] >= 3){
                            num_nn_3S ++;
                        }
                        else if ((bonds[neighbor_idx][BondTypes::staggered] >= 2) && (bonds[neighbor_idx][BondTypes::eclipse] == 1)){
                            num_nn_2S_1E ++;
                        }
                        else if (bonds[neighbor_idx][BondTypes::eclipse] >= 3){
                            num_nn_3E ++;
                        }
                    }

                    if ( (num_nn_3S + num_nn_2S_1E) >= 2 ) {
                        local_ind.push_back(i);
                        is_ice_like_[i] = 1;
                        types_[i] = ChillPlusTypes::Surface;
                    }
                }
            }
        }

        #pragma omp critical
        {
            if (local_ind.size() != 0){
                ice_like_indices_.insert(ice_like_indices_.end(), local_ind.begin(), local_ind.end());
            }
        }
    }

    // before correction number 
    num_ice_like_atoms_before_correction_ = ice_like_indices_.size();

    // find the water indices 
    for (int i=0;i<is_ice_like_.size();i++){
        if ((! is_ice_like_[i]) && IsInsideProbeVolume_[i]){
            water_indices_.push_back(i);
        }
    }

    // correction  
    if (findtrueice_){
        CorrectForTrueIce();
    }

    // correct ice like if necessary
    if (surface_correction_){
        CorrectIceLikeAtomsBasedOnSurface();
    }

    // assign types
    for (int i=0;i<types_.size();i++){
        if (IsInsideProbeVolume_[i]){
            type_nums_[types_[i]] += 1;
        }
    }

    // calculate number of ice like atoms
    num_ice_like_atoms_ = ice_like_indices_.size();
    std::cout << "num ice like atoms = " << num_ice_like_atoms_ << '\n';
}

void ChillPlus::ShiftTriangleWithRef(Real3& A, Real3& B, Real3& C, Real3& ref){
    Real3 centroid = (A+B+C)/3;

    Real3 shift = simstate_.getSimulationBox().calculateShift(centroid, ref);
    A = A + shift;
    B = B + shift;
    C = C + shift;
}

void ChillPlus::CorrectForTrueIce(){
    t.start();
    // correct for atom group
    const auto& ag = getAtomGroup(atomgroup_name_);
    const auto& pos= ag.getAtomPositions();

    // get the ice positions
    std::vector<Real3> IcePos(ice_like_indices_.size());
    for (int i=0;i<ice_like_indices_.size();i++){
        int index = ice_like_indices_[i];
        Real3 shiftedpos = simstate_.getSimulationBox().shiftIntoBox(pos[index]);
        IcePos[i] = shiftedpos;
    }

    // obtain the mesh from the instantaneous interface
    Mesh mesh;
    density_->CalculateInstantaneousField(IcePos, mesh);
    t.end();
    std::cout << "Duration of instaneousInterface = " << t.diff() << "\n";

    if (cutMesh_){
        MeshTools::CutMesh(mesh, volume_);
    }

    // get the triangles, vertices of the mesh
    const auto& tri = mesh.gettriangles();
    const auto& v   = mesh.getvertices();
    Real3 boxLength = simstate_.getSimulationBox().getSides();

    std::vector<int> newIceIndices;
    t.start();
    #pragma omp parallel
    {
        std::vector<int> localIceIndices;
        #pragma omp for
        for (int i=0;i<IcePos.size();i++){
            int num_intersect=0;
            // define the ice position
            Real3 O = IcePos[i];
            for (auto ti : tri){
                // declare the A B C of the triangle
                Real3 A,B,C;
                A = v[ti[0]].position_, B=v[ti[1]].position_, C=v[ti[2]].position_;

                // shift periodic triangle into whole
                MeshTools::ShiftPeriodicTriangle(v, ti.triangleindices_, boxLength, A, B, C);

                Real t,u,v;
                if (MeshTools::MTRayTriangleIntersection(A,B,C,O,Ray_, t,u,v)){
                    if (t > 0){
                        num_intersect += 1;

                        // there only needs to be 1 intersect triangle
                        break;
                    }
                }
            }
            if (num_intersect > 0){
                localIceIndices.push_back(ice_like_indices_[i]);
            }
        }

        #pragma omp critical
        {
            newIceIndices.insert(newIceIndices.end(), localIceIndices.begin(), localIceIndices.end());
        }
    }
    t.end();
    std::cout << "Moller Trumbore took " << t.diff() << "\n";

    ice_like_indices_.clear();
    ice_like_indices_.insert(ice_like_indices_.end(), newIceIndices.begin(), newIceIndices.end());
}

void ChillPlus::CorrectIceLikeAtomsBasedOnSurface(){
    // get all atom positions
    const auto& ag = getAtomGroup(atomgroup_name_);
    const auto& pos= ag.getAtomPositions();
    std::cout << "before correction ice indices = " << ice_like_indices_.size() << "\n";

    int total_atoms = ice_like_indices_.size() + water_indices_.size();

    // make cell list for all the water atoms 
    auto ice_cell_corr_list = cell_ice_corr_->calculateIndices(pos);
    std::vector<std::vector<std::vector<int>>> surface_cell_corr_list;
    for (auto s : surface_atomgroups_){
        const auto& posA = getAtomGroup(s).getAtomPositions();
        surface_cell_corr_list.push_back(cell_surface_corr_->calculateIndices(posA));
    }

    int iter=1;
    while (true){ 
        // obtain water indices  --> find how many ice like indices are around it 
        std::vector<int> ice_neighbors(water_indices_.size(),0);
        std::vector<int> surface_neighbors(water_indices_.size(),0);

        #pragma omp parallel for
        for (int i=0;i<water_indices_.size();i++){
            int index = water_indices_[i];

            std::vector<int> ice_neighbor_indices = cell_ice_corr_->getNeighborIndex(pos[index]);
            std::vector<int> surface_neighbor_indices = cell_surface_corr_->getNeighborIndex(pos[index]);

            // check all the neighbor cells 
            int ice_neighbor=0;
            int surface_neighbor=0;

            // no CELL grid --> for debug purposes
            // check the ice 
            // for (int j=0;j<pos.size();j++){
            //     if ((is_ice_like_[j]) && (j != index)){
            //         // check the distance 
            //         Real3 distance;
            //         Real distsq;
            //         simstate_.getSimulationBox().calculateDistance(pos[index], pos[j], distance, distsq);
            //         if (distsq <= ice_cutoff_sq_){
            //             ice_neighbor += 1;
            //         }
            //     }

            // }

            for (int neighbor_cell_ind : ice_neighbor_indices){
                for (int neighbor_ind : ice_cell_corr_list[neighbor_cell_ind]){
                    // only perform calculation if it's ice like
                    if (is_ice_like_[neighbor_ind]){
                        // check the distance 
                        Real3 distance;
                        Real distsq;
                        simstate_.getSimulationBox().calculateDistance(pos[index], pos[neighbor_ind], distance, distsq);
                        if (distsq <= ice_cutoff_sq_){
                            ice_neighbor += 1;
                        }
                    }
                }
            }

            // no CellGrid --> for debug purposes
            // for (int j=0;j<surface_atomgroups_.size();j++){
            //     const auto& surface_pos = getAtomGroup(surface_atomgroups_[j]).getAtomPositions();
            //     for (int k=0;k<surface_pos.size();k++){
            //         // check the distance
            //         Real3 distance;
            //         Real distsq;
            //         simstate_.getSimulationBox().calculateDistance(pos[index], surface_pos[k], distance, distsq);
            //         if (distsq <= surface_cutoff_sq_){
            //             surface_neighbor += 1;
            //         }
            //     }
            // }


            for (int j=0;j<surface_atomgroups_.size();j++){
                const auto& surface_pos = getAtomGroup(surface_atomgroups_[j]).getAtomPositions();
                for (int neighbor_cell_ind : surface_neighbor_indices){
                    for (int neighbor_ind : surface_cell_corr_list[j][neighbor_cell_ind]){
                        // check the distance
                        Real3 distance;
                        Real distsq;
                        simstate_.getSimulationBox().calculateDistance(pos[index], surface_pos[neighbor_ind], distance, distsq);
                        if (distsq <= surface_cutoff_sq_){
                            surface_neighbor += 1;
                        }
                    }
                }
            }

            surface_neighbors[i] = surface_neighbor;
            ice_neighbors[i] = ice_neighbor;
        }

        // set up the new water indices --> iterate over existing water indices as no ice should turn into water 
        std::vector<int> new_water_indices;
        int added = 0;
        for (int i=0;i<water_indices_.size();i++){
            int index = water_indices_[i];
            // if the water atom meets the criteria
            if ((surface_neighbors[i] >= surface_threshold_) && (ice_neighbors[i] >= ice_threshold_)){
                ASSERT((is_ice_like_[index] == 0), "There must be something wrong.");
                // change this atom to be ice like
                is_ice_like_[index] = 1;
                ice_like_indices_.push_back(index);
                types_[index] = ChillPlusTypes::Surface;
                added += 1;
            }
            else{
                new_water_indices.push_back(index);
            }
        }

        iter++;

        // break if no new water indices has been added
        if (water_indices_.size() == new_water_indices.size()){
            break;
        }

        water_indices_.clear(); 
        water_indices_.insert(water_indices_.end(), new_water_indices.begin(), new_water_indices.end());
    }

    ASSERT((water_indices_.size() + ice_like_indices_.size() == total_atoms), "After correction, the total number of atoms changed.");
    std::cout << "After correction = " << ice_like_indices_.size() << "\n";
}

void ChillPlus::printIceLikeAtoms(std::ofstream& ofs){
    int time = simstate_.getTime();
    ofs << time << " ";

    const auto& ag = getAtomGroup(atomgroup_name_);
    for (int i=0;i<ice_like_indices_.size();i++){
        ofs << ag.AtomGroupIndices2GlobalIndices(ice_like_indices_[i]) + 1 << " ";
    }
    ofs << "\n";
}

void ChillPlus::printIceTypeNum(std::ofstream& ofs){
    int step = simstate_.getStep();

    if (step == 0){
        ofs << "# Cubic Hexagonal Interfacial Clathrate InterfacialC Surface Liquid\n";
    }
    int time = simstate_.getTime();
    ofs << time << " ";
    for (int i=0;i<7;i++){
        ofs << type_nums_[i] << " ";
    }
    ofs << "\n";
}

void ChillPlus::finishCalculate()
{

}

void ChillPlus::update()
{
    cell_->update();
    cell_surface_->update();
    cell_ice_corr_->update();
    cell_surface_corr_->update();

    Ice_Indices_.clear();
    Ice_Indices_.resize(5);

    // clear water and ice_like_indices
    water_indices_.clear();
    ice_like_indices_.clear();
}