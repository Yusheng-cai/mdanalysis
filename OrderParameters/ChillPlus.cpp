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
    cell_ = cellptr(new CellGrid(simstate_, solvation_shell_r_, 2));

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
    for (auto s : surface_atomgroups_){
        addAtomgroup(s);
    }

    // initialize the probevolumes
    initializeProbeVolumes();
    initializeNotInProbeVolumes();

    // register output function 
    registerOutputFunction("ice_types", [this](std::string name)-> void {this -> printIcetypes(name);});

    // register per iter output functions
    registerPerIterOutputFunction("total_ice_indices", [this](std::ofstream& ofs) -> void {this -> printTotalIceIndicesPerIter(ofs);});
    registerPerIterOutputFunction("total_ice_nonClathrate", [this](std::ofstream& ofs) -> void {this -> printNonClathrateIndicesPerIter(ofs);});
    registerPerIterOutputFunction("Hex_Cubic_ice", [this](std::ofstream& ofs) -> void {this -> printHexCubicIce(ofs);});

    // register output file output
    registerOutputFileOutputs("num_ice_like_atoms", [this](void) -> Real {return this -> getNumIceLikeAtoms();});

    // whether we are finding true ice 
    pack_.Readbool("FindTrueIce", ParameterPack::KeyType::Optional, findtrueice_);
    if (findtrueice_){
        pack_.ReadNumber("sigma", ParameterPack::KeyType::Required, sigma_);
        pack_.ReadNumber("n", ParameterPack::KeyType::Required, n_);
        pack_.ReadArrayNumber("nL", ParameterPack::KeyType::Required, nL_);
        pack_.ReadArrayNumber("Ray", ParameterPack::KeyType::Required, Ray_);
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
        ice_cutoff_ = ice_cutoff_ * ice_cutoff_;
        surface_cutoff_ = surface_cutoff_ * surface_cutoff_;
    }
}

void ChillPlus::calculate()
{
    // clear types
    std::vector<int> ice_t(6,0);

    // obtain the atom groups
    const auto& ag = getAtomGroup(atomgroup_name_);
    const auto& pos= ag.getAtomPositions();

    // for the positions --> calculate whether or not they are in the PV
    IsInsideProbeVolume_.clear(); IsInsideProbeVolume_.resize(pos.size(), false);
    #pragma omp parallel
    {
        #pragma omp for
        for (int i=0;i<pos.size();i++){
            if (isInPV(pos[i])){
                IsInsideProbeVolume_[i] = true;
            }
        }
    }

    // define the arrays 
    water_cell_list_ = cell_->calculateIndices(pos);
    std::vector<std::vector<int>> neighbor_indices(pos.size());
    std::vector<std::vector<Real>> neighbor_distance(pos.size());
    std::vector<std::vector<Real3>> neighbor_vector_distance(pos.size());
    std::vector<std::vector<std::complex<Real>>> qlm(pos.size());
    std::vector<std::vector<Real>> cij(pos.size());


    // define the cell indices fo rsurface atom groups
    for (auto s : surface_atomgroups_){
        const auto& posA = getAtomGroup(s).getAtomPositions();
        surface_cell_list_.push_back(cell_->calculateIndices(posA));
    }

    #pragma omp parallel for
    for (int i=0;i<pos.size();i++){
        // first calculate within water itself
        std::vector<int> neighbor_cell_indices = cell_->getNeighborIndex(pos[i]);
        for (int neighbor_cell_ind : neighbor_cell_indices){
            for (int neighbor_ind : water_cell_list_[neighbor_cell_ind]){
                if (neighbor_ind != i){
                    Real3 distance;
                    Real distsq;
                    simstate_.getSimulationBox().calculateDistance(pos[neighbor_ind], pos[i], distance, distsq);

                    if (distsq <= solvation_shell_r_squared_){
                        neighbor_indices[i].push_back(neighbor_ind);
                        neighbor_distance[i].push_back(std::sqrt(distsq));
                        neighbor_vector_distance[i].push_back(distance);
                    }
                }
            }
        }

        // then calculate wrt surface atom group
        for (int j=0;j<surface_atomgroups_.size();j++){
            const auto& PosSurface = getAtomGroup(surface_atomgroups_[j]).getAtomPositions();
            for (int neighbor_cell_ind : neighbor_cell_indices){
                for (int neighbor_ind : surface_cell_list_[j][neighbor_cell_ind]){
                    Real3 distance;
                    Real distsq;
                    simstate_.getSimulationBox().calculateDistance(PosSurface[neighbor_ind], pos[i], distance, distsq);

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
        if (neighbor_indices[i].size() > 4){
            std::vector<int> argsortIndices = Algorithm::argsort(neighbor_distance[i], true);
            std::vector<int> ind(4);
            std::vector<Real> dist(4);
            std::vector<Real3> vec_dist(4); 
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

    #pragma omp parallel for 
    for (int i=0;i<neighbor_indices.size();i++){
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


    // vector of bool keeping track of whether an atom is ice like
    is_ice_like_.clear();
    is_ice_like_.resize(pos.size(), false);

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

                if (isEclipse(cij[i][j])){
                    bond[BondTypes::eclipse] += 1;
                }
            }

            // define type based on the number of bonds
            int type;
            // if there are less than 4 neighbors, then we call it water immediately
            if (neighbor_indices[i].size() < 4){
                ice_t[ChillPlusTypes::Liquid] += 1;
            }
            else{
                if (Algorithm::IsInMap(mapBondToIceType_, bond, type)){
                    if ((type == ChillPlusTypes::Clathrate) || (type == ChillPlusTypes::Interfacial) || (type == ChillPlusTypes::Hexagonal)){
                        ice_like_indices_.push_back(i);
                        is_ice_like_[i] = true;
                    }
                    ice_t[type] += 1;
                    Ice_Indices_[type].push_back(ag.AtomGroupIndices2GlobalIndices(i));
                }
                else{
                    ice_t[ChillPlusTypes::Liquid] += 1;
                    water_indices_.push_back(i);
                }
            }
        }
    }

    // obtain ice types --> number 
    ice_types_.push_back(ice_t);

    // calculate the number of ice like atoms 
    num_ice_like_atoms_ = ice_t[ChillPlusTypes::Cubic] + ice_t[ChillPlusTypes::Hexagonal] + ice_t[ChillPlusTypes::Interfacial];

    // correct ice like if necessary
    if (surface_correction_){
        CorrectIceLikeAtomsBasedOnSurface();
    }

    if (findtrueice_){
        CorrectForTrueIce();
    }
}

void ChillPlus::ShiftTriangleWithRef(Real3& A, Real3& B, Real3& C, Real3& ref){
    Real3 centroid = (A+B+C)/3;

    Real3 shift = simstate_.getSimulationBox().calculateShift(centroid, ref);
    A = A + shift;
    B = B + shift;
    C = C + shift;
}

void ChillPlus::CorrectForTrueIce(){
    // correct for atom group
    const auto& ag = getAtomGroup(atomgroup_name_);
    const auto& pos= ag.getAtomPositions();

    // get the ice positions
    std::vector<Real3> IcePos;
    for (int i=0;i<ice_like_indices_.size();i++){
        int index = ice_like_indices_[i];
        IcePos.push_back(pos[index]);
    }

    // obtain the mesh from the instantaneous interface
    Mesh mesh;
    density_->CalculateInstantaneousField(IcePos, mesh);

    // get the triangles, vertices of the mesh
    const auto& tri = mesh.gettriangles();
    const auto& v   = mesh.getvertices();
    Real3 boxLength = simstate_.getSimulationBox().getSides();

    std::vector<int> newIceIndices;
    #pragma omp parallel
    {
        std::vector<int> localIceIndices;
        #pragma omp for
        for (int i=0;i<IcePos.size();i++){
            int num_intersect=0;
            for (auto ti : tri){
                // shift the ice position into box
                Real3 O = simstate_.getSimulationBox().shiftIntoBox(IcePos[i]);

                // declare the A B C of the triangle
                Real3 A,B,C;
                A = v[ti[0]].position_, B=v[ti[1]].position_, C=v[ti[2]].position_;

                // shift periodic triangle into whole
                MeshTools::ShiftPeriodicTriangle(v, ti.triangleindices_, boxLength, A, B, C);

                // shift triangle with respect to ref
                ShiftTriangleWithRef(A,B,C,O);

                Real t,u,v;
                if (MeshTools::MTRayTriangleIntersection(A,B,C,O,Ray_, t,u,v)){
                    num_intersect += 1;
                }
            }
            if (num_intersect != 0){
                localIceIndices.push_back(ice_like_indices_[i]);
            }
        }

        #pragma omp critical
        {
            newIceIndices.insert(newIceIndices.end(), localIceIndices.begin(), localIceIndices.end());
        }
    }

    ice_like_indices_ = newIceIndices;
}

void ChillPlus::CorrectIceLikeAtomsBasedOnSurface(){
    // get all atom positions
    const auto& ag = getAtomGroup(atomgroup_name_);
    const auto& pos= ag.getAtomPositions();

    while (true){ 
        // obtain water indices  --> find how many ice like indices are around it 
        std::vector<int> ice_neighbors(water_indices_.size(),0);
        std::vector<int> surface_neighbors(water_indices_.size(),0);

        #pragma omp parallel for
        for (int i=0;i<water_indices_.size();i++){
            int index = water_indices_[i];

            std::vector<int> cell_neighbor_indices = cell_->getNeighborIndex(pos[index]);

            // check all the neighbor cells 
            int ice_neighbor=0;
            int surface_neighbor=0;

            // check the ice 
            for (int neighbor_cell_ind : cell_neighbor_indices){
                for (int neighbor_ind : water_cell_list_[neighbor_cell_ind]){
                    if (is_ice_like_[neighbor_ind]){
                        // check the distance 
                        Real3 distance;
                        Real distsq;
                        simstate_.getSimulationBox().calculateDistance(pos[index], pos[neighbor_ind], distance, distsq);
                        if (distsq < ice_cutoff_){
                            ice_neighbor += 1;
                        }
                    }
                }
            }

            for (int j=0;j<surface_atomgroups_.size();j++){
                const auto& surface_pos = getAtomGroup(surface_atomgroups_[j]).getAtomPositions();
                for (int neighbor_cell_ind : cell_neighbor_indices){
                    for (int neighbor_ind : surface_cell_list_[j][neighbor_cell_ind]){
                        // check the distance
                        Real3 distance;
                        Real distsq;
                        simstate_.getSimulationBox().calculateDistance(pos[index], surface_pos[neighbor_ind], distance, distsq);
                        if (distsq < surface_cutoff_){
                            surface_neighbor += 1;
                        }
                    }
                }
            }

            surface_neighbors[i] = surface_neighbor;
            ice_neighbors[i] = ice_neighbor;
        }

        // set up the new water indices 
        std::vector<int> new_water_indices;
        for (int i=0;i<water_indices_.size();i++){
            int index = water_indices_[i];
            if ((surface_neighbors[i] > surface_threshold_) && (ice_neighbors[i] > ice_threshold_)){
                is_ice_like_[index] = true;
                ice_like_indices_.push_back(index);
            }
            else{
                new_water_indices.push_back(index);
            }
        }

        if (water_indices_.size() == new_water_indices.size()){
            break;
        }

        water_indices_.clear(); 
        water_indices_.insert(water_indices_.end(), new_water_indices.begin(), new_water_indices.end());

    }
}

void ChillPlus::printTotalIceIndicesPerIter(std::ofstream& ofs){
    int time = simstate_.getTime();
    ofs << time << " ";
    for (int i=0;i<Ice_Indices_.size();i++){
        for (int j=0;j<Ice_Indices_[i].size();j++){
            ofs << Ice_Indices_[i][j] << " ";
        }
    }
    ofs << "\n";
}

void ChillPlus::printNonClathrateIndicesPerIter(std::ofstream& ofs){
    int time = simstate_.getTime();
    ofs << time << " ";
    for (int i=0;i<Ice_Indices_.size();i++){
        if ((i != ChillPlusTypes::Interfacial_Clathrate) && (i != ChillPlusTypes::Clathrate)){
            for (int j=0;j<Ice_Indices_[i].size();j++){
                ofs << Ice_Indices_[i][j] << " ";
            }
        }
    }
    ofs << "\n";
}

void ChillPlus::finishCalculate()
{

}

void ChillPlus::update()
{
    cell_->update();

    Ice_Indices_.clear();
    Ice_Indices_.resize(5);

    // clear water and ice_like_indices
    water_indices_.clear();
    ice_like_indices_.clear();
}

void ChillPlus::printIcetypes(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    ofs << "# Cubic Hexagonal Interfacial Clathrate Interfacial_Clathrate Liquid\n";

    for (int i=0;i<ice_types_.size();i++)
    {
        for (int j=0;j<ice_types_[i].size();j++)
        {
            ofs << ice_types_[i][j] << " ";
        }
        ofs << "\n";
    }

    ofs.close();
}

void ChillPlus::printHexCubicIce(std::ofstream& ofs){
    int time = simstate_.getTime();
    ofs << time << " ";
    for (int i=0;i<Ice_Indices_.size();i++){
        if ((i==ChillPlusTypes::Hexagonal) || (i==ChillPlusTypes::Cubic)){
            for (int j=0;j<Ice_Indices_[i].size();j++){
                ofs << Ice_Indices_[i][j] << " ";
            }
        }
    }
    ofs << "\n";
}