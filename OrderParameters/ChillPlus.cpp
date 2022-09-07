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
    for (int i=-harmonics_degree_;i<=harmonics_degree_;i++)
    {
        m_vec_.push_back(i);
    }
    ASSERT((m_vec_.size() == num_m_), "The size of m vector is incorrect.");

    // read atom group
    pack_.ReadString("atomgroup", ParameterPack::KeyType::Required, atomgroup_name_);
    addAtomgroup(atomgroup_name_);

    // register output function 
    registerOutputFunction("ice_types", [this](std::string name)-> void {this -> printIcetypes(name);});
}

void ChillPlus::calculate()
{
    // clear types
    std::vector<int> ice_t(5,0);

    // obtain the atom groups
    const auto& ag = getAtomGroup(atomgroup_name_);
    const auto& pos= ag.getAtomPositions();

    // define the arrays 
    std::vector<std::vector<int>> cell_indices = cell_->calculateIndices(pos);
    std::vector<std::vector<int>> neighbor_indices(pos.size());
    std::vector<std::vector<Real>> neighbor_distance(pos.size());
    std::vector<std::vector<Real3>> neighbor_vector_distance(pos.size());
    std::vector<std::vector<std::complex<Real>>> qlm(pos.size());
    std::vector<std::vector<Real>> cij(pos.size());

    #pragma omp parallel for
    for (int i=0;i<pos.size();i++)
    {
        std::vector<int> neighbor_cell_indices = cell_->getNeighborIndex(pos[i]);
        for (int neighbor_cell_ind : neighbor_cell_indices)
        {
            for (int neighbor_ind : cell_indices[neighbor_cell_ind])
            {
                if (neighbor_ind != i)
                {
                    Real3 distance;
                    Real distsq;
                    simstate_.getSimulationBox().calculateDistance(pos[neighbor_ind], pos[i], distance, distsq);

                    if (distsq < solvation_shell_r_squared_)
                    {
                        neighbor_indices[i].push_back(neighbor_ind);
                        neighbor_distance[i].push_back(std::sqrt(distsq));
                        neighbor_vector_distance[i].push_back(distance);
                    }
                }
            }
        }
    }

    #pragma omp parallel for
    for (int i=0;i<neighbor_indices.size();i++)
    {
        // if larger than 4, then we sort and reassign the neighbor indices with the 4 of the closest distances
        if (neighbor_indices[i].size() > 4)
        {
            std::vector<int> argsortIndices = Algorithm::argsort(neighbor_distance[i], true);
            std::vector<int> ind(4);
            std::vector<Real> dist(4);
            std::vector<Real3> vec_dist(4); 
            for (int j=0;j<4;j++)
            {
                ind[j] = neighbor_indices[i][argsortIndices[j]];
                vec_dist[j] = neighbor_vector_distance[i][argsortIndices[j]];
                dist[j] = neighbor_distance[i][argsortIndices[j]];
            }
            neighbor_indices[i] = ind;
            neighbor_distance[i] = dist;
            neighbor_vector_distance[i] = vec_dist;
        }
    }

    #pragma omp parallel for
    for (int i=0;i<neighbor_indices.size();i++)
    {
        std::vector<std::complex<Real>> ql_i;
        for (int m : m_vec_)
        {
            std::complex<Real> qlm_i = {{0,0}};
            Real factor = 1.0 / neighbor_indices[i].size();
            for (int j=0;j<neighbor_indices[i].size();j++)
            {
                Real3 dist = neighbor_vector_distance[i][j];

                // normalized distance
                dist = dist / neighbor_distance[i][j];

                // azimuthal angle phi (0, 2 pi)
                // first quadrant
                Real phi = std::atan2(dist[1],dist[0]);
                if (phi < 0)
                {
                    phi += 2 * Constants::PI;
                }

                // zenithal angle theta (0, pi)
                Real theta = std::acos(dist[2]);

                auto res = boost::math::spherical_harmonic(harmonics_degree_, m, theta, phi);
                qlm_i += res;
            }
            ql_i.push_back(qlm_i * factor);
        }
        qlm[i] = ql_i;
    }

    #pragma omp parallel for 
    for (int i=0;i<neighbor_indices.size();i++)
    {
        cij[i].resize(neighbor_indices[i].size());
        for (int j=0;j<neighbor_indices[i].size();j++)
        {
            int neighborIdx = neighbor_indices[i][j];
            std::complex<Real> top = {{0,0}};
            std::complex<Real> qli = {{0,0}};
            std::complex<Real> qlneighbor = {{0,0}};
            for (int m=0;m<m_vec_.size();m++)
            {
                top += qlm[i][m] * std::conj(qlm[neighborIdx][m]);
                qli += qlm[i][m] * std::conj(qlm[i][m]);
                qlneighbor += qlm[neighborIdx][m] * std::conj(qlm[neighborIdx][m]);
            }
            qli = std::sqrt(qli);
            qlneighbor = std::sqrt(qlneighbor);
            cij[i][j] = top.real()/(qli.real() * qlneighbor.real());
        }
    }

    for (int i=0;i<cij.size();i++)
    {
        INT2 bond = {{0,0}};
        for (int j=0;j<cij[i].size();j++)
        {
            if (isStaggered(cij[i][j]))
            {
                bond[BondTypes::staggered] += 1;
            }

            if (isEclipse(cij[i][j]))
            {
                bond[BondTypes::eclipse] += 1;
            }
        }
        int type;

        if (Algorithm::IsInMap(mapBondToIceType_, bond, type))
        {
            ice_t[type] += 1;
        }
        else
        {
            ice_t[ChillPlusTypes::Liquid] += 1;
        }
    }
    ice_types_.push_back(ice_t);
}

void ChillPlus::finishCalculate()
{

}

void ChillPlus::update()
{
    cell_->update();
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