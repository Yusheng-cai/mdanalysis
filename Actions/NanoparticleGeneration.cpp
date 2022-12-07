#include "NanoparticleGeneration.hpp"

NanoparticleGeneration::NanoparticleGeneration(const std::vector<Real3> &particle_pos, Real De, Real alpha,
                                               Real re, Real sigma, Real epsilon)
    : particle_pos_(particle_pos), De_(De), alpha_(alpha), re_(re), sigma_(sigma), epsilon_(epsilon)
{
    Morse_ = Mptr(new MorsePotential(De_, re_, alpha_));
    LJ_ = LJptr(new LennardJones(sigma_, epsilon_));

    ProcessParticlePositions();
}

void NanoparticleGeneration::ProcessParticlePositions()
{
    ASSERT((particle_pos_.size() != 0), "Passed in empty particle positions to generate nanoparticle");
    numParticles_ = particle_pos_.size();
    std::vector<Real> xpos(numParticles_);
    std::vector<Real> ypos(numParticles_);
    std::vector<Real> zpos(numParticles_);

    for (int i = 0; i < numParticles_; i++){
        xpos[i] = particle_pos_[i][0];
        ypos[i] = particle_pos_[i][1];
        zpos[i] = particle_pos_[i][2];
    }

    Real xmin, xmax, ymin, ymax, zmin, zmax, Lx, Ly, Lz;
    xmin = Algorithm::min(xpos);
    xmax = Algorithm::max(xpos);
    Lx = xmax - xmin;
    ymin = Algorithm::min(ypos);
    ymax = Algorithm::max(ypos);
    Ly = ymax - ymin;
    zmin = Algorithm::min(zpos);
    zmax = Algorithm::max(zpos);
    Lz = zmax - zmin;

    box_min_ = {{xmin - Lx, ymin - Ly, zmin - Lz}};
    box_max_ = {{xmax + Lx, ymax + Ly, zmax + Lz}};

    for (int i = 0; i < 3; i++){
        box_N_[i] = std::ceil((box_max_[i] - box_min_[i]) / box_step_);
        box_size_[i] = box_N_[i] * box_step_;
    }

    potential_grid_.resize(box_N_, 0);

    CreateOffset();
}

void NanoparticleGeneration::CreateOffset()
{
    int off_each_pos = std::ceil(cutoff_ / box_step_);

    offsets_.clear();
    for (int i = -off_each_pos; i <= off_each_pos; i++){
        for (int j = -off_each_pos; j <= off_each_pos; j++){
            for (int k = -off_each_pos; k <= off_each_pos; k++){
                offsets_.push_back({{i, j, k}});
            }
        }
    }
}

NanoparticleGeneration::INT3 NanoparticleGeneration::FindIJKOnGrid(const Real3 &pos)
{
    Real3 diff = pos - box_min_;
    INT3 ijk;

    for (int i = 0; i < 3; i++){
        ijk[i] = std::round(diff[i] / box_step_);
    }

    return ijk;
}

NanoparticleGeneration::Real3 NanoparticleGeneration::GridPositionFromIJK(const INT3 &ijk)
{
    Real3 grid_position;

    for (int i = 0; i < 3; i++)
    {
        grid_position[i] = ijk[i] * box_step_ + box_min_[i];
    }

    return grid_position;
}

void NanoparticleGeneration::FixIndex(INT3 &ijk)
{
    for (int i = 0; i < 3; i++)
    {
        ijk[i] %= box_N_[i];
        if (ijk[i] < 0)
        {
            ijk[i] += box_N_[i];
        }
    }
}

void NanoparticleGeneration::CalculateDistance(const Real3 &v1, const Real3 &v2, Real3 &dist, Real &distsq)
{
    dist = v1 - v2;
    distsq = 0.0;
    for (int i = 0; i < 3; i++)
    {
        if (dist[i] > (box_size_[i] * 0.5)){
            dist[i] -= box_size_[i];
        }
        else if (dist[i] < (-box_size_[i] * 0.5)){
            dist[i] += box_size_[i];
        }

        distsq += (dist[i] * dist[i]);
    }
}

void NanoparticleGeneration::GeneratePotentialGrid()
{
    potential_grid_buffer_.set_master_object(potential_grid_);
    for (auto it = potential_grid_buffer_.beginworker(); it != potential_grid_buffer_.endworker(); it++){
        it->resize(box_N_, 0);
    }

    #pragma omp parallel
    {
        Lattice<Real>& potential_buffer = potential_grid_buffer_.access_buffer_by_id();

        #pragma omp for
        for (int i = 0; i < numParticles_; i++){
            INT3 ijk = FindIJKOnGrid(particle_pos_[i]);

            for (INT3 off : offsets_){
                INT3 off_ijk = ijk + off;
                FixIndex(off_ijk);

                Real3 off_pos = GridPositionFromIJK(off_ijk);

                Real3 dist;
                Real distsq;

                // calculate distance between particle and off grid
                CalculateDistance(particle_pos_[i], off_pos, dist, distsq);

                // check whether distance is within cutoff
                if (distsq < cutoff_sq_){
                    // if so calculate the energy --> add to the lattice grid
                    Real r = std::sqrt(distsq);
                    Real energy = Morse_->calculate(r);

                    potential_buffer(off_ijk) = potential_buffer(off_ijk) + energy;
                }
            }
        }
    }

    for (auto it = potential_grid_buffer_.beginworker(); it != potential_grid_buffer_.endworker(); it++){
        #pragma omp parallel for
        for (int i = 0; i < potential_grid_.getSize(); i++){
            potential_grid_[i] = potential_grid_[i] + it->operator[](i);
        }
    }
}

void NanoparticleGeneration::addSulfur(){
    INT3 min_ijk  = potential_grid_.idx_to_ijk(min_index_);
    Real3 particle_pos = GridPositionFromIJK(min_ijk);
    sulfur_positions_.push_back(particle_pos);

    for (INT3 off : offsets_){
        INT3 off_ijk = off + min_ijk;
        FixIndex(off_ijk);

        if (off_ijk == min_ijk){
            potential_grid_(min_ijk) = std::numeric_limits<Real>::max();
        }
        else{
            Real3 lattice_pos = GridPositionFromIJK(off_ijk);
            Real3 dist; 
            Real distsq;
            CalculateDistance(particle_pos, lattice_pos, dist, distsq);

            if (distsq < cutoff_sq_){
                Real energy = LJ_->calculate(std::sqrt(distsq));
                potential_grid_(off_ijk) = potential_grid_(off_ijk) + energy;
            }
        }
    }

    min_index_ = Algorithm::argmin(potential_grid_.getData());
    min_energy_ = potential_grid_[min_index_];

    // FindMinIndexNearSulfur();
}

void NanoparticleGeneration::FindMinIndexNearSulfur(){
    // curr min index
    int curr_min_idx = min_index_;
    INT3 curr_min_ijk = potential_grid_.idx_to_ijk(curr_min_idx);

    std::vector<int> indices;
    std::vector<Real> energies;

    // find the energies of the grids near the current min index 
    for (INT3 off : offsets_){
        INT3 ijk = curr_min_ijk + off;
        FixIndex(ijk);
        int idx = potential_grid_.ijk_to_idx(ijk);
        Real energy = potential_grid_(ijk);

        energies.push_back(energy);
        indices.push_back(idx);
    }

    int mini = Algorithm::argmin(energies);
    min_index_ = indices[mini];
    min_energy_ = energies[mini];
}

void NanoparticleGeneration::Generate()
{
    GeneratePotentialGrid();
    
    // find the minimum energy as well as minimum index in the first iteration
    min_index_ = Algorithm::argmin(potential_grid_.getData());
    min_energy_ = potential_grid_[min_index_];

    while (true){
        std::cout << "minimum energy = " << min_energy_ << '\n';
        if (min_energy_ > e_thresh_ ){
            break;
        }

        addSulfur();
    }

    MinimizeEnergy();
}

void NanoparticleGeneration::MinimizeEnergy(){
    Real step = 0.01;
    int iteration=0;
    Real last_total_E=0.0, this_total_E=0.0;
    Real max_force = std::numeric_limits<Real>::max();

    while (max_force > 100){
        // calculate total energy
        std::vector<Real3> sulfur_force(sulfur_positions_.size(), {0,0,0});
        this_total_E = 0.0;

        #pragma omp parallel
        {
            std::vector<Real3> local_force(sulfur_positions_.size(), {0,0,0});
            Real local_E=0.0;

            #pragma omp for
            for (int i=0;i<sulfur_positions_.size();i++){
                for (int j=0;j<particle_pos_.size();j++){
                    Real3 dist;
                    Real rsq=0.0;

                    CalculateDistance(sulfur_positions_[i], particle_pos_[j], dist, rsq);

                    if (rsq < cutoff_sq_){
                        Real r = std::sqrt(rsq);

                        Real e = Morse_->calculate(r);
                        Real3 force = Morse_->calculate_force(dist, r);

                        local_force[i] = local_force[i] + force;

                        local_E += e;
                    }
                }

                for (int j=i+1;j<sulfur_positions_.size();j++){
                    Real3 dist;
                    Real rsq=0.0;

                    CalculateDistance(sulfur_positions_[i], sulfur_positions_[j], dist, rsq);

                    if (rsq < cutoff_sq_){
                        Real r = std::sqrt(rsq);

                        Real e = LJ_->calculate(r);
                        Real3 force = LJ_->calculate_force(dist, r);

                        local_force[i] = local_force[i] + force;
                        local_force[j] = local_force[j] - force;

                        local_E += e;
                    }
                }
            }

            #pragma omp critical
            this_total_E += local_E;

            #pragma omp critical 
            for (int i=0;i<sulfur_positions_.size();i++){
                sulfur_force[i] = sulfur_force[i] + local_force[i];
            }
        }


        max_force = - std::numeric_limits<Real>::max();
        for (int i=0;i<sulfur_force.size();i++){
            for (int j=0;j<3;j++){
                if (std::abs(sulfur_force[i][j]) > max_force){
                    max_force = std::abs(sulfur_force[i][j]);
                }
            }
        }

        // if we are on the first iteration --> simply update the last energy
        if (iteration == 0){
            last_total_E = this_total_E;
        }
        else{
            if (this_total_E < last_total_E){
                step = 1.2 * step;
            }
            else{
                step = 0.2 * step;
            }

            last_total_E = this_total_E;
        }
        std::cout << "step = " << step << "\n";

        // perform the step
        for (int i=0;i<sulfur_positions_.size();i++){
            sulfur_positions_[i] = sulfur_positions_[i] + sulfur_force[i] / max_force * step;
        }

        std::cout << "At iteration " << iteration << " energy = " << this_total_E << "\n";

        iteration++;
    }
}

void NanoparticleGeneration::writeGroFile(std::string name){
    FILE* fp;
    fp = fopen(name.c_str(), "w");

    int num_AU = particle_pos_.size();
    int num_S  = sulfur_positions_.size();
    int total_atom = num_AU + num_S;

    fprintf(fp, "%s\n", "Nanoparticle");
    fprintf(fp, "\t%d\n", total_atom);

    for (int i=1;i<particle_pos_.size()+1;i++){
        fprintf(fp, gro_c_string, i, "AU", "AU", i, \
        particle_pos_[i-1][0], particle_pos_[i-1][1], particle_pos_[i-1][2]);
    }

    for (int i=1;i<sulfur_positions_.size()+1;i++){
        fprintf(fp, gro_c_string, i + num_AU, "S", "S", i + num_AU, \
        sulfur_positions_[i-1][0], sulfur_positions_[i-1][1], sulfur_positions_[i-1][2]);
    }

    fprintf(fp, "\t%d %d %d", box_size_[0], box_size_[1], box_size_[2]);

    fclose(fp);
}