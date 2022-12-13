#include "Crystal.hpp"

Crystal::Crystal(std::string fname) 
: cif_fname_(fname){
    gemmi::cif::Document doc = gemmi::cif::read_file(cif_fname_);
    gemmi::cif::Column x,y,z, atypes;
    ASSERT((doc.blocks.size() > 0), "There needs to be at least one block in the cif file, blocks are sections that begin with data_ e.g. data_tomato");
    for (gemmi::cif::Block& bl : doc.blocks){
        x = bl.find_loop("_atom_site_fract_x");
        y = bl.find_loop("_atom_site_fract_y");
        z = bl.find_loop("_atom_site_fract_z");
        atypes = bl.find_loop("_atom_site_type_symbol");

        a_ = 0.1 * (Real)gemmi::cif::as_number(*bl.find_value("_cell_length_a"));
        b_ = 0.1 * (Real)gemmi::cif::as_number(*bl.find_value("_cell_length_b"));
        c_ = 0.1 * (Real)gemmi::cif::as_number(*bl.find_value("_cell_length_c"));
        alpha_ = (Real)gemmi::cif::as_number(*bl.find_value("_cell_angle_alpha")) * Constants::PI / 180.0;
        beta_  = (Real)gemmi::cif::as_number(*bl.find_value("_cell_angle_beta")) * Constants::PI / 180.0;
        gamma_  = (Real)gemmi::cif::as_number(*bl.find_value("_cell_angle_gamma")) * Constants::PI / 180.0;

        ASSERT((x.length() == y.length()), "X length must equal to y length");
        ASSERT((y.length() == z.length()), "Y length must equal to z length");
        ASSERT((atypes.length() == z.length()), "Atom types length must equal to z length");

        int num = x.length();

        for (int i=0;i<num;i++){
            Real3 p;
            p[0] = (Real)gemmi::cif::as_number(x[i]);
            p[1] = (Real)gemmi::cif::as_number(y[i]);
            p[2] = (Real)gemmi::cif::as_number(z[i]);

            cell_positions_.push_back(p);

            unit_cell_atom_types_.push_back(atypes[i]);
        }
    }

    unique_atom_types_ = unit_cell_atom_types_;

    Algorithm::unique(unique_atom_types_);

    construct_matrix();
}

void Crystal::construct_matrix(){
    Real factor = (std::cos(alpha_) - std::cos(beta_) * std::cos(gamma_)) / std::sin(gamma_);
    m_[0][0] = a_;
    m_[0][1] = b_ * std::cos(gamma_);
    m_[0][2] = c_ * std::cos(beta_);
    m_[1][0] = 0.0;
    m_[1][1] = b_ * std::sin(gamma_);
    m_[1][2] = c_ * factor;
    m_[2][0] = 0.0;
    m_[2][1] = 0.0;
    m_[2][2] = c_ * std::sqrt(1 - std::pow(std::cos(beta_),2.0) - factor * factor);
}

void Crystal::calculate(){
    total_cell_pos_.clear();
    total_atom_types_.clear();
    total_resname_.clear();

    // start tiling
    for (int i=0;i<nx_;i++){
        for (int j=0;j<ny_;j++){
            for (int k=0;k<nz_;k++){
                Real3 vc = {{i,j,k}};
                Real3 xc = LinAlg3x3::MatrixDotVector(m_, vc);

                for (auto [ind,vec] : Algorithm::enumerate(cell_positions_)){
                    total_cell_pos_.push_back(xc + vec);
                    total_atom_types_.push_back(unit_cell_atom_types_[ind]);
                    total_resname_.push_back(unit_cell_atom_types_[ind]);
                }
            }
        }
    }
}

void Crystal::writeOutput(const std::string& fname){
    mda_tools::WriteGroFile(fname, total_cell_pos_, total_atom_types_, total_resname_, {{100,100,100}});
}