#include "ResidueGroup.h"

ResidueGroup::ResidueGroup(const ResidueInput& input)
:grofile_(input.grofile_), top_(input.top_)
{
    input.pack_.ReadVectorString("selection", ParameterPack::KeyType::Required, selection_str_);
    input.pack_.ReadString("name", ParameterPack::KeyType::Required, name_);
 
    // the first item of selection describes the selection method --> i.e. atom_index
    AtomGroupParsingInput parsingInput = { grofile_, selection_str_ };
    strat_ = stratptr(AtomGroupParsingRegistry::Factory::instance().create(selection_str_[0],parsingInput));
    
    // This is ensured to be sorted by AtomGroupStrategy
    strat_->Parse(AtomGroupGlobalIndices_);

    // get all the residues from Gro file
    const auto& GroResidues_ = grofile_.getResidues();

    // get all the indices of the residues from parsing strategy
    const auto& residueIndices = strat_->getResidueIndices();
    size_ = residueIndices.size();

    Residues_.resize(residueIndices.size());
    AtomsPerResidue_.resize(residueIndices.size());
    int index=0;

    atomSize_ = 0;

    // a set is always sorted in C++
    for (auto it = residueIndices.begin();it != residueIndices.end();it++)
    {
        Residues_[index] = GroResidues_[*it];


        for (int i=0;i<Residues_[index].atoms_.size();i++)
        {
            std::string residueName = Residues_[index].atoms_[i].residueName_;
            std::string atomname = Residues_[index].atoms_[i].atomName_;
            Residues_[index].atoms_[i].type_ = top_.getAtomTypeFromAtomNameResname(residueName, atomname);
            Residues_[index].atoms_[i].mass_ = top_.getMassFromAtomTypeResname(residueName, Residues_[index].atoms_[i].type_);
            Residues_[index].atoms_[i].charge_ = top_.getChargeFromAtomNameResname(residueName, Residues_[index].atoms_[i].atomName_);

            Molecule::atom a = Residues_[index].atoms_[i];
            TotalResidues_.atoms_.push_back(a);
        }

        atomSize_ += Residues_[index].atoms_.size();

        index++;
    }
}

void ResidueGroup::update(const VectorReal3& total_atoms_)
{
    for (int i=0;i<Residues_.size();i++)
    { 
        #ifdef MY_DEBUG
        std::cout << "------Residue " << i << "-------" << std::endl;
        #endif 

        for (int j=0;j<Residues_[i].atoms_.size();j++)
        {
            int atNum = Residues_[i].atoms_[j].atomNumber_ - 1;
            // be careful, atomNumber is 1 based
            Residues_[i].atoms_[j].positions_ = total_atoms_[atNum];
            #ifdef MY_DEBUG
            std::cout << "Atom " << atNum << ": " << Residues_[i].atoms_[j].positions_[0] << " " << Residues_[i].atoms_[j].positions_[1] << \
            " " << Residues_[i].atoms_[j].positions_[2] << std::endl;
            #endif
        }
    }

    for (int i=0;i<TotalResidues_.atoms_.size();i++)
    {
        int atNum = TotalResidues_.atoms_[i].atomNumber_ - 1;
        // be careful, atomNumber is 1 based
        TotalResidues_.atoms_[i].positions_ = total_atoms_[atNum];
    }
}