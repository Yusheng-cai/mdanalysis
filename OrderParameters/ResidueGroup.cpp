#include "ResidueGroup.h"

ResidueGroup::ResidueGroup(const ResidueInput& input)
:pack_(input.pack_), grofile_(input.grofile_), top_(input.top_)
{
    pack_.ReadVectorString("selection", ParameterPack::KeyType::Required, selection_str_);
    pack_.ReadString("name", ParameterPack::KeyType::Required, name_);
 
    // the first item of selection describes the selection method --> i.e. atom_index
    AtomGroupParsingInput parsingInput = { grofile_, selection_str_ };
    strat_ = stratptr(AtomGroupParsingRegistry::Factory::instance().create(selection_str_[0],parsingInput));
    
    // This is ensured to be sorted by AtomGroupStrategy
    strat_->Parse(AtomGroupGlobalIndices_);

    // get all the indices of the residues from parsing strategy
    const auto& residueIndices = strat_->getResidueIndices();
    atomSize_ = 0;
    Residues_.clear();
    TotalResidues_.atoms_.clear();

    // a set is always sorted in C++
    for (auto it = residueIndices.begin();it != residueIndices.end();it++){
        int resindex = *it;
        auto& r = top_.getResidueByIndex(resindex);
        Residues_.push_back(r);
        int numatoms = r.atoms_.size();
        AtomsPerResidue_.push_back(numatoms);
        atomSize_ += numatoms;

        for (int i=0;i<r.atoms_.size();i++){
            TotalResidues_.atoms_.push_back(r.atoms_[i]);
        }
    }
}

void ResidueGroup::update(const VectorReal3& total_atoms)
{
    for (int i=0;i<Residues_.size();i++){ 
        for (int j=0;j<Residues_[i].atoms_.size();j++){
            // be careful, atomNumber is 1 based
            int atNum = Residues_[i].atoms_[j].atomNumber_ - 1;
            Residues_[i].atoms_[j].positions_ = total_atoms[atNum];
        }
    }

    for (int i=0;i<TotalResidues_.atoms_.size();i++){
        // be careful, atomNumber is 1 based
        int atNum = TotalResidues_.atoms_[i].atomNumber_ - 1;
        TotalResidues_.atoms_[i].positions_ = total_atoms[atNum];
    }
}