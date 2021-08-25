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
    Residues_.resize(residueIndices.size());
    AtomsPerResidue_.resize(residueIndices.size());
    int index=0;

    // a set is always sorted in C++
    for (auto it = residueIndices.begin();it != residueIndices.end();it++)
    {
        Residues_[index] = GroResidues_[*it];

        for (int i=0;i<Residues_.size();i++)
        {
            Residues_[index].atoms_[i].mass_ = top_.getMassFromAtomName(Residues_[index].atoms_[i].atomName_);
        }
    }
}

void ResidueGroup::update(const VectorReal3& total_atoms_)
{
    for (int i=0;i<Residues_.size();i++)
    { 
        for (int j=0;j<Residues_[i].atoms_.size();j++)
        {
            Residues_[i].atoms_[j].positions_ = total_atoms_[Residues_[i].atoms_[j].atomNumber_];
        }
    }
}