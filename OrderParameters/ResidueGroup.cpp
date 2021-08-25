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

    strat_->getResidueIndices();
}