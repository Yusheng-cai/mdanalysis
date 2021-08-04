#pragma once
#include "xdr/GroFile.h"
#include "tools/GenericFactory.h"

#include <vector>
#include <string>
#include <map>
#include <algorithm>

struct AtomGroupParsingInput
{
    GroFile& grofile_;
    std::vector<std::string> selection_str;
};

class AtomGroupParsingStrategy
{
    public:
        AtomGroupParsingStrategy(AtomGroupParsingInput& input):grofile_(input.grofile_), selection_str_(input.selection_str){};

        virtual void Parse(std::vector<int>& indices) = 0;
    protected:
        GroFile& grofile_;
        std::vector<std::string>& selection_str_;
        std::vector<std::string> index_str_;
};

class AtomIndexParsing:public AtomGroupParsingStrategy
{
    public:
        AtomIndexParsing(AtomGroupParsingInput& input):AtomGroupParsingStrategy(input){};

        virtual void Parse(std::vector<int>& indices);

};

class ResidueNumberParsing:public AtomGroupParsingStrategy
{
    public:
        ResidueNumberParsing(AtomGroupParsingInput& input):AtomGroupParsingStrategy(input){};

        virtual void Parse(std::vector<int>& indices);
};

class AtomTypeParsing:public AtomGroupParsingStrategy
{
    public:
        AtomTypeParsing(AtomGroupParsingInput& input):AtomGroupParsingStrategy(input){};

        virtual void Parse(std::vector<int>& indices);
};

class ResidueNameParsing:public AtomGroupParsingStrategy
{
    public:
        ResidueNameParsing(AtomGroupParsingInput& input):AtomGroupParsingStrategy(input){};

        virtual void Parse(std::vector<int>& indices);
};

namespace AtomGroupParsingRegistry
{
    using Base = AtomGroupParsingStrategy;
    using Key  = std::string;

    using Factory = GenericFactory<Base,Key,AtomGroupParsingInput&>;

    template<typename D>
    using registry_ = RegisterInFactory<Base, D, Key, AtomGroupParsingInput&>;
}