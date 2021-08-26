#pragma once
#include "xdr/GroFile.h"
#include "tools/GenericFactory.h"

#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <algorithm>
#include <iostream>
#include <set>

struct AtomGroupParsingInput
{
    GroFile& grofile_;
    std::vector<std::string> selection_str;
};

class AtomGroupParsingStrategy
{
    public:
        AtomGroupParsingStrategy(AtomGroupParsingInput& input):grofile_(input.grofile_), selection_str_(input.selection_str){};
        virtual ~AtomGroupParsingStrategy(){};

        virtual void Parse(std::vector<int>& indices) = 0;
        virtual void update(std::vector<int>& indices){};
        std::string getType() const {return type_;}
        const std::set<int>& getResidueIndices() const {return ResidueIndices_;}

        void SortAndCheckNoDuplicate(std::vector<int>& indices);
    protected:
        GroFile& grofile_;
        std::vector<std::string>& selection_str_;
        std::vector<std::string> index_str_;
        std::set<int> ResidueIndices_;

        std::string type_;
};

class AtomIndexParsing:public AtomGroupParsingStrategy
{
    public:
        AtomIndexParsing(AtomGroupParsingInput& input):AtomGroupParsingStrategy(input){};
        virtual ~AtomIndexParsing(){};

        virtual void Parse(std::vector<int>& indices);

};

class ResidueNumberParsing:public AtomGroupParsingStrategy
{
    public:
        ResidueNumberParsing(AtomGroupParsingInput& input):AtomGroupParsingStrategy(input){};
        virtual ~ResidueNumberParsing(){};

        virtual void Parse(std::vector<int>& indices);
};

class AtomTypeParsing:public AtomGroupParsingStrategy
{
    public:
        AtomTypeParsing(AtomGroupParsingInput& input):AtomGroupParsingStrategy(input){};
        virtual ~AtomTypeParsing(){};

        virtual void Parse(std::vector<int>& indices);
};

class ResidueNameParsing:public AtomGroupParsingStrategy
{
    public:
        ResidueNameParsing(AtomGroupParsingInput& input):AtomGroupParsingStrategy(input){type_ = "ResidueName";};
        virtual ~ResidueNameParsing(){};

        virtual void Parse(std::vector<int>& indices);
};

class IndexFileParsing: public AtomGroupParsingStrategy
{
    public:
        IndexFileParsing(AtomGroupParsingInput& input);
        virtual ~IndexFileParsing(){};

        virtual void update(std::vector<int>& indices);
        virtual void Parse(std::vector<int>& indices);

        bool isOpen();

    private:
        std::string fileName_;

        std::ifstream ifs_;
        std::stringstream ss_;

        std::vector<std::vector<int>> Fileindices_;

        std::string comment_symbol = "#";

        int frame_count = 0;

        int totalFrames_;

        std::vector<int> Frames_;
};

namespace AtomGroupParsingRegistry
{
    using Base = AtomGroupParsingStrategy;
    using Key  = std::string;

    using Factory = GenericFactory<Base,Key,AtomGroupParsingInput&>;

    template<typename D>
    using registry_ = RegisterInFactory<Base, D, Key, AtomGroupParsingInput&>;
}