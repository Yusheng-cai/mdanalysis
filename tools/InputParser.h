#pragma once 
#include "Assert.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <sstream>
#include <utility>
#include <functional>
#include <array>
#include "Algorithm.h"

namespace StringTools
{
    template <typename T>
    T StringToType(std::string str);

    void to_lower(std::string& str);

    // Check if a string can be converted into a numeric number (float)
    bool isNumber(std::string str);

    template<typename T>
    bool VectorStringTransform(std::vector<std::string> vecstr, std::vector<T>& output);

    void RemoveBlankInString(std::string& str);

    // Checke if a line only consists of white spaces 
    bool CheckIfOnlyWhiteSpace(std::string& str);

    // parse the input string to numbers 
    // e.g. 1-2000:10
    void ConvertStringToIndices(const std::vector<std::string>& str, std::vector<int>& Indices);

    // sort a vector of ints
    void SortAndCheckDuplicates(std::vector<int>& indices);

    // strip a string of blank spaces
    std::string strip(const std::string& input);

    // split a sentence into a vector of strings
    std::vector<std::string> split(const std::string& input);

    // split a sentence but skip certain comment strings 
    std::vector<std::string> split(const std::string& input, const std::vector<std::string>& comment_str, bool ignore_after_comment=false);

    // Function that reads tabulated data by specified column
    template <typename T>
    void ReadTabulatedData(std::string filename, int col, std::vector<T>& data);

    // Function that reads tabulated data 
    template <typename T>
    void ReadTabulatedData(std::string filename, std::vector<std::vector<T>>& data);

    // function that reads the file extension 
    std::string ReadFileExtension(std::string filename);

    // read index files --> solid_like_atoms.index
    void ParseIndexFile(std::string filename, std::vector<std::vector<int>>& indices);
}


class ParameterPack
{
    public:
        ParameterPack():packname_("default"){};

        // Instantiate the parameter pack by name 
        ParameterPack(std::string packname):packname_(packname){};
        ~ParameterPack(){};

        enum KeyType{
            Required,
            Optional
        };

        // getters
        std::string get_packname() const {return packname_;};

        // insert into the parameterpack object
        std::string& insert(const std::string& key, const std::string& value);
        std::vector<std::string>& insert(const std::string& key, const std::vector<std::string>& value);
        ParameterPack& insert(const std::string& key, const ParameterPack& parampack);

        // Find all instances of values in Parameter Pack that matches with key
        // const function can be called by const + nonconst
        // non-const function can only be called by non-const 
        std::vector<const std::string*> findValues(const std::string& key, const KeyType) const;
        const std::string* findValue(const std::string& key, const KeyType) const;
        std::vector<const std::vector<std::string>*> findVectors(const std::string& key, const KeyType) const;
        const std::vector<std::string>* findVector(const std::string& key, const KeyType) const;
        std::vector<const ParameterPack*> findParamPacks(const std::string& key, const KeyType) const;
        const ParameterPack* findParamPack(const std::string& key, const KeyType) const;

        template <typename T>
        bool ReadNumber(const std::string& key, const KeyType, T& val) const;
        bool ReadString(const std::string& key, const KeyType, std::string& str) const;
        bool Readbool(const std::string& key, const KeyType, bool& boolean) const;

        template <typename T>
        bool ReadVectorNumber(const std::string& key, const KeyType, std::vector<T>& vecval) const;
        bool ReadVectorString(const std::string& key, const KeyType, std::vector<std::string>& vecstr) const;

        template<typename T, std::size_t dim>
        bool ReadArrayNumber(const std::string& key, const KeyType, std::array<T,dim>& arrval) const;
 
    private:
        std::multimap<std::string, std::string> value_;
        std::multimap<std::string, std::vector<std::string>> vectors_;
        std::multimap<std::string, ParameterPack> parampacks_;
        std::string packname_;
};

class TokenStream
{
    public:
        enum Status
        {
            Success, 
            Close_Brace, // }
            Open_Brace, // {
            Open_bracket, // [
            Close_bracket, // ]
            Failure,
            EndOfFile
        };

        TokenStream(std::ifstream& ifstream):ifstream_(ifstream){};

        // This is the most general way of reading tokens, we should read tokens 1 by 1 instead of reading
        // an entire line in first (former is more general) 
        Status ReadNextToken(std::string& token);
    private:   
        std::ifstream& ifstream_;
        std::stringstream line_stream_;

        std::string comment_str_="#";
};

class InputParser
{
    public:
        InputParser(){};
        ~InputParser(){};
        void ParseFile(const std::string& filename, ParameterPack& parampack);
        TokenStream::Status ParseNextToken(TokenStream& toks, ParameterPack& parampack);
        TokenStream::Status ParseParamPack(TokenStream& toks, ParameterPack& parampack);
        TokenStream::Status ParseVector(TokenStream& toks, std::vector<std::string>& vecvals);
};


template <typename T>
T StringTools::StringToType(std::string str)
{
    std::stringstream ss(str);
    T variable;
    ss >> variable;

    ASSERT((! ss.fail()), "Convert from string " << str <<" to number failed.");

    return variable;
}

template <typename T>
bool StringTools::VectorStringTransform(std::vector<std::string> vecstr, std::vector<T>& output)
{
    output.clear();
    for (auto str: vecstr)
    {
        output.push_back(StringTools::StringToType<T>(str));
    }

    return true;
}

template <typename T>
bool ParameterPack::ReadNumber(const std::string& key, const ParameterPack::KeyType keytype, T& val) const
{
    auto str = findValue(key, keytype);

    if (str != nullptr)
    {
        val = StringTools::StringToType<T>(*str);

        return true;
    }

    return false;
}

template <typename T>
bool ParameterPack::ReadVectorNumber(const std::string& key, const ParameterPack::KeyType keytype, std::vector<T>& vecval) const
{
    auto vecstr = findVector(key, keytype);

    // only clear the vector if the vector is found
    if (vecstr != nullptr)
    {
        vecval.clear();
        for (int i =0; i< vecstr->size();i++)
        {
            T val  = StringTools::StringToType<T>(vecstr->at(i));
            vecval.push_back(val);
        }

        return true;
    }

        
    return false;
}

template <typename T, std::size_t dim>
bool ParameterPack::ReadArrayNumber(const std::string& key, const ParameterPack::KeyType keytype, std::array<T,dim>& arrval) const
{
    std::vector<T> vecval;
    bool vecNum = ReadVectorNumber<T>(key,keytype, vecval);

    if (vecNum == true)
    {
        ASSERT((vecval.size() == dim), "In readArrayNumber, the read vector for key= " << key << " is " << vecval.size() << " while the required size is " << dim);
        for (int i=0;i<dim;i++)
        {
            arrval[i] = vecval[i];
        }

        return true;
    }
    else
    {
        return false;
    }
}

template <typename T>
void StringTools::ReadTabulatedData(std::string filename, int col, std::vector<T>& data)
{
    std::vector<std::vector<T>> total_data;
    ReadTabulatedData(filename, total_data);

    data.clear();
    data.resize(total_data.size());
    for (int i=0;i<total_data.size();i++)
    {
        data[i] = total_data[i][col];
    }
}


template <typename T>
void StringTools::ReadTabulatedData(std::string filename, std::vector<std::vector<T>>& data)
{
    data.clear();

    // read the filename 
    std::ifstream ifs;
    ifs.open(filename);

    ASSERT((ifs.is_open()), "The file with name " << filename << " is not opened.");

    std::string sentence;
    while(std::getline(ifs, sentence))
    {
        if (sentence.find("#") == std::string::npos)
        {
            std::stringstream ss;
            ss.str(sentence);

            std::vector<T> vecNum;
            T number;
            while (ss >> number)
            {
                vecNum.push_back(number);
            }

            data.push_back(vecNum);
        }
    }
}