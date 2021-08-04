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
    vecval.clear();

    if (vecstr != nullptr)
    {
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

    ASSERT((vecval.size() == dim), "In readArrayNumber, the read vector for key= " << key << " is " << vecval.size() << " while the required size is " << dim);

    if (vecNum == true)
    {
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