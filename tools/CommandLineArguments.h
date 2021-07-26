#pragma once
#include "InputParser.h"
#include "Assert.h"
#include "CommonTypes.h"

#include <vector>
#include <string>
#include <map>
#include <array>

class CommandLineArguments
{
    public:
        enum Keys {Required,Optional};

        CommandLineArguments(int argc, char** argv);
        ~CommandLineArguments(){};

        void ParseArgv();

        template<typename T>
        bool readValue(const std::string& key, const Keys, T& value) const;

        template<typename T, std::size_t dim>
        bool readArray(const std::string& key, const Keys, std::array<T,dim>& arr) const;

        template<typename T>
        bool readVector(const std::string& key,const Keys, std::vector<T>& vec) const;

        bool has_key(const std::string& key) const;

        bool readString(const std::string& key, const Keys type, std::string& str) const
        {
            return readValue<std::string>(key, type,str);
        }

        void print() const;

        // getters
        int get_num_keys() const {return num_keys_;}

    private:
        std::vector<std::string> argv_ = {};
        int argc_ = 0;

        // Number of actual keys contained in the command line
        int num_keys_ = 0;

        std::map<std::string, std::vector<std::string>> args_map_;
};
template<typename T>
bool CommandLineArguments::readValue(const std::string& key, const Keys type, T& val) const
{
    std::array<T,1> arr;

    bool Read = readArray<T, 1>(key, type, arr);

    if (Read == false)
    {
        return false;
    }

    val = arr[0];

    return true;
}

template<typename T, std::size_t dim>
bool CommandLineArguments::readArray(const std::string& key, const Keys type, std::array<T,dim>& arr) const
{
    std::vector<T> temp_vec = {}; 
    bool Read = readVector<T>(key, type, temp_vec);

    if (Read == false)
        return false;

    ASSERT((temp_vec.size() == dim), "The size of the array passed in is " << temp_vec.size() << ". However the required amount is " << dim);

    for (int i=0;i<dim;i++)
    {
        arr[i] = temp_vec[i];
    }  

    return true;
}

template<typename T>
bool CommandLineArguments::readVector(const std::string& key, const Keys type, std::vector<T>& vec) const
{
    auto it  = args_map_.find(key);

    // Check for requiredness
    if (type == Keys::Required)
    {
        ASSERT((it != args_map_.end()), "The key " << key << " is not found!");
    }
    if (type == Keys::Optional)
    {
        if (it == args_map_.end())
        {
            return false;
        }
    }

    std::vector<std::string> vecstr = it -> second;

    // Convert everything in a vector string to type T
    StringTools::VectorStringTransform<T>(vecstr, vec);

    return true;
}
