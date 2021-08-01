#pragma once

#include "Assert.h"
#include <map>
#include <string>

template<typename key, typename Value>
class Registry
{
    public:
        using map = std::map<key,Value>;
        using iterator = typename map::iterator;

        Registry(){};
        ~Registry(){};

        // throw if there is already a key in the map
        void insert(const key& key_, const Value& val);

        // throw if the key is not found
        const Value& find(const key& key_) const;
        Value& find(const key& key_);

        const iterator cbegin() const { return MapKeyToValue.cbegin();}
        const iterator cend() const {return MapKeyToValue.cend();}

        iterator begin() {return MapKeyToValue.begin();}
        iterator end() {return MapKeyToValue.end();}

    private:    
        map MapKeyToValue;
};

template<typename key, typename Value>
void Registry<key,Value>::insert(const key& key_, const Value& val)
{
    auto it = MapKeyToValue.find(key_);

    ASSERT(( it == MapKeyToValue.end()), "The key " << key_ << " is already found in the map.");

    MapKeyToValue.insert(std::make_pair(key_, val));
}

template<typename key, typename Value>
const Value& Registry<key,Value>::find(const key& key_) const
{
    auto it = MapKeyToValue.find(key_);

    ASSERT(( it != MapKeyToValue.end()), "The key " << key_ << " is not found in the map.");

    return it -> second;
}

template<typename key, typename Value>
Value& Registry<key,Value>::find(const key& key_)
{
    auto it = MapKeyToValue.find(key_);

    ASSERT(( it != MapKeyToValue.end()), "The key " << key_ << " is not found in the map.");

    return it -> second;
}
