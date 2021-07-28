#pragma once 
#include "CommonTypes.h"
#include "Assert.h"

#include <string>
#include <functional>
#include <map>


//  Useful resource: drdobbs.com/conversations-abstract-factory-template/184403786 
//  Useful resource on variadic templates (whatever that means): https://kevinushey.github.io/blog/2016/01/27/introduction-to-c++-variadic-templates/

template <class Base, typename Key=std::string, typename... Derived_Args>
class GenericFactory
{
    public:
        using createfunc = std::function<Base*(Derived_Args...)>;
        using FnRegistry   = std::map<Key, createfunc>;

        static GenericFactory& instance() {static GenericFactory _instance;return _instance;};
        void RegisterInFactory(const Key& key, createfunc func);
        Base* create(const Key& key, Derived_Args... args);
        const FnRegistry& get_registry(){return registry;}

    private:
        GenericFactory(){};

        // delete copy constructor and operator=
        GenericFactory(const GenericFactory& other) = delete;
        GenericFactory& operator=(const GenericFactory& other)=delete;

        // The static instance 
        static GenericFactory _instance;
        // The registry
        FnRegistry registry;
};

template <class Base, typename Key, typename... Derived_Args>
void GenericFactory<Base, Key, Derived_Args...>::RegisterInFactory(const Key& key, createfunc func)
{
    auto it = registry.find(key);
    ASSERT((it == registry.end()), "The key " << key << " is already registered in this factory.");

    registry.insert(std::pair<Key, createfunc>(key, func));
}

template <class Base, typename Key, typename... Derived_Args>
Base* GenericFactory<Base,Key,Derived_Args...>::create(const Key& key, Derived_Args... args)
{
    auto it = registry.find(key);
    ASSERT((it != registry.end()), "The key " << key << " is not found in this factory");

    return registry.find(key)->second(args...);
}

template<class Base, class Derived, typename Key=std::string, typename... Derived_Args>
class RegisterInFactory
{
    public:
        static Base* CreatInstance(Derived_Args... Args)
        {
            return new Derived(Args...);
        }

        RegisterInFactory(const Key& key)
        {
            GenericFactory<Base,Key,Derived_Args...>::instance().RegisterInFactory(key, CreatInstance);
        };
};