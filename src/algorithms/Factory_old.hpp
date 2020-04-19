#ifndef FACTORY_HPP
#define FACTORY_HPP
#include <map>
#include <vector>
#include <memory>
#include <functional>

#include "Neal2.hpp"
#include "Neal8.hpp"


template<template <class> class Hierarchy, class Hypers, class Mixture>
class Factory{
private:
    // Aliases
    using AlgoBuilder = std::function<
        std::unique_ptr<Algorithm<Hierarchy, Hypers, Mixture>>() >;
    using Identifier = std::string;

    // Deleted constructors
    Factory() = default;
    Factory(const Factory &f) = delete;
    Factory& operator =(const Factory &f) = delete;

    // Storage for algorithm builders
    std::map<Identifier, AlgoBuilder> storage;

public:
    // Destructor
    ~Factory() = default;
    
    static Factory& Instance(){
        static Factory factory;
        return factory;
    }

    template <typename ...Args>
    auto create_algorithm_object(const Identifier &name, Args&&... args) const {
        switch(name){
            case "neal2":
                return std::make_unique<Neal2<Hierarchy, Hypers, Mixture>>(
                std::forward<Args>(args)...);
            case "neal8":
                return std::make_unique<Neal8<Hierarchy, Hypers, Mixture>>(
                std::forward<Args>(args)...);
            default:
                return std::unique_ptr<Algorithm<Hierarchy, Hypers, Mixture>>(
                std::forward<Args>(args)...);
        }
    }

    void add_builder(const Identifier &name, const AlgoBuilder &builder){
        auto f = storage.insert(std::make_pair(name, builder));
        if(f.second == false){
            std::cout <<
                "Warning: new duplicate builder was not added to storage" <<
                std::endl;
        }
    }

    std::vector<Identifier> list_of_known_builders() const {
        std::vector<Identifier> tmp;
        tmp.reserve(storage.size());
        for(auto i = storage.begin(); i != storage.end(); i++){
            tmp.push_back(i->first);
        }
        return tmp;
    }

};


#endif // FACTORY_HPP
