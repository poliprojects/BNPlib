#ifndef FACTORY_HPP
#define FACTORY_HPP
#include <map>
#include <vector>
#include <memory>
#include <functional>
#include <stdexcept>
#include <type_traits>
#include <sstream>

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

    std::map<Identifier, AlgoBuilder> storage;

public:
    // Destructor
    ~Factory() = default;
    
    static Factory& Instance(){
        static Factory factory;
        return factory;
    }
  
    auto create_algorithm_object(const Identifier &name) const {
        switch(name){
            case "neal2":
                return std::make_unique<Neal2<Hierarchy, Hypers, Mixture>>();
            case "neal8":
                return std::make_unique<Neal8<Hierarchy, Hypers, Mixture>>();
            default:
                return std::unique_ptr<Algorithm<Hierarchy, Hypers, Mixture>>();
        }
    }

  //! Register the given rule.
    void add_builder(const Identifier &name, const AlgoBuilder &builder){
        auto f = storage.insert(std::make_pair(name, builder));
        if(f.second == false){
            std::stringstream idAsString;
            idAsString <<name;
            std::string message=std::string(
                "Double registration in Factory of id: ") + idAsString.str() +
                std::string(" is not allowed");
          throw std::invalid_argument(message);
        }
    }

    std::vector<Identifier> registered()const { // returns a list of reg. rules
        std::vector<Identifier> tmp;
        tmp.reserve(storage.size());
        for(auto i=storage.begin(); i!=storage.end();++i)
            tmp.push_back(i->first);
        return tmp;
    }

};


#endif // FACTORY_HPP
