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

template <typename AbstractProduct,
    typename Builder=std::function<std::unique_ptr<AbstractProduct> ()> >
class Factory{
private:
    // Aliases
    using Builder_type = Builder;
    using Identifier = std::string;

    // Deleted constructors
    Factory() = default;
    Factory(const Factory &f) = delete;
    Factory& operator =(const Factory &f) = delete;

    std::map<Identifier,Builder_type> storage;

public:
    // Destructor
    ~Factory() = default;
    
    static Factory& Instance(){
        static Factory factory;
        return factory;
    }
  
    template<typename... Args>
    std::unique_ptr<AbstractProduct> create_agorithm_object(const Identifier &name,
        Args&&... args) const {
        auto f = storage.find(name);
        if(f == storage.end()){
            std::stringstream idAsString;
            idAsString<<name;
            std::string out="Identifier " + idAsString.str() +
            " is not stored in the factory";
            throw std::invalid_argument(out);
        }
        else {
            // Use of std::forward to forward arguments to the constructor
            //  return std::make_unique<AbstractProduct>(f->second(
                //std::forward<Args>(args)...));
            return f->second(std::forward<Args>(args)...);
        }
    }

  //! Register the given rule.
    void add(const Identifier &name, const Builder_type &func){
        auto f = storage.insert(std::make_pair(name, func));
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
