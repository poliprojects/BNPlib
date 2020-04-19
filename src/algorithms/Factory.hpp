#ifndef FACTORY_HPP
#define FACTORY_HPP
#include <map>
#include <vector>
#include <memory>
#include <functional>
//#include <stdexcept>, #include <type_traits>, #include <sstream>


template<class AbstractProduct>
class Factory{
private:
    // Aliases
    using Builder = std::function< std::unique_ptr<AbstractProduct>() >;

    using Identifier = std::string;

    // Deleted constructors
    Factory() = default;
    Factory(const Factory &f) = delete;
    Factory& operator=(const Factory &f) = delete;

    // Storage for algorithm builders
    std::map<Identifier, Builder> storage;

public:
    // Destructor
    ~Factory() = default;

    // Method to create the factory
    static Factory& Instance(){
        static Factory factory;
        return factory;
    }

    template<typename... Args>
    std::unique_ptr<AbstractProduct> create_object(const Identifier &name,
        Args&&... args) const {
        auto f = storage.find(name);
        if(f == storage.end()){
            throw std::invalid_argument("Error: factory identifier not found");
        }
        else{
            return std::make_unique<AbstractProduct>(f->second( // TODO shared?
                std::forward<Args>(args)...));
            //return f->second(std::forward<Args>(args)...);
        }
    }

    void add_builder(const Identifier &name, const Builder &builder){
        auto f = storage.insert(std::make_pair(name, builder));
        if(f.second == false){
            std::cout <<
                "Warning: new duplicate builder was not added to factory" <<
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
