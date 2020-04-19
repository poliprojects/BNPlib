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

  static Factory& Instance(){
  // We use the Meyer's trick to istantiate the factory as Singleton
  static Factory theFactory;
  return theFactory;
}

  template<typename... Args>
  std::unique_ptr<AbstractProduct> create(Identifier const & name, Args&&... args) 
  const {
  auto f = storage.find(name);
  if (f == storage.end())
    {
      std::stringstream idAsString;
      //I am assuming that identifier has a output streaming operator
      idAsString<<name;
td::string out="Identifier " + idAsString.str() +
 " is not stored in the factory";
hrow std::invalid_argument(out);
    }
  else
    {
      return std::make_unique<AbstractProduct>(f->second(
        std::forward<Args>(args)...));
      //return f->second(std::forward<Args>(args)...);
    }
}
  //! Register the given rule.
  void add(Identifier const & name, 
                                                 Builder const & func){
  auto f = 
    storage.insert(std::make_pair(name, func));
  if (f.second == false)
    {
      std::stringstream idAsString;
      idAsString <<name;
      std::string message=std::string(
        "Double registration in Factory of id: ")+idAsString.str()+
        std::string(" is not allowed");
      throw std::invalid_argument(message);
    }
}
  //! Returns a list of registered rules.
  std::vector<Identifier> registered()const {
  std::vector<Identifier> tmp;
  tmp.reserve(storage.size());
  for(auto i=storage.begin(); i!=storage.end();++i)
    tmp.push_back(i->first);
  return tmp;
}
  //! Unregister a rule.
  void unregister(Identifier const & name){ storage.erase(name);}
  //! Destructor
  ~Factory()=default;
};


#endif // FACTORY_HPP
