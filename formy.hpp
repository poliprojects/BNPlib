#ifndef H_FACTORY_HPP_
#define H_FACTORY_HPP_
#include <map>
#include <vector>
#include <memory>
#include <functional>
#include <stdexcept>
#include <type_traits>
#include <sstream>
namespace GenericFactory{
  
  /*! @brief A generic object factory.
    
    Freely inspired from the Alexandrescu Book.
    It is implemented as a Singleton. The compulsory way to 
    access a method is Factory::Instance().method().
    Typycally to access the factory one does
    \code
    auto&  myFactory = Factory<A,I,B>::Instance();
    myFactory.add(...)
    \endcode
    @tparam AbstractProduct The type of the base abstract class whose object are produced by the factory
    @tparam Identifier The identifier used to select the concrete object produced by the factory. It must have
            an ordering relation since they will be stored as key in a std::map
    @tparam Builder The type of the callable object that produces the the concrete object derived from AbstractProduct
    It may have arguments but its return type must be std::unique_ptr<AbstractProduct>.
  */
  template
  < typename AbstractProduct,
    typename Identifier,
    typename Builder=std::function<std::unique_ptr<AbstractProduct> ()>
    >
  class Factory{
  public:
    //! The container for the rules.
    using AbstractProduct_type=AbstractProduct;
    //! The identifier.
    /*  
	We must have an ordering since we use a map with 
	the identifier as key.
    */
    using  Identifier_type=Identifier;
    //! The builder type.
    /*!
      The default is a function wrapper that returns a unique_ptr<AbstractProduct>
      and with no arguments. 
    */
    using  Builder_type=Builder;
    //! Method to access the only instance of the factory
    static Factory & Instance();
    //! Get the rule with given name . 
    /*!
      The pointer is null if no rule is present.
      \tparam Args Automatically deduced parameter pack, it is the
      list of types of the builder arguments. It may be empty.
      \param name The identifier
      \param args The arguments pack to be passed to the builder. It may be empty
      
      \note The defualt builder does not take any argument. This version is rovided in case 
      the template is specialised with a builder that takes arguments.
    */
    template<typename... Args>
    std::unique_ptr<AbstractProduct> create(Identifier const & name, Args&&... args) const;
    //! Register the given rule.
    void add(Identifier const &, Builder_type const &);
    //! Returns a list of registered rules.
    std::vector<Identifier> registered()const;
    //! Unregister a rule.
    void unregister(Identifier const & name){ _storage.erase(name);}
    //! Destructor
    ~Factory()=default;
  private:
    typedef std::map<Identifier,Builder_type> Container_type;
    //! Made private since it is a Singleton
    Factory()=default;
    //! Deleted since it is a Singleton
    Factory(Factory const &)=delete;
    //! Deleted since it is a Singleton
    Factory & operator =(Factory const &)=delete;
    //! It contains the actual object factory.
    Container_type _storage;
  };


  
  template
  <
    typename AbstractProduct,
    typename Identifier,
    typename Builder
    >
  Factory<AbstractProduct,Identifier,Builder> & 
  Factory<AbstractProduct,Identifier,Builder>::Instance() {
    // We use the Meyer's trick to istantiate the factory as Singleton
    static Factory theFactory;
    return theFactory;
  }
  
  template
  <
    typename AbstractProduct,
    typename Identifier,
    typename Builder
    >
  template<typename... Args>
   std::unique_ptr<AbstractProduct>
  //  auto
  Factory<AbstractProduct,Identifier,Builder>::create(Identifier const & name, Args&&... args) 
    const {
    auto f = _storage.find(name);
    if (f == _storage.end())
      {
        std::stringstream idAsString;
        //I am assuming that identifier has a output streaming operator
        idAsString<<name;
	std::string out="Identifier " + idAsString.str() +
	  " is not stored in the factory";
	throw std::invalid_argument(out);
      }
    else
      {
        // Use of std::forward to forward arguments to the constructor
        //	return std::make_unique<AbstractProduct>(f->second(std::forward<Args>(args)...));
        return f->second(std::forward<Args>(args)...);
      }
  }
  
  template
  <
    typename AbstractProduct,
    typename Identifier,
    typename Builder
    >
  void  
  Factory<AbstractProduct,Identifier,Builder>::add(Identifier const & name, 
                                                   Builder_type const & func){
    auto f = 
      _storage.insert(std::make_pair(name, func));
    if (f.second == false)
      {
        std::stringstream idAsString;
        idAsString <<name;
        std::string message=std::string("Double registration in Factory of id: ")+idAsString.str()+
          std::string(" is not allowed");
        throw std::invalid_argument(message);
      }
  }
  
  template
  <
    typename AbstractProduct,
    typename Identifier,
    typename Builder
    >
  std::vector<Identifier>  Factory<AbstractProduct,Identifier,Builder>::registered() 
    const {
    std::vector<Identifier> tmp;
    tmp.reserve(_storage.size());
    for(auto i=_storage.begin(); i!=_storage.end();++i)
      tmp.push_back(i->first);
    return tmp;
  }
  
  
}// end namespace


#endif /* BC_FACTORY1_HPP_ */
