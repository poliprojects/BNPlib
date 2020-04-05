#ifndef FACTORY_IMP_HPP
#define FACTORY_IMP_HPP

#include "Factory.hpp"

template <class AbstractProduct, class Builder_type>
class Factory{
public:
	void add(const std::string &name, Builder_type const &);
	auto create(const std::string &name) const; // get object
	static Factory &Instance();

private:
	using Container_type = std::map<std::string, Builder_type>;
	using QBuilder = std::function< std::unique_ptr<QuadratureRule>() >;
	Factory() = default;
	Factory(Factory const &) = delete;
	Factory& operator=(Factory const &) = delete;
	Container_type storage;
};

using RulesFactory = Factory<QuadratureRule, QBuilder>;
__attribute__((constructor))
static void registerRules(){
	RulesFactory &factory = RulesFactory::Instance();
	factory.add("Simpson", QBuilder<Simpson>() );
	factory.add("Trapezoidal", QBuilder<Trapezoidal>() );
}

#endif // FACTORY_IMP_HPP
