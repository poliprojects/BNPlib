#ifndef FACTORY_UTILS_HPP
#define FACTORY_UTILS_HPP

#include "Factory.hpp"

template<class AlgoType>
std::unique_ptr<Algorithm<Hierarchy, Hypers, Mixture>> AlgoBuilder(){
	return std::make_unique<AlgoType>();
}

__attribute__((constructor))
static void add_known_algorithms(){
	Factory &factory = Factory::Instance();
	factory.add_builder("neal2",
		AlgoBuilder< Neal2<Hierarchy, Hypers, Mixture> >() );
	factory.add_builder("neal8",
		AlgoBuilder< Neal8<Hierarchy, Hypers, Mixture> >() );
}

#endif // FACTORY_UTILS_HPP
