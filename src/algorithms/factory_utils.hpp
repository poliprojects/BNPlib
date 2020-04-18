#ifndef FACTORY_UTILS_HPP
#define FACTORY_UTILS_HPP

#include "Factory.hpp"

using AlgoBuilderType = std::function<
	std::unique_ptr<Algorithm<Hierarchy, Hypers, Mixture>>() >;

template<class AlgoType>
std::unique_ptr<Algorithm<Hierarchy, Hypers, Mixture>> AlgoBuilder(){
	return std::make_unique<AlgoType>();
}

__attribute__((constructor))
static void add_known_algorithms(){
	Factory<Hierarchy, Hypers, Mixture> &factory = Factory::Instance();
	factory.add_builder("neal2",
		AlgoBuilder< Neal2<Hierarchy, Hypers, Mixture> >() );
	factory.add_builder("neal8",
		AlgoBuilder< Neal8<Hierarchy, Hypers, Mixture> >() );
}

void load_algo_factory(){
	auto &factory = Factory::Instance();
}

#endif // FACTORY_UTILS_HPP
