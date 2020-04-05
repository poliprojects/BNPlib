#ifndef FACTORY_UTILS_HPP
#define FACTORY_UTILS_HPP

#include "Factory.hpp"

template<class AlgoType>
std::unique_ptr<Algorithm> AlgoBuilder(){
	return std::make_unique<AlgoType>();
}

__attribute__((constructor))
static void add_known_algorithms(){
	Factory &factory = Factory::Instance();
	factory.add("neal2", AlgoBuilder<Neal2>() );
	factory.add("neal8", AlgoBuilder<Neal8>() );
}

#endif // FACTORY_UTILS_HPP
