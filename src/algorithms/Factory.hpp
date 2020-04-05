#ifndef FACTORY_HPP
#define FACTORY_HPP

#include "Neal2.hpp"
#include "Neal8.hpp"

using AlgoBuilder = std::function< std::unique_ptr<Algorithm>() >;

class Factory{
private:
	std::map<std::string, AlgoBuilder> storage;

	// Constructors
	Factory() = default;
	Factory(const Factory &f) = delete;
	Factory& operator=(const Factory &f) = delete;

public:
	static Factory &Instance();

	void add(const std::string &name, const AlgoBuilder &builder) {
		...
	}

	auto create(const std::string &name) const { // get object
		switch(name){
			case "neal2":
				return std::make_unique<Neal2>();
				break;
			case "neal8":
				return std::make_unique<Neal8>();
				break;
			default:
				return std::unique_ptr<Algorithm>();
		}
	}
};

#endif // FACTORY_HPP
