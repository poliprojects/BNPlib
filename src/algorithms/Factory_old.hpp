#ifndef FACTORY_HPP
#define FACTORY_HPP

template<template <class> class Hierarchy, class Hypers, class Mixture>
class Factory{
private:
	using AlgoBuilderType = std::function<
		std::unique_ptr<Algorithm<Hierarchy, Hypers, Mixture>>() >;

	std::map<std::string, AlgoBuilderType> storage;

public:
	static Factory& Instance(){
		static Factory factory;
		return factory;
	}

	void add_builder(const std::string &name, const AlgoBuilderType &builder) {
		storage[name] = builder;
	}

	auto create_algorithm(const std::string &name) const { // get object
		switch(name){
			case "neal2":
				return std::make_unique<Neal2<Hierarchy, Hypers, Mixture>>();
			case "neal8":
				return std::make_unique<Neal8<Hierarchy, Hypers, Mixture>>();
			default:
				return std::unique_ptr<Algorithm<Hierarchy, Hypers, Mixture>>();
		}
	}
};

#endif // FACTORY_HPP
