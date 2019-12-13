#include <iostream>
#include <random>
#include <tuple>
#include <vector>
#include "DistributionWrapper.hpp"



int main()
{
	std::tuple<Normal,Normal> args = {Normal(2,1),Normal(2,9)};
	DistributionWrapper<Gamma,Normal,Normal> d(args);
	return 0;
}
