#include <iostream>
#include <random>
#include <tuple>
#include <vector>
//#include "DistributionWrapper.hpp"


class B{
	private:
		double b;
		std::string str;
	public:
		B(int x) : b(2*x) {std::cout << "aaa" << std::endl;}
		double get_b(){return b;}
};


int main()
{
	B b(3);
	std::tuple<int,int> args = {3,1};
	std::tuple<B,int> foo(args);
	
	std::cout << std::get<0>(foo).get_b() << " " << std::get<1>(foo) << std::endl;
	return 0;
}


int main2()
{
	std::tuple<Normal,Normal> args = {Normal(2,1),Normal(2,9)};
	DistributionWrapper<Gamma,Normal,Normal> d(args);
	return 0;
}
