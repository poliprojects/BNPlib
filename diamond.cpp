#include <iostream>
#include <math.h>

class A {
public:
	virtual void foo() = 0;
};

class B: public A {
public:
	void foo() override {std::cout << "b" << std::endl;}
};

class C: public A {
public:
	void foo() override {std::cout << "c" << std::endl;}
};

class D: public B, public C {
public:
	void foo() override {
		B::foo();
		C::foo();
	}
};

int main(int argc, char const *argv[])
{
	D d;
	d.foo();
	return 0;
}
