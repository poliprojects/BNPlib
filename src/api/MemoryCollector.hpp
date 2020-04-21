#ifndef MEMORYCOLLECTOR_HPP
#define MEMORYCOLLECTOR_HPP

#include "BaseCollector.hpp"


class MemoryCollector: public BaseCollector {
protected:
    std::deque<State> chains;

public:
    void collect(State iteration_state) override {
    	chains.push_back(iteration_state);
    }
    std::deque<State> get_chains() override {return chains;};

    // Destructor
    virtual ~MemoryCollector() = default;
    // TODO constructor?
};


#endif // MEMORYCOLLECTOR_HPP
