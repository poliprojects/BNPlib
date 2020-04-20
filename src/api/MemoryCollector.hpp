#ifndef MEMORYCOLLECTOR_HPP
#define MEMORYCOLLECTOR_HPP

#include "BaseCollector.hpp"


class MemoryCollector: public BaseCollector {
protected:
    std::deque<IterationOutput> chains;

public:
    void collect(IterationOutput iteration_state) override {
    	chains.push_back(iteration_state);
    }
    std::deque<IterationOutput> get_chains() override {return chains;};

    // Destructor
    virtual ~MemoryCollector() = default;
    // TODO constructor?
};


#endif // MEMORYCOLLECTOR_HPP
