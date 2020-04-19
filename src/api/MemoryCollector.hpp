#ifndef MEMORYCOLLECTOR_HPP
#define MEMORYCOLLECTOR_HPP

#include "BaseCollector.hpp"


class MemoryCollector: public BaseCollector {
protected:
    std::deque<IterationOutput> chains;

   public:
    MemoryCollector() {}
    void collect(IterationOutput iteration_state) override {
    chains.push_back(iteration_state);
    }
    std::deque<IterationOutput> get_chains() override {return chains;};
    virtual ~MemoryCollector() = default;
};


#endif // MEMORYCOLLECTOR_HPP
