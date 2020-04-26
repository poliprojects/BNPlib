#ifndef MEMORYCOLLECTOR_HPP
#define MEMORYCOLLECTOR_HPP

#include "BaseCollector.hpp"


class MemoryCollector: public BaseCollector {
protected:
    std::deque<State> chains;

    State next_state() override;

public:
    void collect(State iteration_state) override;
    void finish() override;
    State get_state(unsigned int i) override;

    // Destructor and constructor
    virtual ~MemoryCollector() = default;
    // TODO constructor?
};

#include "MemoryCollector.imp.hpp"

#endif // MEMORYCOLLECTOR_HPP
