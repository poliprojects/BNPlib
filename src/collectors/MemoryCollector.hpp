#ifndef MEMORYCOLLECTOR_HPP
#define MEMORYCOLLECTOR_HPP

#include "BaseCollector.hpp"


class MemoryCollector: public BaseCollector {
protected:
    //! Deque that contains all states in Protobuf-object form
    std::deque<State> chain;

    //! Reads the next state, based on the curr_iter curson
    State next_state() override;

public:
    // DESTRUCTOR AND CONSTRUCTORS
    ~MemoryCollector() = default;
    MemoryCollector() = default;

    //! Initializes collector
    void start() override {return;}
    //! Closes collector
    void finish() override {return;}

    //! Writes the given state to the collector
    void collect(State iter_state) override;

    // GETTERS AND SETTERS
    //! Returns i-th state in the collector
    State get_state(unsigned int i) override {return chain[i];}
    //! Returns the whole chain in form of a deque of States
    std::deque<State> get_chain() override {return chain;}

};


#endif // MEMORYCOLLECTOR_HPP
