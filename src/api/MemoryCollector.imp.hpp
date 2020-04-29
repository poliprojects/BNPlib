#ifndef MEMORYCOLLECTOR_IMP_HPP
#define MEMORYCOLLECTOR_IMP_HPP

#include "MemoryCollector.hpp"


State MemoryCollector::next_state() {
    if(curr_iter == size-1){
        curr_iter = -1;
        return chains[size-1];
    }
    else return chains[curr_iter];
}


void MemoryCollector::collect(State iteration_state) {
    chains.push_back(iteration_state);
    size++;
}


void MemoryCollector::finish() {
	return;
}


State MemoryCollector::get_state(unsigned int i) {
	return chains[i];
}


#endif // MEMORYCOLLECTOR_IMP_HPP
