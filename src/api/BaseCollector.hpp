#ifndef BASECOLLECTOR_HPP
#define BASECOLLECTOR_HPP

#include <fstream>
#include <string>
#include <deque>
#include <vector>
#include "output.pb.h"

#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>

#include <google/protobuf/text_format.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/util/delimited_message_util.h>


class BaseCollector {
protected:
	unsigned int size;
	unsigned int curr_iter = 0;

	virtual State get_next_state() = 0;

public:
	// Destructor
	// TODO constructor?
	virtual ~BaseCollector() = default;

    virtual void collect(State iteration_state) = 0;
    virtual std::deque<State> get_chains() = 0;

    State next(){
    	curr_iter += 1;
    	if(curr_iter >= size){
    		throw std::out_of_range("Error: curr_iter >= size in collector");
    	}
    	return get_next_state();
    }
};

#endif // BASECOLLECTOR_HPP
