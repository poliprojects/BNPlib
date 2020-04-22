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
	unsigned int size = 0;
	unsigned int curr_iter = -1;

	virtual State next_state() = 0;

public:
	// Destructor
	// TODO constructor?
	virtual ~BaseCollector() = default;

    virtual void collect(State iteration_state) = 0;
    virtual State get_state(unsigned int i) = 0; 

    virtual void finish() = 0;
    // Getters and setters
    State get_next_state(){
        curr_iter++;
        
    	if(curr_iter > size){
    		throw std::out_of_range("Error: curr_iter > size in collector");
    	}
    	return next_state();
    }

    unsigned int get_size() const {return size;}
};


#endif // BASECOLLECTOR_HPP
