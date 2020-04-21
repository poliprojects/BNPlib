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
public:
	// Destructor
	// TODO constructor?
	virtual ~BaseCollector() = default;

    virtual void collect(State iteration_state) = 0;
    virtual std::deque<State> get_chains() = 0;
};

#endif // BASECOLLECTOR_HPP
