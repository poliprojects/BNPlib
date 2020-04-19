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
    BaseCollector(){}
    virtual void collect(IterationOutput iteration_state) = 0;
    virtual ~BaseCollector() = default;
    virtual std::deque<IterationOutput> get_chains() = 0;
};

#endif // BASECOLLECTOR_HPP
