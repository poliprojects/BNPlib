#include <fstream>


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

  virtual void collect(IterationOutput iteration_state) = 0;
};



class MemoryCollector : public BaseCollector {
protected:
std::deque<IterationOutput> chains;

public:
void collect(IterationOutput iteration_state) {
    //chains.push_back(iteration_state);
    }
};

class FileCollector: public BaseCollector {
protected:
int outfd;
//google::protobuf::io::FileOutputStream fout;

public:
FileCollector(std::string filename) {
   int outfd = open(filename.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0777);
   //fout=google::protobuf::io::FileOutputStream(outfd); 
}

 ~FileCollector() {
   //fout.Close();
   //close(outfd);
}

void collect(IterationOutput iteration_state) {
    //write_to_file(iteration_state) --> pseudocodice, da implementare
}

};






