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
};



class MemoryCollector : public BaseCollector {
    protected:
    std::deque<IterationOutput> chains;

    public:
    MemoryCollector(){}
    void collect(IterationOutput iteration_state) override {
    chains.push_back(iteration_state);
    }
    std::deque<IterationOutput> get_chains(){return chains;};
    virtual ~MemoryCollector() = default;
};





class FileCollector: public BaseCollector {
protected:
int outfd;
google::protobuf::io::FileOutputStream *fout;

public:
FileCollector(std::string filename) {
   int outfd = open(filename.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0777);
   fout=new google::protobuf::io::FileOutputStream(outfd);

}



void collect(IterationOutput iteration_state) override {
        bool success;
        success = google::protobuf::util::SerializeDelimitedToZeroCopyStream(iteration_state, fout);
        if (! success)
            std::cout << "Writing Failed" << std::endl;
}

 virtual ~FileCollector() {
   delete fout;
   close(outfd);
}

};






