#ifndef FILECOLLECTOR_HPP
#define FILECOLLECTOR_HPP

#include "BaseCollector.hpp"


class FileCollector: public BaseCollector {
protected:
    int outfd;
    google::protobuf::io::FileOutputStream *fout;

public:
    FileCollector(std::string filename) {
        int outfd = open(filename.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0777);
        fout = new google::protobuf::io::FileOutputStream(outfd);
    }

    std::deque<State> get_chains() override { // TODO ?
        std::cerr << "error" << std::endl;
        return std::deque<State>();
    };

    void collect(State iteration_state) override {
        bool success;
        success = google::protobuf::util::SerializeDelimitedToZeroCopyStream(
            iteration_state, fout);
        if (!success){
            std::cout << "Writing Failed" << std::endl;
        }
    }

    // Destructor
    virtual ~FileCollector() {
        delete fout;
        close(outfd);
    }

    // TODO constructor?
};


#endif // FILECOLLECTOR_HPP
