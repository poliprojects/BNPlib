#ifndef FILECOLLECTOR_HPP
#define FILECOLLECTOR_HPP

#include "BaseCollector.hpp"


class FileCollector: public BaseCollector {
protected:
    int infd;
    int outfd;
    google::protobuf::io::FileInputStream* fin;
    google::protobuf::io::FileOutputStream *fout;
    bool is_open_read;
    bool is_open_write;

public:
    FileCollector(std::string filename) {
        int outfd = open(filename.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0777);
        fout = new google::protobuf::io::FileOutputStream(outfd);
    }

    void open_for_reading() { // TODO private?
        infd = open(filename.c_str(), O_RDWR); // TODO filename?
        fin = new google::protobuf::io::FileInputStream(infd);
        is_open_read = true;
    }

    std::deque<State> get_chains() override { // TODO ?
        std::cerr << "error" << std::endl;
        return std::deque<State>();
    }

    State get_next_state() override {
        if (!is_open_read){
            open_for_reading();
        }

        bool keep = true;
        bool clean_eof = true; // TODO WTF?
        State out;
        keep = google::protobuf::util::ParseDelimitedFromZeroCopyStream(
            &out, fin, nullptr);
     
        if (!keep){
            std::out_of_range("Error: surpassed EOF in FileCollector");
        }

        return out;
    }


    void collect(State iteration_state) override {
        bool success;
        success = google::protobuf::util::SerializeDelimitedToZeroCopyStream(
            iteration_state, fout);
        if (!success){
            std::cout << "Writing in FileCollector failed" << std::endl;
        }
    }

    // Destructor
    virtual ~FileCollector() {
        if(is_open_read){
            fin->Close();
            close(infd);
        }
        if(is_open_write){
            fout->Close();
            close(outfd);
        }
    }

    // TODO constructor?
};


#endif // FILECOLLECTOR_HPP
