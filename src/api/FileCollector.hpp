#ifndef FILECOLLECTOR_HPP
#define FILECOLLECTOR_HPP

#include "BaseCollector.hpp"


class FileCollector: public BaseCollector {
protected:
    // TODO docs on all of these:
    int infd;
    int outfd;

    // ...
    google::protobuf::io::FileInputStream *fin;
    google::protobuf::io::FileOutputStream *fout;
    std::string filename;
    // ...
    bool is_open_read=false;
    bool is_open_write;

    void open_for_reading() {
        infd = open(filename.c_str(), O_RDONLY); // TODO temp: filename?
        fin = new google::protobuf::io::FileInputStream(infd);
        is_open_read = true;
    }

    State next_state() override {
      
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

public:
    // Constructor and destructor
    FileCollector(std::string filename) : filename(filename){
        int outfd = open(filename.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0777);
        fout = new google::protobuf::io::FileOutputStream(outfd);
        is_open_write=true;
    }

    virtual ~FileCollector() {
        if(is_open_read){
            fin->Close();
            close(infd);
        }
        
    }

    void finish() override {

        if(is_open_write){
            fout->Close();
            close(outfd);
        }

    }

    std::deque<State> get_chains() override { // TODO ?
        std::cerr << "error" << std::endl;
        return std::deque<State>();
    }

    void collect(State iteration_state) override {
        bool success;
        success = google::protobuf::util::SerializeDelimitedToZeroCopyStream(
            iteration_state, fout);
        size++;
        if (!success){
            std::cout << "Writing in FileCollector failed" << std::endl;
        }
    }

};


#endif // FILECOLLECTOR_HPP
