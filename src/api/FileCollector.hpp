#ifndef FILECOLLECTOR_HPP
#define FILECOLLECTOR_HPP

#include "BaseCollector.hpp"


class FileCollector: public BaseCollector {
protected:
    // Docs on all of these:
    int infd;
    int outfd;

    // ...
    google::protobuf::io::FileInputStream *fin;
    google::protobuf::io::FileOutputStream *fout;
    std::string filename;

    // ...
    bool is_open_read = false;
    bool is_open_write;

    void open_for_reading();

    State next_state() override;

public:
    // Constructor and destructor
    FileCollector()=default;
    FileCollector(std::string filename) : filename(filename){
        int outfd = open(filename.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0777);
        fout = new google::protobuf::io::FileOutputStream(outfd);
        is_open_write = true;
    }

     ~FileCollector() {
        if (is_open_write){
            fout->Close();
            close(outfd);
        }
        if(is_open_read){
            fin->Close();
            close(infd);
        } 
    }

    void finish() override;

    State get_state(unsigned int i) override;

    void collect(State iteration_state) override;

};

#include "FileCollector.imp.hpp"

#endif // FILECOLLECTOR_HPP
