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
    void close_reading();
    State next_state() override;

public:
    // Constructor and destructor
    FileCollector() = default;
    FileCollector(std::string filename) : filename(filename){}

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

    void start() override;
    void finish() override;
    std::deque<State> get_chain() override;
    State get_state(unsigned int i) override;

    void collect(State iteration_state) override;
};


#endif // FILECOLLECTOR_HPP
