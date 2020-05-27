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

    //! Flag that indicates if the collector is open in read-mode
    bool is_open_read = false;
    //! Flag that indicates if the collector is open in write-mode
    bool is_open_write;

    void open_for_reading();
    void close_reading();
    State next_state() override;

public:
    // DESTRUCTOR AND CONSTRUCTORS
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
    FileCollector() = default; // TODO serve?
    FileCollector(std::string filename) : filename(filename){} // TODO &
    //! Initializes collector
    void start() override;
    //! Closes collector
    void finish() override;

    //! Writes the given state to the collector
    void collect(State iteration_state) override;
    //! Returns i-th state in the collector
    State get_state(unsigned int i) override;
    //! Returns the whole chain in form of a deque of States
    std::deque<State> get_chain() override;
};


#endif // FILECOLLECTOR_HPP
