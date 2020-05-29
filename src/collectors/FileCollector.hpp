#ifndef FILECOLLECTOR_HPP
#define FILECOLLECTOR_HPP

#include "BaseCollector.hpp"


class FileCollector: public BaseCollector {
protected:
    //! Unix file descriptor for reading mode
    int infd;
    //! Unix file descriptor for writing mode
    int outfd;
    //! Pointer to a reading file stream
    google::protobuf::io::FileInputStream *fin;
    //! Pointer to a writing file stream
    google::protobuf::io::FileOutputStream *fout;
    //! Name of file from which read/write
    std::string filename;
    //! Flag that indicates if the collector is open in read-mode
    bool is_open_read = false;
    //! Flag that indicates if the collector is open in write-mode
    bool is_open_write;

    //! Opens collector in reading mode
    void open_for_reading();
    //! Terminates reading mode for the collector
    void close_reading();
    //! Reads the next state, based on the curr_iter curson
    State next_state() override;

public:
    // DESTRUCTOR AND CONSTRUCTORS
     ~FileCollector() {
        if(is_open_write){
            fout->Close();
            close(outfd);
        }
        if(is_open_read){
            fin->Close();
            close(infd);
        }
    }
    FileCollector() = default; // TODO serve?
    FileCollector(const std::string &filename_) : filename(filename_){}
    //! Initializes collector
    void start() override;
    //! Closes collector
    void finish() override;

    //! Writes the given state to the collector
    void collect(State iter_state) override;
    //! Returns i-th state in the collector
    State get_state(unsigned int i) override;
    //! Returns the whole chain in form of a deque of States
    std::deque<State> get_chain() override;
};


#endif // FILECOLLECTOR_HPP
