#ifndef FILECOLLECTOR_IMP_HPP
#define FILECOLLECTOR_IMP_HPP

#include "FileCollector.hpp"


void FileCollector::open_for_reading() {
    infd = open(filename.c_str(), O_RDONLY);
    fin = new google::protobuf::io::FileInputStream(infd);
    is_open_read = true;
}


State FileCollector::next_state() {
    if (!is_open_read){
        open_for_reading();
    }

    State out;
    bool keep = google::protobuf::util::ParseDelimitedFromZeroCopyStream(
        &out, fin, nullptr);

    if (!keep){
        std::out_of_range("Error: surpassed EOF in FileCollector");
    }
    if(curr_iter == size-1){
        curr_iter = -1;
        fin->Close();
        close(infd);
        is_open_read = false;
    }
    return out;
}


void FileCollector::finish() {
    if(is_open_write){
        fout->Close();
        close(outfd);
        is_open_write = false;
    }
}


State FileCollector::get_state(unsigned int i) {
    State state;
    for(unsigned int k = 0; k < i+1; k++){
        state = get_next_state();
    }
    if(i<size-1){
        curr_iter=-1;
        fin->Close();
        close(infd);
        is_open_read = false;}
    return state;
}


void FileCollector::collect(State iteration_state) {
    bool success =
        google::protobuf::util::SerializeDelimitedToZeroCopyStream(
        iteration_state, fout);
    size++;
    if (!success){
        std::cout << "Writing in FileCollector failed" << std::endl;
    }
}


#endif // FILECOLLECTOR_IMP_HPP
