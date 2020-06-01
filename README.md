# BNPlib: A Nonparametric C++ Library
Authors: Bruno Guindani, Elena Zazzetti

## Requirements
* Command line, to run all following commands
* Python programming language (versions 3.7 and 3.8 were tested)

## Installation
* Download: ```git clone https://github.com/poliprojects/BNPlib.git```
* Install auxiliary libraries: ```bash/install_libs.sh```
* Compile a .proto file with protoc to produce the classes needed: ```bash/produce_protos.sh```
* Compile Python interface: ```make pybind_generate```
* Compile any test file: in Makefile change EXE value to the filename, then run ```make```
* Produce documentation (it will be saved in the ```doc``` dir): ```doxygen```

## Usage examples
* ```./maintest_uni csv/data_uni.csv neal2 memory```
* ```./maintest_multi csv/data_multi_2cl.csv neal8 file collector.recordio```
