# BNPlib: A Nonparametric C++ Library

## Installation
1) Download: ```git clone https://github.com/poliprojects/BNPlib.git```
2) Install auxiliary libraries: ```bash/install_libs.sh```
3) Compile a .proto file with protoc to produce the classes needed:
  ```bash/produce_protos.sh```
4) Compile library: ```make```
5) Produce documentation: ```doxygen``` (it will be saved in the ```doc``` dir)

## Usage examples
* ```./maintest_uni csv/data_uni.csv neal2 memory```
* ```./maintest_multi csv/data_multi_2cl.csv neal2 memory```

## Code instructions (TODO?)
* Change number of total iterations (default 10000) and burn-in iterations
  (default 1000) with the appropriate setter methods of the algorithm object:
  ```sampler.set_maxiter(20000); sampler.set_burnin(2000);```
* In Neal8 algorithm, change number of auxiliary blocks (default 3) with the
  appropriate setter method of the algorithm object: ```sampler.set_n_aux(4);```
