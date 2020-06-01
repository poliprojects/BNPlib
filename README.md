# BNPlib: A Nonparametric C++ Library

## Installation
* Download: ```git clone https://github.com/poliprojects/BNPlib.git```
* Install auxiliary libraries: ```bash/install_libs.sh```
* Compile a .proto file with protoc to produce the classes needed: ```bash/produce_protos.sh```
* Compile library: ```make```
* Produce documentation (it will be saved in the ```doc``` dir): ```doxygen```

## Usage examples
* ```./maintest_uni csv/data_uni.csv neal2 memory```
* ```./maintest_multi csv/data_multi_2cl.csv neal2 memory```

## Code instructions (TODO?)
* Change number of total iterations (default 1000) and burn-in iterations
  (default 100) with the appropriate setter methods of the algorithm object: ```sampler.set_maxiter(2000); sampler.set_burnin(200);```
* In Neal8 algorithm, change number of auxiliary blocks (default 3) with the
  appropriate setter method of the algorithm object: ```sampler.set_n_aux(4);```
