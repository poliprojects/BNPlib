# BNPlib: A Nonparametric C++ Library

## Preparation
* Download: ```git clone https://github.com/poliprojects/BNPlib.git```
* Install libraries: ```bash/install_libs.sh```
* Compile library: ```make```
* Compile a .proto file with protoc to produce the classes needed (not
  necessary, files are already there): ```bash/produce_protos.sh```

### Documentation
* TODO

### Usage examples
* ```./maintest_uni csv/data_uni.csv neal2 memory```

### Code instructions
* Change number of total iterations (default 10000) and burn-in iterations
  (default 1000) with the appropriate setter methods of the algorithm object:
  ```sampler.set_maxiter(20000); sampler.set_burnin(2000);```
* In Neal8 algorithm, change number of auxiliary blocks (default 3) with the
  appropriate setter method of the algorithm object: ```sampler.set_n_aux(4);```
