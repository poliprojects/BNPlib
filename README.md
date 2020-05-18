# BNPlib: A Nonparametric C++ Library
## Requirements
* Command line
* Python

## Preparation
* Download from Github by running:
  ```git clone https://github.com/poliprojects/BNPlib.git```
* Python packages: pybind11 (through ```pip install pybind11```)
* Compile a .proto file with protoc: first ```cd src/api``` then run:
  ```../../lib/protobuf/src/protoc -I=. --python_out=. output.proto```
* same with Python (TODO)

### Installing
* Compile by running ```make```
* Compile protobuf: run ```./install_protobuf.sh``` (or just
```install_protobuf.sh```, depending on your OS)
* Compile a .proto file with protoc: first ```cd src/api``` then run
```../../lib/protobuf/src/protoc -I=. --cpp_out=. output.proto```

### Documentation
* TBD

### Usage examples
* ```./maintest_uni csv/data_uni.csv neal2 memory```

#### main_factory_run
* ```./main_factory_run csv/data.csv "neal2" "memory"```
* ```./main_factory_run csv/data.csv "neal2" "file" ("filename")```


#### main_factory_dataless
* ```./main_factory_dataless "neal2_dataless" "memory"```
* ```./main_factory_dataless "neal2" "file" ("filename")```

### Code instructions
* Change number of total iterations (default 10000) and burn-in iterations
  (default 1000) with the appropriate setter methods of the algorithm object:
  ```sampler.set_maxiter(20000); sampler.set_burnin(2000);```
* In Neal8 algorithm, change number of auxiliary blocks (default 3) with the
  appropriate setter method of the algorithm object: ```sampler.set_n_aux(4);```
