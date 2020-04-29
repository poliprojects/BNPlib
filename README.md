# BNPlib: A Nonparametric C++ Library
## Instructions
You will need a terminal to run the following commands.

### Downloading from GitHub
* Run ```git clone https://github.com/poliprojects/BNPlib.git```

### Installing
* Compile by running ```make```
* Compile protobuf: run ```./install_protobuf.sh``` (or just
```install_protobuf.sh```, depending on your OS)
* Compile a .proto file with protoc: first ```cd src/api``` then run
```../../lib/protobuf/src/protoc -I=. --cpp_out=. output.proto```

### Documentation
* TBD

### Usage examples
* ```./main path/to/filename.csv Collector``` and input required values to
  stdin (choose Collector: "file" or "memory")
* ```./main path/to/filename.csv "FileCollector" "filename"``` in case of FileCollector specify filename
* ```./main csv/data.csv <<< "5.0 0.1 2.0 2.0 1.0 3"``` to automatically fill
  in all values

#### main_factory_run
./main_factory_run csv/data.csv "neal2" "memory"
./main_factory_run csv/data.csv "neal2" "file" ("filename")


#### main_factory_dataless
./main_factory_dataless "neal2_dataless" "memory"
./main_factory_dataless "neal2" "file" ("filename")


### Code instructions
* Change number of total iterations (default 10000) and burn-in iterations
  (default 1000) with the appropriate setter methods of the algorithm object:
  ```sampler.set_maxiter(20000); sampler.set_burnin(2000);```
* In Neal8 algorithm, change number of auxiliary blocks (default 3) with the
  appropriate setter method of the algorithm object: ```sampler.set_n_aux(4);```
 
