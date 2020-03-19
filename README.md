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
* ```./main path/to/filename.csv``` and input required values to stdin
* ```./main csv/data.csv <<< "5.0 0.1 2.0 2.0 1.0 3"``` to automatically fill in
  all values
