Compile a .proto file with protoc:
```lib/protobuf/src/protoc -I=. --cpp_out=. output.proto```

Compile protobuf:
```cd lib/protobuf```
```./configure && make```

sudo apt install autoconf automake libtool curl make g++ unzip -y
git clone https://github.com/google/protobuf.git
cd protobuf
git submodule update --init --recursive
./autogen.sh
./configure
make
make check
sudo make install
sudo ldconfig
