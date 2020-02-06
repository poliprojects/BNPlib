Compile a .proto file with protoc:
```lib/protobuf/src/protoc -I=. --cpp_out=. output.proto```

Compile protobuf:
```
sudo apt install autoconf automake libtool curl make g++ unzip -y
git clone --branch v3.11.0 https://github.com/protocolbuffers/protobuf.git
cd lib/protobuf
./autogen.sh
./configure
make
make check
sudo make install
sudo ldconfig
cd ../..
```
