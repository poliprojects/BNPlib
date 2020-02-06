Compile a .proto file with protoc:
```lib/protobuf/src/protoc -I=. --cpp_out=. output.proto```

Compile protobuf:
```cd lib/protobuf```
```./configure && make```

sudo apt install autoconf automake libtool curl make g++ unzip -y
git clone ...
cd protobuf
./autogen.sh
./configure
sudo make
make check
make install
ldconfig
