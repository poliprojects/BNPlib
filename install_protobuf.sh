#!/usr/bin/env bash

# Install necessary tools
sudo apt install autoconf automake libtool curl make g++ unzip -y

# Get protobuf library
git clone --branch v3.11.0 https://github.com/protocolbuffers/protobuf.git

# Compile protobuf
cd lib/protobuf
./autogen.sh
./configure
make
make check
sudo make install
sudo ldconfig
cd ../..
