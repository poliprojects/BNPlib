#!/usr/bin/env bash

cd src/api
../../lib/protobuf/src/protoc -I=. --cpp_out=. output.proto
# TODO protoc for Python
