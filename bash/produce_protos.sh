#!/usr/bin/env bash

lib/protobuf/src/protoc    --cpp_out=src/collectors chain_state.proto
lib/protobuf/src/protoc --python_out=src/python     chain_state.proto
