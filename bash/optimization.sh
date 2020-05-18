#!/usr/bin/env bash

valgrind --tool=callgrind ./maintest_uni csv/data_uni.csv neal2 memory
