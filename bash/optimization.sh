#!/usr/bin/env bash

#valgrind --tool=callgrind ./maintest_uni csv/data_uni.csv neal2 memory
valgrind --tool=callgrind ./maintest_multi csv/data_multi_2cl.csv neal2 memory
