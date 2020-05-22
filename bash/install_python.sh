#!/usr/bin/env bash

sudo apt update
sudo apt install software-properties-common
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt install python3.7

sudo apt install python3-pip
pip3 install pybind11

sudo apt install python3-matplotlib

make
python3 src/python/bnp_interface.py
python3 src/python/console.py
