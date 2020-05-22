
sudo apt install git
git clone https://github.com/poliprojects/BNPlib.git


(sudo apt-get install libprotobuf-dev protobuf-compiler) mi sa inutile

chmod +x install_protobuf.sh
./bash/install_protobuf.sh


 git clone https://github.com/stan-dev/math.git


##################
PYTHON

sudo apt update
sudo apt install software-properties-common
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt install python3.7

sudo apt install python3-pip
pip3 install pybind11

sudo apt-get install python3-matplotlib

make
python3 src/python/bnp_interface.py
python3 src/python/console.py

