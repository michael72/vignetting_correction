pip3 install setuptools
pip3 install --user conan
source ~/.profile
mkdir -p build
cd build
conan install ..
# sudo apt install cmake
