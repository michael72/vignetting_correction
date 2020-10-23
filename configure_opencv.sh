pip3 install setuptools
pip3 install --user conan
source ~/.profile
mkdir -p build
cd build
conan install ../conanfile_opencv.txt --build=opencv --build=protobuf
#conan install ../conanfile_opencv.txt 
# sudo apt install cmake
