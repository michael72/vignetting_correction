#mkdir -p build
cd build

cmake -E env CXXFLAGS="-DUSE_OPENCV" cmake .. 

cmake --build . --config Release
