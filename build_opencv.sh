mkdir -p build
mkdir -p bin
cd build

cmake -E env CXXFLAGS="-DUSE_OPENCV" cmake .. 

cmake --build . --config Release

cp cmd/bin/vgncorr_cmd ../bin/vgncorr
cd ..
