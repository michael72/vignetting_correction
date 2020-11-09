mkdir -p build
mkdir -p bin
cd build
cmake .. 
cmake --build . --config Debug
cp cmd/bin/vgncorr_cmd ../bin/vgncorr
cd ..
