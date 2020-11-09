if not exist build mkdir build
if not exist bin mkdir bin
cd build
set USE_OPENCV=1
conan install ..\conanfile_opencv.txt --settings arch=x86_64
cmake -E env CXXFLAGS="-DUSE_OPENCV" cmake .. -DCMAKE_GENERATOR_PLATFORM=x64 
cmake --build . --config Release

copy cmd\bin\vgncorr_cmd.exe ..\bin\vgncorr.exe
cd ..
