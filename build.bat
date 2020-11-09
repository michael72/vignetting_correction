if not exist build mkdir build
if not exist bin mkdir bin
cd build
set USE_OPENCV=
conan install .. --settings arch=x86_64
cmake .. -DCMAKE_GENERATOR_PLATFORM=x64
cmake --build . --config Release
copy cmd\bin\vgncorr_cmd.exe ..\bin\vgncorr.exe
cd ..
