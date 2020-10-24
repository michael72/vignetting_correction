if not exist build mkdir build
cd build
set USE_OPENCV=
conan install .. --settings arch=x86_64
cmake .. -DCMAKE_GENERATOR_PLATFORM=x64
cmake --build . --config Release
cd ..