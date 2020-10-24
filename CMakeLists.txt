# Configure:
# Set environment variable FETCHCONTENT_BASE_DIR to the desired download location for external libraries and headers
# Otherwise ~/.cmake-deps will be used if user directory is found - otherwise project/_deps is used.

cmake_minimum_required(VERSION 3.13) # >= 3.14 needed for FetchContent_MakeAvailable

project( vgncorr )

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

set ( SOURCE_FILES
	  VignettingCorrection.cpp
	  VignettingCorrection.h 
	  ImageAlgo.h 
	  ImageAlgo_gil.h
	  ImageAlgo_opencv.h )

add_executable( vgncorr ${SOURCE_FILES} )
#add_compile_definitions("USE_OPENCV")
target_include_directories( vgncorr PRIVATE "${boost_SOURCE_DIR}/include" )
target_link_libraries( vgncorr ${CONAN_LIBS} )
