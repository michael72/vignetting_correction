
cmake_minimum_required(VERSION 3.12) 

project( vgncorr_lib )

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# include directories retrieved from conan
include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

set ( SOURCE_FILES
	  VignettingCorrection.cpp
	  VignettingCorrection.h 
	  GaussianBlur.h
	  GaussianBlur.cpp
	  ImageAlgo.h 
	  ImageAlgoBase.h
	  ImageAlgo_gil.h
	  ImageAlgo_opencv.h
	  vgncorr.h
	  vgncorr.cpp
	)

add_library( vgncorr_lib ${SOURCE_FILES} )
if(DEFINED ENV{USE_OPENCV}) 
  add_compile_definitions("USE_OPENCV")
endif()
target_include_directories( vgncorr_lib PUBLIC "${boost_SOURCE_DIR}/include" )
