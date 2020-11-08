cmake_minimum_required(VERSION 3.12) 

project( vgncorr_cmd )

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

set ( SOURCE_FILES
	  Main.cpp )

add_executable( vgncorr_cmd ${SOURCE_FILES} )

target_link_libraries( vgncorr_cmd LINK_PUBLIC vgncorr_lib ${CONAN_LIBS})
