# Configure:
# Set environment variable FETCHCONTENT_BASE_DIR to the desired download location for external libraries and headers
# Otherwise ~/.cmake-deps will be used if user directory is found - otherwise project/_deps is used.

cmake_minimum_required(VERSION 3.12) 

project( vgncorr )

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory (lib)
add_subdirectory (cmd)
