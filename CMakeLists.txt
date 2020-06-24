project(woodseer)
cmake_minimum_required(VERSION 2.4.6)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)




set(PROJECT_NAME woodseer)
set($ENV{CPATH} "")


find_package(PkgConfig)
find_package(OpenCV REQUIRED)

#set the default path for built executables to the "bin" directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
#set the default path for built libraries to the "lib" directory
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)


list(APPEND CMAKE_CXX_FLAGS "-fmax-errors=10 -pthread -std=c++11 ${CMAKE_CXX_FLAGS} -g -ggdb -Wall")


add_definitions(-g -ggdb -Wall -O1 -pthread -std=c++11 -DUSE_GSL)

# include all libs so far
include_directories(
    ${PROJECT_SOURCE_DIR}/include 
    ${OpenCV_INCLUDE_DIR}
)

add_executable(generate_surface_layers src/generate_surface_layers.cpp)
target_link_libraries(generate_surface_layers ${OpenCV_LIBRARIES})


