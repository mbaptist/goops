project(goops)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules/")

#Set local paths of third-party libraries
#cat++
set(CAT_INCLUDE_DIRS ${CAT_INCLUDE_DIRS} $ENV{HOME}/fcup/projects/cat/install/include/cat)
set(CAT_LIBRARY_DIRS ${CAT_LIBRARY_DIRS} $ENV{HOME}/fcup/projects/cat/install/lib)
#goops
set(GOOPS_INCLUDE_DIRS ${GOOPS_INCLUDE_DIRS} $ENV{HOME}/fcup/projects/goops/install/include/goops)
set(GOOPS_LIBRARY_DIRS ${GOOPS_LIBRARY_DIRS} $ENV{HOME}/fcup/projects/goops/install/lib)

add_subdirectory(goops)
add_subdirectory(test)

