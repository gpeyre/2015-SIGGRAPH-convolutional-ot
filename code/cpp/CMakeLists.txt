project(ot)

cmake_minimum_required(VERSION 2.6)
cmake_policy(VERSION 2.6)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#---------#
# TO EDIT #
#---------#

option(USE_OPENMP "OpenMP?"      ON)
option(USE_DEBUG "Debug mode?"   ON)
option(USE_64BIT_ARCH "64 arch?" ON)

#---------#
# OPTIONS #
#---------#

if(USE_DEBUG)
    set(CMAKE_BUILD_TYPE Debug)
else(USE_DEBUG)
    set(CMAKE_BUILD_TYPE Release)
endif(USE_DEBUG)

if(USE_64BIT_ARCH)
    set(LINK_FLAGS -m64)
    add_definitions(-m64)
endif(USE_64BIT_ARCH)

#----------#
# Packages #
#----------#

if(USE_OPENMP)
    find_package(OpenMP)
endif(USE_OPENMP)

find_package(GLUT)
find_package(OpenGL)

find_package(Eigen3)
#-------#
# Files #
#-------#

set(OT_SRC  ./ot/ConvSolver.cpp ./ot/LinearSolver.cpp ./ot/TriMesh.cpp)
set(GUI_SRC ./gui/Camera.cpp ./gui/Image.cpp ./gui/Shader.cpp ./gui/Viewer.cpp)
set(SRCS ${GUI_SRC} ${OT_SRC} main.cpp)

#----------#
# INCLUDES #
#----------#

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wextra -Wno-unused-parameter -Wno-unused-variable -Wno-long-long")

include_directories(${GLUT_INCLUDE_DIR})
include_directories(${OPENGL_INCLUDE_DIR})
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/ot)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/gui)
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

#-------#
# Build #
#-------#

if(APPLE)
    set(CMAKE_SHARED_LIBRARY_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CXX_FLAGS} -framework Accelerate")
endif(APPLE)

if(USE_OPENMP AND OPENMP_FOUND)
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   ${OpenMP_C_FLAGS}"  )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_executable(${PROJECT_NAME} ${SRCS})
target_link_libraries(${PROJECT_NAME} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
