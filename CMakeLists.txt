cmake_minimum_required (VERSION 3.1)

project(Phoenix)

# compiler options
set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_CXX_COMPILER "mpic++")
set(CMAKE_VERBOSE_MAKEFILE off)
set(CMAKE_BUILD_TYPE Debug)

# standard flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -o3 -fopenmp")
set(CMAKE_INSTALL_PREFIX /home/brbass/research/owl/bin)

# gsl flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -L/home/software/rhel6/gsl/1.15/lib -lgsl -lgslcblas -I/home/software/rhel6/gsl/1.15/include")

# trilinos flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lepetra -lamesos -laztecoo -lepetraext -ltriutils -lepetra -lteuchosremainder -lteuchosnumerics -lteuchoscomm -lteuchosparameterlist -lteuchoscore -I/usr/cac/rhel6/trilinos/12.0.1/include -L/usr/cac/rhel6/trilinos/12.0.1/lib")

# enable_testing()

file(GLOB SOURCES "src/*.cc" "lib/*.cpp")

set(MAIN main/Phoenix.cc)

include_directories(src lib)

add_executable(phoenix ${SOURCES} ${MAIN})

install(TARGETS phoenix DESTINATION .)

