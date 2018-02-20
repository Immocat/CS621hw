project(non_rigid)

cmake_minimum_required(VERSION 2.8)

# setup output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# additional CMake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# setup GLUT
find_package(GLUT)
include_directories(${GLUT_INCLUDE_DIR})
if(NOT GLUT_FOUND)
    message(ERROR " GLUT not found")
endif()


#setup ceres solver
include_directories(${EIGEN_INCLUDE_DIR}) #this line is very important and necesary
find_package(Ceres REQUIRED)
message("CERES_LIBRARIES = ${CERES_LIBRARIES}")


# setup OpenGL
find_package(OpenGL)
include_directories(${OpenGL_INCLUDE_DIRS})
if(NOT OPENGL_FOUND)
    message(ERROR " OpenGL not found")
endif()

# setup OpenMesh
find_package(OpenMesh)
include_directories(${OPENMESH_INCLUDE_DIRS})
if(NOT OPENMESH_FOUND)
    message(ERROR " OpenMesh not found")
endif()
set_property(
    DIRECTORY
    APPEND PROPERTY COMPILE_DEFINITIONS _USE_MATH_DEFINES
)

# collect sources
file(GLOB all_sources *.cc)
file(GLOB all_headers *.hh)

file(GLOB ann_sources ANN/*.cpp)
file(GLOB ann_headers ANN/*.h)
#file(GLOB geo_headers geodesic/*.h)
#file(GLOB geo_sources geodesic/*.cpp)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")



include(FindOpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else(OPENMP_FOUND)
  message("ERROR: OpenMP could not be found.")
endif(OPENMP_FOUND)

# compile and link
add_library(ann ${ann_sources} ${ann_headers})
#add_library(geodesic ${geo_headers} ${geo_sources})
#set_target_properties(geodesic PROPERTIES LINKER_LANGUAGE CXX)
add_executable(non_rigid ${all_sources} ${all_headers})
#SET( CMAKE_CXX_FLAGS  "-pthread -fopenmp -std=c++11" )
target_link_libraries(non_rigid ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} ${OPENMESH_LIBRARIES} ann ${CERES_LIBRARIES})
