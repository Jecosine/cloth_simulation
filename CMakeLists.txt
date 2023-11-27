cmake_minimum_required(VERSION 3.12)
project(cloth VERSION 1.0.0)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()
message("-- CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

# ##############################################################################
# enable C++  ##
# ##############################################################################
enable_language(CXX)

# ##############################################################################
# setup general C++ compilation flags#
# ##############################################################################

# add explictly debug symbols 
set(CMAKE_CXX_FLAGS " -g -O3")
# add g++ specific flags
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS}  -fpic -Wall -Wextra -Wno-deprecated-declarations -Wpedantic"
  )
endif()

# add icc specific flags
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -qopt-report=5 -qopt-report-phase=vec -qopt-report-file=stdout"
  )
endif()

# OpenGL
find_package(OpenGL REQUIRED)
message(STATUS "OpenGL included at ${OPENGL_INCLUDE_DIR}")

# GLUT
find_package(GLUT REQUIRED)
message(STATUS "GLUT included at ${GLUT_INCLUDE_DIR}")
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -framework GLUT")
endif()

# PAPI
find_program(PAPI_FOUND papi_avail)
if(PAPI_FOUND)
  message(STATUS "PAPI found")
  link_libraries(papi)
  add_compile_definitions(PAPI)
endif()

# OpenMP
find_package(OpenMP)
message(STATUS "OpenMP included at ${OPENMP_INCLUDE_DIR}")

include_directories(${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS})

add_executable(
  kernel_main
  kernel_main.cpp
  cloth_code.cpp
)

add_executable(
  kernel_opt
  kernel_main.cpp
  cloth_code_opt.cpp
)

add_executable(
  kernel_sse
  kernel_main.cpp
  cloth_code_sse.cpp
)

add_executable(
  kernel_vect_omp
  kernel_main.cpp
  cloth_code_vect_omp.cpp
)

add_executable(
  kernel_omp
  kernel_main_omp.cpp
  cloth_code_omp.cpp
)


set_target_properties(kernel_omp PROPERTIES COMPILE_FLAGS " -O3 " )
set_target_properties(kernel_sse PROPERTIES COMPILE_FLAGS " -O3 " )
set_target_properties(kernel_vect_omp PROPERTIES COMPILE_FLAGS " -O3 " )
set_target_properties(kernel_omp PROPERTIES COMPILE_FLAGS " -O3 " )

if (PAPI_FOUND) 
  add_executable(
    kernel_main_papi
    kernel_main_papi.cpp
    cloth_code.cpp)
  
  add_executable(
    kernel_opt_papi
    kernel_main_papi.cpp
    cloth_code_opt.cpp)

  add_executable(
    kernel_sse_papi
    kernel_main_papi.cpp
    cloth_code_sse.cpp
  )
  add_executable(
    kernel_vect_omp_papi
    kernel_main_papi.cpp
    cloth_code_vect_omp.cpp
  )
  add_executable(
    kernel_omp_papi
    kernel_main_omp_papi.cpp
    cloth_code_omp.cpp
  )
  set_target_properties(kernel_omp_papi PROPERTIES COMPILE_FLAGS " -O3 " )
  set_target_properties(kernel_sse_papi PROPERTIES COMPILE_FLAGS " -O3 " )
  set_target_properties(kernel_vect_omp_papi PROPERTIES COMPILE_FLAGS " -O3 " )
endif()



if(OpenMP_FOUND)
  target_link_libraries(kernel_main m OpenMP::OpenMP_CXX)
  target_link_libraries(kernel_opt m OpenMP::OpenMP_CXX)
  target_link_libraries(kernel_vect_omp m OpenMP::OpenMP_CXX)
  target_link_libraries(kernel_omp m OpenMP::OpenMP_CXX)
  target_link_libraries(kernel_omp_papi m OpenMP::OpenMP_CXX)
  target_link_libraries(kernel_vect_omp_papi m OpenMP::OpenMP_CXX)
  target_link_libraries(opengl_main m ${OPENGL_LIBRARIES} ${GLUT_LIBRARY}
                        OpenMP::OpenMP_CXX)
else()
  target_link_libraries(kernel_main m)
  target_link_libraries(kernel_opt m)
  target_link_libraries(kernel_sse m)
  target_link_libraries(kernel_sse_papi m)
  target_link_libraries(opengl_main m ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O3") 