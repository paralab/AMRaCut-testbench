# CMakeLists.txt
cmake_minimum_required(VERSION 3.20)  # Minimum CMake version required

project(Partitioning CXX)  # Project name

set(CMAKE_COLOR_MAKEFILE ON)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20) 
set(CMAKE_CXX_STANDARD_REQUIRED ON)


if(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE Release)

endif()



# compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -fopenmp -pedantic -Wall -Wextra -Wshadow -fdiagnostics-color=always")

set(CMAKE_CXX_FLAGS_DEBUG "-g -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(MAIN_TARGET_NAME "testbench_main")

add_executable(${MAIN_TARGET_NAME} src/testbench_main.cpp )


add_subdirectory(src)

include_directories(src)


set(DEPS_INSTALL_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../install")

find_library(AMRACUT_LIBRARY
    NAMES libamracut amracut
    PATHS "${DEPS_INSTALL_ROOT}/AMRaCut/lib"
    REQUIRED
    NO_DEFAULT_PATH)

find_library(GMSH_LIBRARY
    NAMES libgmsh gmsh
    PATHS "${DEPS_INSTALL_ROOT}/gmsh/lib"
    REQUIRED
    NO_DEFAULT_PATH)

find_package(MPI REQUIRED)


target_include_directories(${MAIN_TARGET_NAME} PRIVATE
    ${DEPS_INSTALL_ROOT}/AMRaCut/include
    ${DEPS_INSTALL_ROOT}/gmsh/include
)

target_link_libraries(${MAIN_TARGET_NAME} PRIVATE
    ${AMRACUT_LIBRARY}
    ${GMSH_LIBRARY}
    MPI::MPI_CXX
)


if(ENABLE_VTK_FEATURES)
    if(NOT VTK_INSTALL_DIR_PATH)
        message(FATAL_ERROR "The VTK_INSTALL_DIR_PATH is required when ENABLE_VTK_FEATURES is on")
    endif()
    message("VTK features enabled")
    add_compile_definitions(ENABLE_VTK_FEATURES=1)
    set(VTK_DIR "${VTK_INSTALL_DIR_PATH}")

    find_package(
        VTK REQUIRED
        COMPONENTS
            CommonCore
            CommonDataModel
            IOXML
    )

    target_link_libraries(${MAIN_TARGET_NAME} PRIVATE ${VTK_LIBRARIES})

    vtk_module_autoinit(
        TARGETS ${MAIN_TARGET_NAME}
        MODULES ${VTK_LIBRARIES}
      )

endif()