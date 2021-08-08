# CMake version.
cmake_minimum_required(VERSION 3.10)

# Project.
project(
    Precept
    VERSION 1.0.0 
    DESCRIPTION "C++ header only utilities"
    LANGUAGES CXX
    )

# Set C++17.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

option(BUILD_NATIVE "Build and tune for native hardware?" ON)

# Set release.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# If necessary, include CTest.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    include(CTest)
endif()

# Include GNUInstallDirs.
include(GNUInstallDirs)
message(STATUS "${PROJECT_NAME} ${PROJECT_VERSION}")
message(STATUS "INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "INSTALL_INCLUDEDIR: ${CMAKE_INSTALL_INCLUDEDIR}")
message(STATUS "INSTALL_LIBDIR: ${CMAKE_INSTALL_LIBDIR}")
message(STATUS "INSTALL_BINDIR: ${CMAKE_INSTALL_BINDIR}")

# Add Precept library.
add_library(
    ${PROJECT_NAME}
    INTERFACE
    )
target_include_directories(
    ${PROJECT_NAME}
    INTERFACE
    $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )

# Add Precept-graphics library.
add_library(
    ${PROJECT_NAME}-graphics
    STATIC
    src/AdjacencyLists.cpp
    src/AnimationCurve.cpp
    src/Color.cpp
    src/DynamicBbTree.cpp
    src/DynamicKdTree.cpp
    src/Image.cpp
    src/ImmutableBbTree.cpp
    src/ImmutableKdTree.cpp
    src/Microsurface.cpp
    src/WavefrontObj.cpp
    )
set_target_properties(
    ${PROJECT_NAME}-graphics
    PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED YES
    )
if (BUILD_NATIVE)
    target_compile_options(
        ${PROJECT_NAME}-graphics
        PRIVATE 
        "$<$<CONFIG:RELEASE>:-O3;-DNDEBUG;-march=native;-mtune=native>"
        )
endif()
target_link_libraries(
    ${PROJECT_NAME}-graphics
    ${PROJECT_NAME}
    )

install(
    TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}_Targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
install(
    TARGETS ${PROJECT_NAME}-graphics
    EXPORT ${PROJECT_NAME}_Targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    COMPONENT graphics
    )

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${PROJECT_NAME}ConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
    )
set(INCLUDE_DIR include)
configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    PATH_VARS INCLUDE_DIR
    INSTALL_DESTINATION
    ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake
    )

install(
    EXPORT ${PROJECT_NAME}_Targets
    FILE ${PROJECT_NAME}Targets.cmake
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake
    )
install(
    FILES
    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    DESTINATION
    "${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake"
    )

# Install headers.
install(
    DIRECTORY "${PROJECT_SOURCE_DIR}/include/"
    DESTINATION include
    )

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
    add_subdirectory(ctest)
endif()

