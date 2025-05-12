cmake_minimum_required(VERSION 3.3)

project(TCPolygonMatching)

set(CMAKE_CXX_STANDARD 17)

option(CXX "enable C++ compilation" ON)

if(CXX)
  enable_language(CXX)
endif()

#activate address sanitizing
if(CMAKE_BUILD_TYPE MATCHES Debug)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
  set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fsanitize=address")
endif()


list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Visual Studio compiler with static runtime libraries
if(MSVC AND MT)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT")
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /MTd")
endif()

find_package(GUROBI REQUIRED)

find_package(CGAL 6.0.1 EXACT REQUIRED)
set(CGAL_BUILD_TYPE Release)
message(STATUS "Found CGAL version ${CGAL_VERSION}")

find_package(Boost 1.83.0 EXACT REQUIRED COMPONENTS system thread)
if(Boost_FOUND)
  message(STATUS "Found Boost version ${Boost_VERSION}")
else()
  message(FATAL_ERROR "Boost 1.83 not found")
endif()

find_package(Threads REQUIRED)

include_directories(${GUROBI_INCLUDE_DIRS})

# ignore BOOST deprecated headers
add_definitions("-DBOOST_ALLOW_DEPRECATED_HEADERS")
add_definitions("-DBOOST_BIND_GLOBAL_PLACEHOLDERS")

# do not show 'copy' warning on build
add_compile_options(-Wno-nonnull)

# list source files here
set(sources
        src/main.cpp
        src/shapefile_io_operations.cpp
        src/localization.cpp
        src/graph_computations.cpp
        src/solution.cpp
        src/linear_program.cpp
        include/cgal_includes.h
        include/shapefile_io_operations.h
        include/localization.h
        include/graph_computations.h
        include/solution.h
        include/linear_program.h
        include/command_line_parser.h
        src/command_line_parser.cpp
        include/tree_computations.h
        src/tree_computations.cpp
        include/logger.h
        src/logger.cpp
        include/binary_io.h
        src/binary_io.cpp
        include/polygon_operations.h
        src/polygon_operations.cpp
        include/polygon_wh.h
        src/polygon_wh.cpp
        include/threading.h
        src/threading.cpp
)

add_executable(${CMAKE_PROJECT_NAME} ${sources})

if(CXX)
  set(CMAKE_CXX_STANDARD 11)
  target_link_libraries(${CMAKE_PROJECT_NAME} optimized ${GUROBI_CXX_LIBRARY}
                        debug ${GUROBI_CXX_DEBUG_LIBRARY})
endif()

target_link_libraries(${CMAKE_PROJECT_NAME} ${GUROBI_LIBRARY})
target_link_libraries(${CMAKE_PROJECT_NAME} CGAL::CGAL)
target_link_libraries(${CMAKE_PROJECT_NAME} Threads::Threads)
target_link_libraries(${CMAKE_PROJECT_NAME} Boost::system Boost::thread)

#uncomment if shapelib is not pre-built
#add_subdirectory(lib/shapelib EXCLUDE_FROM_ALL)
#target_link_libraries(${CMAKE_PROJECT_NAME} shp)

# link pre-built shapelib
set(SHAPELIB_DIR "${CMAKE_SOURCE_DIR}/lib/shapelib")  # Adjust if needed
find_library(SHAPELIB_LIB NAMES shp PATHS ${SHAPELIB_DIR}/build)

if(NOT SHAPELIB_LIB)
  message(FATAL_ERROR "shapelib not found")
endif()

include_directories(${SHAPELIB_DIR})  # Include shapelib headers
target_link_libraries(${CMAKE_PROJECT_NAME} ${SHAPELIB_LIB})
# end linking pre-built shapelib

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  include(FeatureSummary)
  feature_summary(WHAT ALL)
endif()
