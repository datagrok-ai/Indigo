cmake_minimum_required(VERSION 3.6)

project(indigo-atomic LANGUAGES CXX)
set(DATA_PATH ${CMAKE_SOURCE_DIR}/data)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/common.h.in ${CMAKE_CURRENT_BINARY_DIR}/common.h)
file(GLOB_RECURSE ${PROJECT_NAME}_SOURCES CONFIGURE_DEPENDS 
${CMAKE_CURRENT_SOURCE_DIR}/**/*.cpp
${CMAKE_CURRENT_SOURCE_DIR}/**/*.h
)
add_executable (indigo-atomic 
 "run.cpp"
 "common.cpp"
 "common.h"
)

target_link_libraries (indigo-atomic 
                                 PRIVATE indigo-core
)
target_include_directories(${PROJECT_NAME}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
    PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
if(MSVC)
    target_link_options(${PROJECT_NAME}
            PRIVATE -force:multiple)
elseif(APPLE)
    target_link_options(${PROJECT_NAME}
            PRIVATE -Wl,-m)
elseif(MINGW OR UNIX OR MSYS OR CYGWIN)
    target_link_options(${PROJECT_NAME}
            PRIVATE -Wl,--allow-multiple-definition)
endif()