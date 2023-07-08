cmake_minimum_required(VERSION 3.20)
project(libsparsematrix VERSION 1.0.0 LANGUAGES CXX)

# git version
execute_process(
    COMMAND git rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
    COMMAND git log -1 --format=%h
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
    COMMAND pwd
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE SOURCE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
configure_file(${PROJECT_INC_DIR}/${PROGRAM_NAME}_git.hpp.in ${PROJECT_INC_DIR}/${PROGRAM_NAME}_git.hpp)

# program name/build dir
set(RUNTIME_OUTPUT_DIRECTORY                    "${CMAKE_SOURCE_DIR}/build/release")
set(ARCHIVE_OUTPUT_DIRECTORY                    "${CMAKE_SOURCE_DIR}/build/release")
if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    message("-- Debug build.")
    set(RUNTIME_OUTPUT_DIRECTORY                "${CMAKE_SOURCE_DIR}/build/debug")
    set(ARCHIVE_OUTPUT_DIRECTORY                "${CMAKE_SOURCE_DIR}/build/debug")
    set(PROJECT_EXE_NAME                        "${PROGRAM_NAME}.dbg.exe")
    set(PROJECT_LIB_NAME                        "${PROGRAM_NAME}.dbg")
elseif(${CMAKE_BUILD_TYPE} STREQUAL "Release")
    message("-- Release build.")
    set(RUNTIME_OUTPUT_DIRECTORY                "${CMAKE_SOURCE_DIR}/build/release")
    set(ARCHIVE_OUTPUT_DIRECTORY                "${CMAKE_SOURCE_DIR}/build/release")
    set(PROJECT_EXE_NAME                        "${PROGRAM_NAME}.exe")
    set(PROJECT_LIB_NAME                        "${PROGRAM_NAME}")
endif()

# sources/headers
set(PROJECT_MAIN_SOURCE                         "${PROJECT_SRC_DIR}/${PROGRAM_NAME}_main.cpp")
set(PROJECT_MAIN_HEADER                         "${PROJECT_INC_DIR}/${PROGRAM_NAME}_main.hpp")
file(GLOB PROJECT_SOURCES LIST_DIRECTORIES true "${PROJECT_SRC_DIR}/${PROGRAM_NAME}_*.cpp")
file(GLOB PROJECT_HEADERS LIST_DIRECTORIES true "${PROJECT_INC_DIR}/${PROGRAM_NAME}_*.hpp")
list(REMOVE_ITEM PROJECT_SOURCES                "${PROJECT_MAIN_SOURCE}")
list(REMOVE_ITEM PROJECT_HEADERS                "${PROJECT_MAIN_HEADER}")

# build : make -j
include_directories(${PROJECT_INC_DIR})
add_library(${PROJECT_LIB_NAME} STATIC           ${PROJECT_SOURCES})
add_executable(${PROJECT_EXE_NAME}               ${PROJECT_MAIN_SOURCE})
target_link_libraries(${PROJECT_EXE_NAME}        ${PROJECT_LIB_NAME})
set_target_properties(${PROJECT_EXE_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${RUNTIME_OUTPUT_DIRECTORY})
set_target_properties(${PROJECT_LIB_NAME} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${RUNTIME_OUTPUT_DIRECTORY})

# copy compile_commands.json to project root directory
add_custom_target(
    copy-compile-commands ALL
    ${CMAKE_COMMAND} -E copy_if_different 
        ${RUNTIME_OUTPUT_DIRECTORY}/compile_commands.json 
        ${CMAKE_SOURCE_DIR}
)

# yh::util
set(YHUTIL_DIR                                  "${EXTERN_DIR}/libutil")
set(YHUTIL_INC_DIR                              "${YHUTIL_DIR}/include")
set(YHUTIL_LIBRARIES                            "${YHUTIL_DIR}/lib/libutil.a")
include_directories(${YHUTIL_INC_DIR})
target_link_libraries(${PROJECT_EXE_NAME}        ${YHUTIL_LIBRARIES})

## openmp
#find_package(OpenMP)
#if(OPENMP_FOUND)
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#endif(OPENMP_FOUND)

# spdlog
set(SPDLOG_DIR                              "${EXTERN_DIR}/spdlog")
set(SPDLOG_INC_DIR                          "${SPDLOG_DIR}/include")
set(SPDLOG_LIBRARIES                        "${SPDLOG_DIR}/build/libspdlog.a")
include_directories(${SPDLOG_INC_DIR})
target_link_libraries(${PROJECT_EXE_NAME}        ${SPDLOG_LIBRARIES})

# gzstream
#find_package(GZSTREAM)
#if (GZSTREAM_FOUND)
#    include_directories(${GZSTREAM_INCLUDE_DIRS})
#    target_link_libraries(${PROJECT_EXE_NAME} ${GZSTREAM_LIBRARIES})
##    target_link_libraries(${PROJECT_TEST_EXE_NAME} ${GZSTREAM_LIBRARIES})
#    target_link_libraries(${PROJECT_EXE_NAME} z)
##    target_link_libraries(${PROJECT_TEST_EXE_NAME} z)
#endif(GZSTREAM_FOUND)

# fmt
find_package(fmt)
if(fmt_FOUND)
    target_link_libraries(${PROJECT_EXE_NAME} fmt::fmt)
endif(fmt_FOUND)

## tcl/tk
#find_package(TCL)
#if(TCL_FOUND)
#    include_directories(${TCL_INCLUDE_PATH})
#    target_link_libraries(${PROJECT_EXE_NAME} ${TCL_LIBRARY})
#endif(TCL_FOUND)

## tomlplusplus
#include_directories("../${TOMLPLUSPLUS_DIR}/include")

# install : make install
install(TARGETS     ${PROJECT_EXE_NAME} RUNTIME     DESTINATION "${INSTALL_DIR}/bin")
install(TARGETS     ${PROJECT_LIB_NAME} ARCHIVE     DESTINATION "${INSTALL_DIR}/lib")
install(FILES       ${PROJECT_HEADERS}              DESTINATION "${INSTALL_DIR}/include")
install(DIRECTORY   ${CMAKE_SOURCE_DIR}/docs/html   DESTINATION "${INSTALL_DIR}/docs")
