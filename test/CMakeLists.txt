cmake_minimum_required(VERSION 3.20)
project(libsparsematrix.test VERSION 1.0.0 LANGUAGES CXX)

#set(PROGRAM_NAME                                "sparsematrix")

# program name/build dir
set(RUNTIME_OUTPUT_DIRECTORY                    "${CMAKE_SOURCE_DIR}/build/release")
set(ARCHIVE_OUTPUT_DIRECTORY                    "${CMAKE_SOURCE_DIR}/build/release")
if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    message("-- Debug build.")
    set(RUNTIME_OUTPUT_DIRECTORY                "${CMAKE_SOURCE_DIR}/build/debug")
    set(ARCHIVE_OUTPUT_DIRECTORY                "${CMAKE_SOURCE_DIR}/build/debug")
    set(PROJECT_TEST_EXE_NAME                   "${PROGRAM_NAME}.dbg.test.exe")
    set(PROJECT_TEST_LIB_NAME                   "${PROGRAM_NAME}.dbg.test")
elseif(${CMAKE_BUILD_TYPE} STREQUAL "Release")
    message("-- Release build.")
    set(RUNTIME_OUTPUT_DIRECTORY                "${CMAKE_SOURCE_DIR}/build/release")
    set(ARCHIVE_OUTPUT_DIRECTORY                "${CMAKE_SOURCE_DIR}/build/release")
    set(PROJECT_TEST_EXE_NAME                   "${PROGRAM_NAME}.test.exe")
    set(PROJECT_TEST_LIB_NAME                   "${PROGRAM_NAME}.test")
endif()

# sources/headers
set(PROJECT_TEST_MAIN_SOURCE                    "${PROGRAM_NAME}_main.test.cpp")
file(GLOB PROJECT_TEST_SOURCES LIST_DIRECTORIES true "${PROGRAM_NAME}_*.test.cpp")
list(REMOVE_ITEM PROJECT_TEST_SOURCES "{PROJECT_TEST_MAIN_SOURCE}")

# googletest : make & ./program_name.test.exe
include_directories("${EXTERN_DIR}/googletest/googletest/include")
include_directories("${CMAKE_SOURCE_DIR}/include")
set(GOOGLETEST_LIBRARIES                        "${EXTERN_DIR}/googletest/build/lib/libgtest.a")
add_executable(${PROJECT_TEST_EXE_NAME}         ${PROJECT_TEST_MAIN_SOURCE} ${PROJECT_TEST_SOURCES})
target_link_libraries(${PROJECT_TEST_EXE_NAME}  ${PROJECT_LIB_NAME})
target_link_libraries(${PROJECT_TEST_EXE_NAME}  ${GOOGLETEST_LIBRARIES})
set_target_properties(${PROJECT_TEST_EXE_NAME}  PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${RUNTIME_OUTPUT_DIRECTORY})
set_target_properties(${PROJECT_TEST_EXE_NAME}  PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${RUNTIME_OUTPUT_DIRECTORY})

add_custom_target(test
    COMMENT "Run unit test"
    COMMAND ${RUNTIME_OUTPUT_DIRECTORY}/${PROJECT_TEST_EXE_NAME} --gtest_output=xml:${RUNTIME_OUTPUT_DIRECTORY}/test.result.xml
)

# yh::sparsematrix
if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(SPF_LIBRARIES                           "sparsematrix.dbg")
elseif(${CMAKE_BUILD_TYPE} STREQUAL "Release")
    set(SPF_LIBRARIES                           "sparsematrix")
endif()
include_directories(${PROJECT_INC_DIR})
target_link_libraries(${PROJECT_TEST_EXE_NAME}  ${SPF_LIBRARIES})

# yh::util
set(YHUTIL_DIR                                  "${EXTERN_DIR}/libutil")
set(YHUTIL_INC_DIR                              "${YHUTIL_DIR}/include")
set(YHUTIL_LIBRARIES                            "${YHUTIL_DIR}/lib/libutil.a")
include_directories(${YHUTIL_INC_DIR})
target_link_libraries(${PROJECT_TEST_EXE_NAME}  ${YHUTIL_LIBRARIES})

# spdlog
set(SPDLOG_DIR                                  "${EXTERN_DIR}/spdlog")
set(SPDLOG_INC_DIR                              "${SPDLOG_DIR}/include")
set(SPDLOG_LIBRARIES                            "${SPDLOG_DIR}/build/libspdlog.a")
include_directories(${SPDLOG_INC_DIR})
target_link_libraries(${PROJECT_TEST_EXE_NAME}  ${SPDLOG_LIBRARIES})

# gzstream
find_package(GZSTREAM)
if (GZSTREAM_FOUND)
    include_directories(${GZSTREAM_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_TEST_EXE_NAME} ${GZSTREAM_LIBRARIES})
    target_link_libraries(${PROJECT_TEST_EXE_NAME} z)
endif(GZSTREAM_FOUND)

# fmt
find_package(fmt)
if(fmt_FOUND)
    #    include_directories(${fmt_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_TEST_EXE_NAME} fmt::fmt)
endif(fmt_FOUND)

#enable_testing()
#if (NOT (TARGET gtest AND TARGET gtest_main))
#    message(FATAL_ERROR "gtest and gtest_main targets NOT found")
#endif()
#
#add_executable(${TEST_PROJECT_NAME} ${LIBRARY_TESTS_SOURCE})
#target_link_libraries(${TEST_PROJECT_NAME} 
#    gtest 
#    gtest_main 
#    ${LIBRARY_NAME} # NOTE: This is defined from project above
#)
#add_test(${TEST_PROJECT_NAME} ${TEST_PROJECT_NAME})
