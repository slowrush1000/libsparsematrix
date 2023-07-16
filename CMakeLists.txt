cmake_minimum_required(VERSION 3.20)
project(libsparsematrix VERSION 1.0.0 LANGUAGES CXX)

set(PROJECT_NAME  "libsparsematrix")
set(PROGRAM_NAME  "sparsematrix")

# Location of additional cmake modules
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" "${CMAKE_MODULE_PATH}")

# ycm
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# c++ standard : 20
set(CMAKE_CXX_STANDARD          20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# directory setting
set(EXTERN_DIR                  "/media/PROJECT02/extern")
set(INSTALL_DIR                 "${EXTERN_DIR}/${PROJECT_NAME}")

# source/include directories
set(PROJECT_SRC_DIR             "${CMAKE_SOURCE_DIR}/src")
set(PROJECT_INC_DIR             "${CMAKE_SOURCE_DIR}/include")

# ctags : make ctags
add_custom_target(ctags ALL
    COMMENT "Generate tags file with ctags"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMAND ctags -R --c++-kinds=+p --fields=+iaS --extra=+q ${CMAKE_SOURCE_DIR}
    COMMAND cp -f tags ${CMAKE_SOURCE_DIR}/src
)

# clangformat : make clang-format
add_custom_target(clang-format ALL
    COMMENT "Format source codes"
    COMMAND /usr/bin/clang-format -i ${CMAKE_SOURCE_DIR}/src/*.cpp
    COMMAND /usr/bin/clang-format -i ${CMAKE_SOURCE_DIR}/include/*.hpp
)

#
add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(docs)
