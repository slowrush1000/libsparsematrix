cmake_minimum_required(VERSION 3.20)
project(libsparsematrix VERSION 1.0.0 LANGUAGES CXX)

set(DOCS_DIR            "${CMAKE_SOURCE_DIR}/docs")
set(DOXYGEN_CONF_FILE   "doxygen.conf")

# doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
    configure_file(${DOCS_DIR}/${DOXYGEN_CONF_FILE} @ONLY)
    add_custom_target(docs
        ${DOXYGEN_EXECUTABLE} "${DOCS_DIR}/${DOXYGEN_CONF_FILE}"
        WORKING_DIRECTORY ${DOCS_DIR}
        COMMENT "Generate API documentation with Doxygen" VERBATIM
   )
endif(DOXYGEN_FOUND)
