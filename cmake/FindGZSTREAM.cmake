# - Find gzstream
# Find the gzstream library
#
#  This module defines the following variables:
#     GZSTREAM_FOUND        - True if GZSTREAM_INCLUDE_DIR & GZSTREAM_LIBRARY are found
#     GZSTREAM_LIBRARIES    - Set when GZSTREAM_LIBRARY is found
#     GZSTREAM_INCLUDE_DIRS - Set when GZSTREAM_INCLUDE_DIR is found
#
#     GZSTREAM_INCLUDE_DIR  - where to find gzstream.h etc.
#     GZSTREAM_LIBRARY      - the gzstream library
#
#
# This file based on FindALSA.cmake
#
#=============================================================================
# Copyright 2009 Kitware, Inc.
# Copyright 2009 Philip Lowman <philip@yhbt.com>
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distributed this file outside of CMake, substitute the full
#  License text for the above reference.)

find_path( GZSTREAM_INCLUDE_DIR NAMES gzstream.h
    #           PATHS ${CMAKE_INCLUDE_PATH}/gzstream /usr/include/gzstream
               PATHS ${CMAKE_INCLUDE_PATH}/gzstream /usr/include
           DOC "The gzstream include directory"
)

find_library( GZSTREAM_LIBRARY NAMES gzstream
              DOC "The gzstream library"
)

find_library( GZSTREAM_LIBRARY_DEBUG gzstream_d )

# handle the QUIETLY and REQUIRED arguments and set GZSTREAM_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GZSTREAM DEFAULT_MSG GZSTREAM_LIBRARY GZSTREAM_INCLUDE_DIR)

if(GZSTREAM_FOUND)
  if ( GZSTREAM_LIBRARY_DEBUG )
    set(GZSTREAM_LIBRARIES optimized ${GZSTREAM_LIBRARY} debug ${GZSTREAM_LIBRARY_DEBUG})
  else ()
    set(GZSTREAM_LIBRARIES ${GZSTREAM_LIBRARY})
  endif ()
  set(GZSTREAM_INCLUDE_DIRS ${GZSTREAM_INCLUDE_DIR})
endif()

mark_as_advanced(GZSTREAM_INCLUDE_DIR GZSTREAM_LIBRARY GZSTREAM_LIBRARY_DEBUG)
