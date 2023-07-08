#!/usr/bin/csh -f

mkdir build docs cmake data src test
cp ../libspf/.clang-format .
cp ../libspf/.gitignore .
cp ../libspf/CMakeLists.txt .
cp ../libspf/README.md .
cp ../libspf/build/build_setup.csh build
cp ../libspf/docs/CMakeLists.txt docs
cp ../libspf/docs/doxygen.conf docs
cp ../libspf/cmake/FindGZSTREAM.cmake cmake
cp ../libspf/data/star.spf data
cp ../libspf/src/CMakeLists.txt src
cp ../libspf/test/CMakeLists.txt test
