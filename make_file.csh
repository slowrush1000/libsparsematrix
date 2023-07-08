#!/usr/bin/csh -f

set file = $argv[1]

touch include/sparsematrix_${file}.hpp
touch     src/sparsematrix_${file}.cpp
touch    test/sparsematrix_${file}.test.cpp
