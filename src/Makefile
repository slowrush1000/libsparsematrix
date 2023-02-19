
all 	: sparse_main.exe

sparse_main.exe	: sparse_sparsematrix.o sparse_main.o
	g++ -g -std=c++17 -o sparse_main.exe sparse_main.o sparse_sparsematrix.o -lgzstream -lz -lfmt

sparse_sparsematrix.o	: sparse_sparsematrix.cpp
	g++ -g -std=c++17 -c -o sparse_sparsematrix.o sparse_sparsematrix.cpp
sparse_main.o	: sparse_main.cpp
	g++ -g -std=c++17 -c -o sparse_main.o sparse_main.cpp

clean	:
	rm -rf *.o *.exe
