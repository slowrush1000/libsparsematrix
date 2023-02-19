/**
 * @file    sparse_main.cpp
 * @author  Cheon Younghoe
 */

#include "sparse_main.hpp"
#include "sparse_sparsematrix.hpp"

/*
void func(int**   a);
*/

int
main(int argc, char* argv[])
{
    sparse::Sparsematrix    my_sparsematrix;
    my_sparsematrix.run(argc, argv);

    /*
    int*    a;
    func(&a);
    printf("a[0] : %d\n", a[0]);
    printf("a[1] : %d\n", a[1]);
    */

    return 0;
}

/*
void
func(int**   a)
{
    int*    b   = new int[100];
    *a          = b;
    b[0]        = 100;
    b[1]        = 1000;
}
*/
