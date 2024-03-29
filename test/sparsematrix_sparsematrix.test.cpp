/**
 * @file sparsematrix_sparsematrix.test.cpp
 * @author Cheon Younghoe (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2023-07-08
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "sparsematrix_sparsematrix.hpp"
#include "gtest/gtest.h"

TEST(SPARSEMATRIX_TEST, init)
{
    yh::sparsematrix::SparseMatrix  sm;
    sm.ReadIJVFile("/media/PROJECT02/project/libsparsematrix/data/01.ijv");
    //sm.WriteIJVFile("01.out.ijv");

    EXPECT_NEAR(100.0, sm.GetValue(1,1), 1e-6);
    EXPECT_NEAR(-50.0, sm.GetValue(1,2), 1e-6);
    EXPECT_NEAR(-30.0, sm.GetValue(2,1), 1e-6);
    EXPECT_NEAR(200.0, sm.GetValue(2,2), 1e-6);

    EXPECT_EQ(4, sm.GetNNZ());
}

TEST(SPARSEMATRIX_TEST, GetNNZ)
{
    yh::sparsematrix::SparseMatrix  sm;
    sm.ReadIJVFile("/media/PROJECT02/project/libsparsematrix/data/01.ijv");

    EXPECT_EQ(4, sm.GetNNZ());
}

TEST(SPARSEMATRIX_TEST, SymmetricTest)
{
    yh::sparsematrix::SparseMatrix  sm;
    sm.SetSymmetric(true);
    sm.ReadIJVFile("/media/PROJECT02/project/libsparsematrix/data/01.ijv");

    EXPECT_NEAR(100.0, sm.GetValue(1,1), 1e-6);
    EXPECT_NEAR(-50.0, sm.GetValue(1,2), 1e-6);
    EXPECT_NEAR(0.0, sm.GetValue(2,1), 1e-6);
    EXPECT_NEAR(200.0, sm.GetValue(2,2), 1e-6);

    EXPECT_EQ(3, sm.GetNNZ());
}

TEST(SPARSEMATRIX_TEST, GetPardisoNIaJaA)
{
    yh::sparsematrix::SparseMatrix  sm;
    sm.ReadIJVFile("/media/PROJECT02/project/libsparsematrix/data/01.ijv");

    const auto [n, ia, ja, a] = sm.GetPardisoNIaJaA<long long, double>();

    std::cout << "n  : " << *n << "\n";
    std::cout << "ia :";
    for (auto pos = 0; pos <= *n; ++pos)
    {
        std::cout << " " << ia[pos];
    }
    std::cout << "\n";

    std::cout << "ja :";
    for (auto pos = 0; pos < *n; ++pos)
    {
        for (auto pos_1 = ia[pos]; pos_1 < ia[pos + 1]; ++pos_1)
            std::cout << " " << ja[pos_1];
    }
    std::cout << "\n";

    std::cout << "a  :";
    for (auto pos = 0; pos < *n; ++pos)
    {
        for (auto pos_1 = ia[pos]; pos_1 < ia[pos + 1]; ++pos_1)
            std::cout << " " << ((double *)a)[pos_1];
    }
    std::cout << "\n";
}

TEST(SPARSEMATRIX_TEST, GetPardiso64NIaJaA)
{
    yh::sparsematrix::SparseMatrix  sm;
    sm.ReadIJVFile("/media/PROJECT02/project/libsparsematrix/data/01.ijv");

    const auto [n, ia, ja, a] = sm.GetPardiso64NIaJaA();

    std::cout << "n  : " << *n << "\n";
    std::cout << "ia :";
    for (auto pos = 0; pos <= *n; ++pos)
    {
        std::cout << " " << ia[pos];
    }
    std::cout << "\n";

    std::cout << "ja :";
    for (auto pos = 0; pos < *n; ++pos)
    {
        for (auto pos_1 = ia[pos]; pos_1 < ia[pos + 1]; ++pos_1)
            std::cout << " " << ja[pos_1];
    }
    std::cout << "\n";

    std::cout << "a  :";
    for (auto pos = 0; pos < *n; ++pos)
    {
        for (auto pos_1 = ia[pos]; pos_1 < ia[pos + 1]; ++pos_1)
            std::cout << " " << ((double *)a)[pos_1];
    }
    std::cout << "\n";
}

TEST(SPARSEMATRIX_TEST, WriteBinReadBin)
{
    yh::sparsematrix::SparseMatrix  sm;
    sm.ReadIJVFile("/media/PROJECT02/project/libsparsematrix/data/01.ijv");

    //const auto [n, ia, ja, a] = sm.GetPardiso64NIaJaA();

    sm.WriteBinIaJaAFile(std::string("01"));

    yh::sparsematrix::SparseMatrix  sm_1;

    auto    n_1             = sm.GetRowSize();
    auto    ia_filename     = std::string("/media/PROJECT02/project/libsparsematrix/build/release/test/01.2.ia.bin");
    auto    ja_filename     = std::string("/media/PROJECT02/project/libsparsematrix/build/release/test/01.2.ja.bin");
    auto    a_filename      = std::string("/media/PROJECT02/project/libsparsematrix/build/release/test/01.2.a.bin");

    sm_1.ReadBinIaJaAFile(n_1, ia_filename, ja_filename, a_filename);

    const auto [n, ia, ja, a] = sm_1.GetPardiso64NIaJaA();

    std::cout << "n  : " << *n << "\n";
    std::cout << "nnz  : " << sm_1.GetNNZ() << "\n";
    std::cout << "ia :";
    for (auto pos = 0; pos <= *n; ++pos)
    {
        std::cout << " " << ia[pos];
    }
    std::cout << "\n";

    std::cout << "ja :";
    for (auto pos = 0; pos < *n; ++pos)
    {
        for (auto pos_1 = ia[pos]; pos_1 < ia[pos + 1]; ++pos_1)
            std::cout << " " << ja[pos_1];
    }
    std::cout << "\n";

    std::cout << "a  :";
    for (auto pos = 0; pos < *n; ++pos)
    {
        for (auto pos_1 = ia[pos]; pos_1 < ia[pos + 1]; ++pos_1)
            std::cout << " " << ((double *)a)[pos_1];
    }
    std::cout << "\n";
}
