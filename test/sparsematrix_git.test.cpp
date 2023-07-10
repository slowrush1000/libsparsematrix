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

#include "sparsematrix_git.hpp"
#include "gtest/gtest.h"
#include <iostream>

TEST(GIT_TEST, init)
{
    std::cout << "k_GIT_BRANCH      : " << yh::sparsematrix::k_GIT_BRANCH << "\n";
    std::cout << "k_GIT_COMMIT_HASH : " << yh::sparsematrix::k_GIT_COMMIT_HASH << "\n";
    std::cout << "k_SOURCE_DIR      : " << yh::sparsematrix::k_SOURCE_DIR << "\n";
}
