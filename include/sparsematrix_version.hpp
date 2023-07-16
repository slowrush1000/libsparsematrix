/**
 * @file sparsematrix_version.hpp
 * @author Cheon Younghoe (you@dosparsematrix.com)
 * @brief
 * @version 0.1
 * @date 2023-07-08
 *
 * @copyright Copyright (c) 2023
 *
 */

#ifndef YH_SPARSEMATRIX_VERSION_H_
#define YH_SPARSEMATRIX_VERSION_H_

#include "sparsematrix_git.hpp"
#include "util_version.hpp"

namespace yh::sparsematrix
{
    const yh::util::Version k_VERSION(std::string("sparsematrix"),
                                      20230709,
                                      0,
                                      0,
                                      k_BUILD_DATE,
                                      k_BUILD_TIME);
}  // namespace yh::sparsematrix

#endif  // YH_SPARSEMATRIX_VERSION_H_
