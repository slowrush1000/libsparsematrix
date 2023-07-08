/**
 * @file sparsematrix_sparsematrix.hpp
 * @author Cheon Younghoe (you@dosparsematrix.com)
 * @brief
 * @version 0.1
 * @date 2023-07-08
 *
 * @copyright Copyright (c) 2023
 *
 */

#ifndef YH_SPARSEMATRIX_SPARSEMATRIX_H_
#define YH_SPARSEMATRIX_SPARSEMATRIX_H_

#include <cstdint>
#include <unordered_map>
#include <string>
#include <vector>

namespace yh::sparsematrix
{
    using value_tt  = double;
    using index_tt  = std::size_t;

    // CCS
    using column_tt = std::unordered_map<index_tt, value_tt>;
    using row_tt    = std::unordered_map<index_tt, column_tt>;

    class SparseMatrix
    {
        public:
            SparseMatrix() = default;
            ~SparseMatrix();

            void AddValue(const index_tt row_index,
                          const index_tt column_index,
                          const value_tt value);
            value_tt GetValue(const index_tt row_index,
                              const index_tt column_index);
            void SetOneIndex(const bool one_index);
            bool GetOneIndex() const;

            std::vector<index_tt> GetRowIndex();
            std::vector<index_tt> GetColumnIndex(const index_tt row_index);

            void WriteIJVFile(const std::string& filename);

        private:
            row_tt m_row;
            bool m_one_index = true;
    };

}  // namespace yh::sparsematrix

#endif  // YH_SPARSEMATRIX_SPARSEMATRIX_H_
