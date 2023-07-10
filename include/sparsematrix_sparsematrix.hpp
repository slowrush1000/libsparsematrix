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
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace yh::sparsematrix
{
    using value_tt = double;
    using index_tt = std::size_t;

    // CCS
    using column_tt =
        std::unordered_map<index_tt, value_tt>;  // column_index, value
    using row_tt =
        std::unordered_map<index_tt, column_tt>;  // row_index, column_tt

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
            void SetSymmetric(const bool symmetric);
            bool GetSymmetric() const;

            std::vector<index_tt> GetRowIndex();
            std::vector<index_tt> GetColumnIndex(const index_tt row_index);

            index_tt GetRowSize() const;
            index_tt GetColumnSize() const;
            index_tt GetNNZ();

            void ReadIJVFile(const std::string& filename);
            void WriteIJVFile(const std::string& filename);

            template <typename T1, typename T2>
            std::tuple<T1*, T1*, T1*, void*>
            GetPardisoNIaJaA()
            {
                auto nnz       = this->GetNNZ();

                T1* pardiso_n  = new T1;
                T1* pardiso_ia = new T1[m_row_size + 1];
                T1* pardiso_ja = new T1[nnz];
                T2* pardiso_a  = new T2[nnz];

                //
                *pardiso_n     = m_row_size;

                T1 a_count     = 0;
                T1 ia_count    = 0;

                pardiso_ia[0]  = 0;

                auto row_index = this->GetRowIndex();
                for (const auto& row : row_index)
                {
                    if (1 == row)
                        pardiso_ia[0] = 0;
                    else
                        pardiso_ia[row - 1] = a_count - pardiso_ia[row - 2];

                    auto column_index = this->GetColumnIndex(row);
                    for (const auto& column : column_index)
                    {
                        pardiso_a[a_count]  = this->GetValue(row, column);
                        pardiso_ja[a_count] = column;
                        ++a_count;
                    }
                }
                pardiso_ia[m_row_size] = nnz;

                return {pardiso_n, pardiso_ia, pardiso_ja, pardiso_a};
            }

            std::tuple<long long*, long long*, long long*, void*>
            GetPardiso64NIaJaA();

        private:
            row_tt m_rows;
            index_tt m_row_size    = 0;
            index_tt m_column_size = 0;
            index_tt m_nnz         = 0;
            bool m_one_index       = true;
            bool m_symmetric       = false;
    };

}  // namespace yh::sparsematrix

#endif  // YH_SPARSEMATRIX_SPARSEMATRIX_H_
