/**
 * @file    sparse_sparsematrix.hpp
 * @author  Cheon Younghoe
 */

#ifndef SPARSE_SPARSEMATRIX_H
#define SPARSE_SPARSEMATRIX_H

#include <unordered_map>
#include <optional>
#include <string>
#include <limits>

namespace sparse
{
    using index_tt                  = std::size_t;
    const index_tt  k_index_tt_init = 0;
    const index_tt  k_index_tt_max  = std::numeric_limits<index_tt>::max();
    const index_tt  k_index_tt_min  = std::numeric_limits<index_tt>::lowest();

    using   value_tt = double;
    using   Cols     = std::unordered_map<index_tt, double>;  // col, value
    using   Rows     = std::unordered_map<index_tt, Cols>;    // row, <col, value>
    
    class Sparsematrix
    {
        public:
            Sparsematrix();
            virtual ~Sparsematrix();

            void                        read_ijv_file(const std::string& file_name, const bool zero_index = false, const bool sym = true);

            void                        add_value(const index_tt row, const index_tt col, const value_tt value);
            std::optional<value_tt>     get_value(const index_tt row, const index_tt col);
            void                        print();
            void                        make_csr_ia_ja_a(index_tt* n, index_tt* nnz, index_tt** ia, index_tt** ja, value_tt** a);
            void                        make_cam_irn_jcn_a(index_tt* n, index_tt* nnz, index_tt** irn, index_tt** jcn, value_tt** a);
            index_tt                    get_nnz();

            void                        run(int targc, char* targv[]);

            void                set_zero_index(const bool zero_index)
            {
                m_zero_index    = zero_index;
            }
            bool                get_zero_index() const
            {
                return m_zero_index;
            }
            void                set_sym(const bool sym)
            {
                m_sym    = sym;
            }
            bool                get_sym() const
            {
                return m_sym;
            }
            void                set_max_row(const bool max_row)
            {
                m_max_row    = max_row;
            }
            index_tt            get_max_row() const
            {
                return m_max_row;
            }
            void                set_max_col(const bool max_col)
            {
                m_max_col    = max_col;
            }
            index_tt            get_max_col() const
            {
                return m_max_col;
            }

        private:
            Rows                m_rows;                     // zero-index based
            bool                m_zero_index    = true;
            bool                m_sym           = true;

            index_tt            m_max_row       = k_index_tt_init;
            index_tt            m_max_col       = k_index_tt_init;
    };
}

#endif  // SPARSE_SPARSEMATRIX_H
