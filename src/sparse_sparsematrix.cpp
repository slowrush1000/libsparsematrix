/**
 * @file    sparse_sparsematrix.cpp
 * @author  Cheon Younghoe
 */

#include "sparse_sparsematrix.hpp"
#include <iostream>
#include <vector>
#include <sstream>
#include <gzstream.h>
#include <fmt/printf.h>

sparse::Sparsematrix::Sparsematrix()
{
}

sparse::Sparsematrix::~Sparsematrix()
{
}

// ijv file : one indexed
// row col value
// i j value
void                        
sparse::Sparsematrix::read_ijv_file(const std::string& file_name, const bool zero_index, const bool sym)
{
    std::cout << "# read ijv file(" << file_name << ") start.\n";

    igzstream   file(file_name.c_str());
    if (false == file.good())
    {
        std::cout << "# error : ijv file(" << file_name << ") open failed.\n";
        return;
    }

    std::string         line;
    sparse::index_tt    row;
    sparse::index_tt    col;
    sparse::value_tt    value;
    std::size_t         n_lines = 0;
    const std::size_t   k_line_step = 100'000'000;

    while (std::getline(file, line))
    {
        ++n_lines;
        if (0 == (n_lines%k_line_step))
        {
            std::cout << "    " << n_lines << " lines\n";
        }

        if (1 == n_lines)
            continue;

        std::stringstream   ss(line);
        ss >> row >> col >> value;

        if (true == zero_index)
        {
            this->add_value(row, col, value);
            if (false == sym)
                this->add_value(col, row, value);
        }
        else
        {
            this->add_value(row - 1, col - 1, value);
            if (false == sym)
                this->add_value(col - 1, row - 1, value);
        }
    }
    std::cout << "    " << n_lines << " lines\n";

    file.close();

    std::cout << "# read ijv file(" << file_name << ") end.\n";
}

void                
sparse::Sparsematrix::add_value(const sparse::index_tt row, const sparse::index_tt col, const sparse::value_tt value)
{
    auto    found_row   = m_rows.find(row);
    if (m_rows.end() == found_row)
    {
        m_max_row   = std::max(m_max_row, row);

        sparse::Cols    cols;
        cols[col]       = value;
        m_rows[row]     = cols;
    }
    else
    {
        auto&   cols        = found_row->second;
        auto    found_col   = cols.find(col);

        m_max_col   = std::max(m_max_col, col);

        if (cols.end() == found_col)
        {
            cols[col]       = value;
        }
        else
        {

            auto&   value_t = found_col->second;
            value_t         += value;
        }
    }
}

std::optional<double>
sparse::Sparsematrix::get_value(const index_tt row, const index_tt col)
{
    auto    found_row   = m_rows.find(row);
    if (m_rows.end() == found_row)
    {
        return {};
    }
    else
    {
        auto    cols        = found_row->second;
        auto    found_col   = cols.find(col);
        if (cols.end() == found_col)
        {
            return {};
        }
        else
        {
            return found_col->second;
        }
    }
}

void
sparse::Sparsematrix::print()
{
    std::cout << "# print start.\n";

    if (m_max_row == m_max_col)
    {
        std::cout << "    row == row : " << m_max_row << " " << m_max_col << "\n";
    }
    else
    {
        std::cout << "    row != row : " << m_max_row << " " << m_max_col << "\n";
    }

    auto    msg = std::string("");

    for (auto row = 0; row < m_rows.size(); ++row)
    {
        //
        std::vector<sparse::index_tt>   col_keys;
        for (const auto& iter : m_rows[row])
        {
            col_keys.emplace_back(iter.first);
        }
        std::sort(col_keys.begin(), col_keys.end());

        //
        for (const auto& col : col_keys)
        {
            auto    value       = m_rows[row][col];

            if (true == this->get_sym())
            {
                if (row > col)
                {
                    continue;
                }
            }

            if (false == m_zero_index)
            {
                msg = fmt::sprintf("%lld %lld %e\n", row + 1, col + 1, value);
            }
            else
            {
                msg = fmt::sprintf("%lld %lld %e\n", row, col, value);
            }
            std::cout << msg;
        }
    }

    std::cout << "# print end.\n";
}

// oneapi mkl pardiso : csr
// zero index
void                        
sparse::Sparsematrix::make_csr_ia_ja_a(sparse::index_tt* n, sparse::index_tt* nnz, sparse::index_tt** ia, sparse::index_tt** ja, sparse::value_tt** a)
{
    std::cout << "# make ia/ja/a(mkl pardiso) start.\n";

    //
    *n      = m_rows.size();
    *nnz    = this->get_nnz();

    //
    sparse::index_tt*   t_ia    = new sparse::index_tt[*n + 1];
    sparse::index_tt*   t_ja    = new sparse::index_tt[*nnz];
    sparse::value_tt*   t_a     = new sparse::value_tt[*nnz];

    //
    sparse::index_tt    t_nnz   = 0;

    //
    for (auto row = 0; row < m_rows.size(); ++row)
    {
        t_ia[row]         = t_nnz;

        //
        std::vector<sparse::index_tt>   col_keys;
        for (const auto& iter : m_rows[row])
        {
            col_keys.emplace_back(iter.first);
        }
        std::sort(col_keys.begin(), col_keys.end());

        //
        for (const auto& col : col_keys)
        {
            auto    value       = m_rows[row][col];

            if (true == this->get_sym())
            {
                if (row > col)
                {
                    continue;
                }
            }

            t_ja[t_nnz]      = col;
            t_a[t_nnz]        = value;

            ++t_nnz;
        }
    }
    t_ia[m_rows.size()]   = t_nnz;

    *ia     = t_ia;
    *ja     = t_ja;
    *a      = t_a;

    std::cout << "# make ia/ja/a(mkl pardiso) end.\n";
}

// mumps : Centralized assembled matrix
// one index
void                        
sparse::Sparsematrix::make_cam_irn_jcn_a(sparse::index_tt* n, sparse::index_tt* nnz, sparse::index_tt** irn, sparse::index_tt** jcn, sparse::value_tt** a)
{
    std::cout << "# make irn/jcn/a(mumps) start.\n";

    //
    *n      = m_rows.size();
    *nnz    = this->get_nnz();

    //
    sparse::index_tt*    t_irn    = new sparse::index_tt[*nnz];
    sparse::index_tt*    t_jcn    = new sparse::index_tt[*nnz];
    sparse::value_tt*    t_a      = new sparse::value_tt[*nnz];

    //
    sparse::index_tt    t_nnz   = 0;

    //
    for (auto row = 0; row < m_rows.size(); ++row)
    {
        //
        std::vector<sparse::index_tt>   col_keys;
        for (const auto& iter : m_rows[row])
        {
            col_keys.emplace_back(iter.first);
        }
        std::sort(col_keys.begin(), col_keys.end());

        //
        for (const auto& col : col_keys)
        {
            auto    value       = m_rows[row][col];

            if (true == this->get_sym())
            {
                if (row > col)
                {
                    continue;
                }
            }

            t_irn[t_nnz]      = row + 1;
            t_jcn[t_nnz]      = col + 1;
            t_a[t_nnz]        = value;

            ++t_nnz;
        }
    }

    *irn    = t_irn;
    *jcn    = t_jcn;
    *a      = t_a;

    std::cout << "# make irn/jcn/a(mumps) end.\n";
}

sparse::index_tt                    
sparse::Sparsematrix::get_nnz()
{
    sparse::index_tt    nnz = 0;

    for (auto row = 0; row < m_rows.size(); ++row)
    {
        if (false == this->get_sym())
        {
            nnz = m_rows[row].size();
        }
        else
        {
            for (const auto& iter : m_rows[row])
            {
                auto    col = iter.first;
                if (row <= col)
                    ++nnz;
            }
        }
    }

    return nnz;
}

void                        
sparse::Sparsematrix::run(int targc, char* targv[])
{
    if (2 != targc)
    {
        std::cout << "sparsematrix usage :\n";
        std::cout << "    % sparsematrix.exe ijv_file(or tri file)\n";
        return;
    }

    this->read_ijv_file(std::string(targv[1]), false, true);

    this->print();

    //
    sparse::index_tt    n;
    sparse::index_tt    nnz;

    sparse::index_tt*   ia;
    sparse::index_tt*   ja;
    sparse::value_tt*   a;

    this->make_csr_ia_ja_a(&n, &nnz, &ia, &ja, &a);

    std::cout << "n     : " << n << "\n";
    std::cout << "nnz   : " << nnz << "\n";
    std::cout << "ia    :";
    for (auto i = 0; i < (n + 1); ++i)
    {
        std::cout << " " << ia[i];
    }
    std::cout << "\n";
    std::cout << "ja    :";
    for (auto i = 0; i < (nnz); ++i)
    {
        std::cout << " " << ja[i];
    }
    std::cout << "\n";
    std::cout << "a     :";
    for (auto i = 0; i < (nnz); ++i)
    {
        std::cout << " " << a[i];
    }
    std::cout << "\n";
    
    delete []   ia;
    delete []   ja;
    delete []   a;

    //
    sparse::index_tt*   irn;
    sparse::index_tt*   jcn;

    this->make_cam_irn_jcn_a(&n, &nnz, &irn, &jcn, &a);

    std::cout << "n     : " << n << "\n";
    std::cout << "nnz   : " << nnz << "\n";
    std::cout << "irn   :";
    for (auto i = 0; i < nnz; ++i)
    {
        std::cout << " " << irn[i];
    }
    std::cout << "\n";
    std::cout << "jcn   :";
    for (auto i = 0; i < nnz; ++i)
    {
        std::cout << " " << jcn[i];
    }
    std::cout << "\n";
    std::cout << "a     :";
    for (auto i = 0; i < nnz; ++i)
    {
        std::cout << " " << a[i];
    }
    std::cout << "\n";

    delete []   irn;
    delete []   jcn;
    delete []   a;
}
