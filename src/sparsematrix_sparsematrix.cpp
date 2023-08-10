/**
 * @file sparsematrix_sparsematrix.cpp
 * @author Cheon Younghoe (you@dosparsematrix.com)
 * @brief
 * @version 0.1
 * @date 2023-07-08
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "sparsematrix_sparsematrix.hpp"

#include <fmt/core.h>
#include <gzstream.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

yh::sparsematrix::SparseMatrix::~SparseMatrix() {}

void
yh::sparsematrix::SparseMatrix::AddValue(
    const yh::sparsematrix::index_tt row_index,
    const yh::sparsematrix::index_tt column_index,
    const yh::sparsematrix::value_tt value)
{
    if (m_row_size < row_index) m_row_size = row_index;
    if (m_column_size < column_index) m_column_size = column_index;

    // if matrix is symmetric, upper triangle matrix is generated.(For Intel
    // Pardiso Solver)
    if ((true == m_symmetric) && (row_index > column_index)) return;

    auto found = m_rows.find(row_index);
    if (m_rows.end() == found)
    {
        yh::sparsematrix::column_tt column;
        column[column_index] = value;
        m_rows[row_index]    = column;
        ++m_nnz;
        /*
        std::cout << __func__ << "row column value nnz : " << row_index << " "
                  << column_index << " " << value << " " << m_nnz << "\n";
                  */
    }
    else
    {
        auto& column = found->second;
        auto found_1 = column.find(column_index);
        if (column.end() == found_1)
        {
            column[column_index] = value;
            ++m_nnz;
            /*
            std::cout << __func__ << "row column value nnz : " << row_index
                      << " " << column_index << " " << value << " " << m_nnz
                      << "\n";
                      */
        }
        else
        {
            auto& value_1 = found_1->second;
            value_1 += value;
        }
    }
}

yh::sparsematrix::value_tt
yh::sparsematrix::SparseMatrix::GetValue(
    const yh::sparsematrix::index_tt row_index,
    const yh::sparsematrix::index_tt column_index)
{
    auto found = m_rows.find(row_index);
    if (m_rows.end() == found)
    {
        return 0.0;
    }
    else
    {
        auto& column = found->second;
        auto found_1 = column.find(column_index);
        if (column.end() == found_1)
        {
            return 0.0;
        }
        else
        {
            return found_1->second;
        }
    }
}

void
yh::sparsematrix::SparseMatrix::SetOneIndex(const bool one_index)
{
    m_one_index = one_index;
}

bool
yh::sparsematrix::SparseMatrix::GetOneIndex() const
{
    return m_one_index;
}

void
yh::sparsematrix::SparseMatrix::SetSymmetric(const bool symmetric)
{
    m_symmetric = symmetric;
}

bool
yh::sparsematrix::SparseMatrix::GetSymmetric() const
{
    return m_symmetric;
}

std::vector<yh::sparsematrix::index_tt>
yh::sparsematrix::SparseMatrix::GetRowIndex()
{
    std::vector<yh::sparsematrix::index_tt> row_index;

    for (auto iter : m_rows) row_index.emplace_back(iter.first);

    std::sort(row_index.begin(), row_index.end());

    return row_index;
}

std::vector<yh::sparsematrix::index_tt>
yh::sparsematrix::SparseMatrix::GetColumnIndex(
    const yh::sparsematrix::index_tt row_index)
{
    std::vector<yh::sparsematrix::index_tt> column_index;

    auto found_row = m_rows.find(row_index);
    if (m_rows.end() != found_row)
    {
        const auto& column = found_row->second;

        for (auto iter : column) column_index.emplace_back(iter.first);
    }

    std::sort(column_index.begin(), column_index.end());

    return column_index;
}

yh::sparsematrix::index_tt
yh::sparsematrix::SparseMatrix::GetRowSize() const
{
    return m_row_size;
}

yh::sparsematrix::index_tt
yh::sparsematrix::SparseMatrix::GetColumnSize() const
{
    return m_column_size;
}

yh::sparsematrix::index_tt
yh::sparsematrix::SparseMatrix::GetNNZ()
{
    return m_nnz;
}

void
yh::sparsematrix::SparseMatrix::ReadIJVFile(const std::string& filename)
{
    igzstream file;
    file.open(filename.c_str());
    if (false == file.good())
    {
        std::cout << "# error : file(" << filename << ") open failed("
                  << __func__ << ").\n";
        exit(0);
    }

    std::string line;
    std::size_t n_lines = 0;

    while (std::getline(file, line))
    {
        ++n_lines;
        if (0 == (n_lines % 1'000'000))
            std::cout << "    " << n_lines << " lines\n";

        yh::sparsematrix::index_tt row;
        yh::sparsematrix::index_tt column;
        yh::sparsematrix::value_tt value;

        std::stringstream sstream(line);
        sstream >> row >> column >> value;

        this->AddValue(row, column, value);
    }
    std::cout << "    " << n_lines << " lines\n";

    file.close();
}

void
yh::sparsematrix::SparseMatrix::WriteIJVFile(const std::string& filename)
{
    std::ofstream file;
    file.open(filename);
    if (false == file.is_open())
    {
        std::cout << "# error : file(" << filename << ") open failed("
                  << __func__ << ").\n";
        exit(0);
    }

    auto row_index = this->GetRowIndex();
    for (const auto& row : row_index)
    {
        auto column_index = this->GetColumnIndex(row);
        for (const auto& column : column_index)
            file << fmt::format(
                "{} {} {:e}\n", row, column, this->GetValue(row, column));
    }

    file.close();
}

void
yh::sparsematrix::SparseMatrix::WriteBinIaJaAFile(const std::string& prefix)
{
    std::cout << "# write bin ia, ja, a file(prefix : " << prefix << ").\n";
    const auto [n, ia, ja, a] = this->GetPardiso64NIaJaA();

    auto ia_filename          = fmt::format("{}.{}.ia.bin", prefix, *n);
    long long nnz             = this->WriteBinIaFile(ia_filename, n, ia, ja, a);
    auto ja_filename          = fmt::format("{}.{}.ja.bin", prefix, *n);
    this->WriteBinJaFile(ja_filename, n, ia, ja, a, nnz);
    auto a_filename = fmt::format("{}.{}.a.bin", prefix, *n);
    this->WriteBinAFile(a_filename, n, ia, ja, a, nnz);
}

long long
yh::sparsematrix::SparseMatrix::WriteBinIaFile(const std::string& filename,
                                               long long* n,
                                               long long* ia,
                                               long long* ja,
                                               void* a)
{
    auto file = std::ofstream(filename, std::ios::binary);
    if (false == file.is_open())
    {
        std::cout << "# error : file(" << filename << ") open failed("
                  << __func__ << ").\n";
        exit(0);
    }

    file.write(reinterpret_cast<const char*>(ia), sizeof(long long) * (*n + 1));

    file.close();

    return ia[*n];
}

void
yh::sparsematrix::SparseMatrix::WriteBinJaFile(const std::string& filename,
                                               long long* n,
                                               long long* ia,
                                               long long* ja,
                                               void* a,
                                               const long long nnz)
{
    auto file = std::ofstream(filename, std::ios::binary);
    if (false == file.is_open())
    {
        std::cout << "# error : file(" << filename << ") open failed("
                  << __func__ << ").\n";
        exit(0);
    }

    file.write(reinterpret_cast<const char*>(ja), sizeof(long long) * nnz);

    file.close();
}

void
yh::sparsematrix::SparseMatrix::WriteBinAFile(const std::string& filename,
                                              long long* n,
                                              long long* ia,
                                              long long* ja,
                                              void* a,
                                              const long long nnz)
{
    auto file = std::ofstream(filename, std::ios::binary);
    if (false == file.is_open())
    {
        std::cout << "# error : file(" << filename << ") open failed("
                  << __func__ << ").\n";
        exit(0);
    }

    file.write(reinterpret_cast<const char*>(a), sizeof(double) * nnz);

    file.close();
}

void
yh::sparsematrix::SparseMatrix::ReadBinIaJaAFile(const long long n,
                                                 const std::string& ia_filename,
                                                 const std::string& ja_filename,
                                                 const std::string& a_filename)
{
    std::cout << "# read bin ia, ja, a file.\n";
    //
    long long* ia       = this->ReadBinIaFile(n, ia_filename);
    long long* ja       = this->ReadBinJaFile(n, ia[n], ja_filename);
    void* a             = this->ReadBinAFile(n, ia[n], a_filename);

    //
    m_nnz               = 0;

    //
    /*
    std::cout << "# add value to sparsermatrix\n";
    std::cout << "n : " << n << "\n";
    std::cout << "nnz : " << m_nnz << "\n";
    */
    long long value_pos = 0;
    for (auto row = 1; row < (n + 1); ++row)
    {
        auto column_start_pos = ia[row - 1];
        auto column_end_pos   = ia[row];

        for (auto column_pos = column_start_pos; column_pos < column_end_pos;
             ++column_pos)
        {
            auto column  = ja[column_pos];
            double value = ((double*)a)[value_pos];

            this->AddValue(row, column, value);

            /*
            std::cout << "row column value value_pos " << row << " " << column
                      << " " << value << " " << value_pos << "\n";
                      */

            ++value_pos;
        }
    }
}

long long*
yh::sparsematrix::SparseMatrix::ReadBinIaFile(const long long n,
                                              const std::string& ia_filename)
{
    long long* ia = new long long[n + 1];

    auto file     = std::ifstream(ia_filename, std::ios::binary);
    if (false == file.is_open())
    {
        std::cout << "# error : file(" << ia_filename << ") open failed("
                  << __func__ << ").\n";
        exit(0);
    }

    file.read(reinterpret_cast<char*>(ia), sizeof(long long) * (n + 1));

    /*
    for (auto pos = 0; pos < (n + 1); ++pos)
    {
        std::cout << "ia[" << pos << "] : " << ia[pos] << "\n";
    }
    */

    file.close();

    return ia;
}

long long*
yh::sparsematrix::SparseMatrix::ReadBinJaFile(const long long n,
                                              const long long nnz,
                                              const std::string& ja_filename)
{
    long long* ja = new long long[nnz];

    auto file     = std::ifstream(ja_filename, std::ios::binary);
    if (false == file.is_open())
    {
        std::cout << "# error : file(" << ja_filename << ") open failed("
                  << __func__ << ").\n";
        exit(0);
    }

    file.read(reinterpret_cast<char*>(ja), sizeof(long long) * nnz);

    /*
    for (auto pos = 0; pos < nnz; ++pos)
    {
        std::cout << "ja[" << pos << "] : " << ja[pos] << "\n";
    }
    */

    file.close();

    return ja;
}

void*
yh::sparsematrix::SparseMatrix::ReadBinAFile(const long long n,
                                             const long long nnz,
                                             const std::string& a_filename)
{
    double* a = new double[nnz];

    auto file = std::ifstream(a_filename, std::ios::binary);
    if (false == file.is_open())
    {
        std::cout << "# error : file(" << a_filename << ") open failed("
                  << __func__ << ").\n";
        exit(0);
    }

    file.read(reinterpret_cast<char*>(a), sizeof(double) * nnz);

    /*
    for (auto pos = 0; pos < nnz; ++pos)
    {
        std::cout << "a[" << pos << "] : " << (double)a[pos] << "\n";
    }
    */

    file.close();

    return a;
}

std::tuple<long long*, long long*, long long*, void*>
yh::sparsematrix::SparseMatrix::GetPardiso64NIaJaA()
{
    return this->GetPardisoNIaJaA<long long, double>();
    // const auto [n, a, ia, ja] = this->GetPardisoAIaJa<long long, double>();
}
