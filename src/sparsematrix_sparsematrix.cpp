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

    // if matrix is symmetric, upper triangle matrix is generated.(For Intel Pardiso Solver)
    if ((true == m_symmetric) && (row_index > column_index))
        return;

    auto found = m_rows.find(row_index);
    if (m_rows.end() == found)
    {
        yh::sparsematrix::column_tt column;
        column[column_index] = value;
        m_rows[row_index]    = column;
    }
    else
    {
        auto &column = found->second;
        auto found_1 = column.find(column_index);
        if (column.end() == found_1)
        {
            column[column_index] = value;
        }
        else
        {
            auto &value_1 = found_1->second;
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
        auto &column = found->second;
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

    for (auto iter : m_rows)
        row_index.emplace_back(iter.first);

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
        const auto &column = found_row->second;

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
    yh::sparsematrix::index_tt nnz = 0;

    auto row_index = this->GetRowIndex();
    for (const auto &row : row_index)
    {
        auto column_index = this->GetColumnIndex(row);
        for (const auto &column : column_index)
            ++nnz;
    }

    return nnz;
}

void
yh::sparsematrix::SparseMatrix::ReadIJVFile(const std::string &filename)
{
    igzstream file;
    file.open(filename.c_str());
    if (false == file.good())
    {
        std::cout << "# error : file(" << filename << ") open failed.\n";
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
yh::sparsematrix::SparseMatrix::WriteIJVFile(const std::string &filename)
{
    std::ofstream file;
    file.open(filename);
    if (false == file.is_open())
    {
        std::cout << "# error : file(" << filename << ") open failed.\n";
        exit(0);
    }

    auto row_index = this->GetRowIndex();
    for (const auto &row : row_index)
    {
        auto column_index = this->GetColumnIndex(row);
        for (const auto &column : column_index)
            file << fmt::format(
                "{} {} {:e}\n", row, column, this->GetValue(row, column));
    }

    file.close();
}

std::tuple<long long*, long long*, long long*, void*> 
yh::sparsematrix::SparseMatrix::GetPardisoNIaJaA64()
{
    return this->GetPardisoNIaJaA<long long, double>();
    //const auto [n, a, ia, ja] = this->GetPardisoAIaJa<long long, double>();
}
