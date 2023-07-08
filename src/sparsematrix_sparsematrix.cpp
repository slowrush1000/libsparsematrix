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
#include <fstream>
#include <iostream>
#include <algorithm>

yh::sparsematrix::SparseMatrix::~SparseMatrix() {}

void
yh::sparsematrix::SparseMatrix::AddValue(
    const yh::sparsematrix::index_tt row_index,
    const yh::sparsematrix::index_tt column_index,
    const yh::sparsematrix::value_tt value)
{
    auto found_row = m_row.find(row_index);
    if (m_row.end() == found_row)
    {
        yh::sparsematrix::column_tt column;
        column[column_index] = value;
        m_row[row_index]     = column;
    }
    else
    {
        auto& row         = found_row->second;
        auto found_column = row.find(column_index);
        if (row.end() == found_column)
        {
            yh::sparsematrix::column_tt column;
            column[column_index] = value;
            m_row[row_index]     = column;
        }
        else
        {
            auto& t_value = found_column->second;
            t_value += value;
        }
    }
}

yh::sparsematrix::value_tt
yh::sparsematrix::SparseMatrix::GetValue(
    const yh::sparsematrix::index_tt row_index,
    const yh::sparsematrix::index_tt column_index)
{
    auto found_row = m_row.find(row_index);
    if (m_row.end() == found_row)
    {
        return 0.0;
    }
    else
    {
        auto& row         = found_row->second;
        auto found_column = row.find(column_index);
        if (row.end() == found_column)
            return 0.0;
        else
            return found_column->second;
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

std::vector<yh::sparsematrix::index_tt> 
yh::sparsematrix::SparseMatrix::GetRowIndex()
{
    std::vector<yh::sparsematrix::index_tt>     row_index;

    for (auto iter : m_row)
        row_index.emplace_back(iter.first);

    std::sort(row_index.begin(), row_index.end());

    return row_index;
}

std::vector<yh::sparsematrix::index_tt> 
yh::sparsematrix::SparseMatrix::GetColumnIndex(const yh::sparsematrix::index_tt row_index)
{
    std::vector<yh::sparsematrix::index_tt>     column_index;

    for (auto iter : m_row)
        row_index.emplace_back(iter.first);

    std::sort(row_index.begin(), row_index.end());

    return row_index;
}

void 
yh::sparsematrix::SparseMatrix::WriteIJVFile(const std::string& filename)
{
    std::fstream   file;
    file.open(filename);
    if (false == file.is_open())
    {
        std::cout << "# error : file(" << filename << ") open failed.\n";
        exit(0);
    }


    file.close();
}
