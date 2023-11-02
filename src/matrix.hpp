#pragma once

#include <vector>
#include <string>

class Matrix {
public:
    Matrix(uint32_t row, uint32_t col): row_(row), col_(col) {
        data_.resize(row);
        for(auto &r: data_) r.resize(col);
    }
    Matrix(uint32_t size): row_(size), col_(size) {
        data_.resize(size);
        for(auto &r: data_) r.resize(size);
    }
    Matrix(): row_(0), col_(0) {}

    Matrix &operator=(const Matrix &rhs) {
        row_ = rhs.row_;
        col_ = rhs.col_;
        data_ = rhs.data_;
        return *this;
    }

    void Resize(uint32_t row, uint32_t col) {
        row_ = row;
        col_ = col;
        data_.resize(row);
        for(auto &r: data_) r.resize(col);
    }

    void Resize(uint32_t size) {
        row_ = size;
        col_ = size;
        data_.resize(size);
        for(auto &r: data_) r.resize(size);
    }

    void Dump(const std::string &fname) const;
    void Load(const std::string &fname);
    int AddLink(uint32_t row, uint32_t col);

    void Print() const;
    void Print(uint32_t row, uint32_t col) const;

    void Set(uint32_t row, uint32_t col, double val) { data_[row][col] = val; }
    void SetCol(uint32_t row, std::vector<double> &col) { data_[row] = col; }

    uint32_t RowSize() const { return row_; }
    uint32_t ColSize() const { return col_; }

    double Get(uint32_t row, uint32_t col) const { return data_[row][col]; }

    std::vector<std::vector<double>> GetMatrix() { return data_; }
    std::vector<std::vector<double>> GetMatrix() const { return data_; }

private:
    uint32_t row_;
    uint32_t col_;
    std::vector<std::vector<double>> data_;

    const char *magic_ = "PHASING";
    const int mlen = sizeof(magic_);
};