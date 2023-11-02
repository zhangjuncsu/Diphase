#include "matrix.hpp"
#include "utility.hpp"

#include <cassert>
#include <cstring>
#include <iostream>

void Matrix::Dump(const std::string &fname) const {
    FILE *fw = fopen(fname.c_str(), "wb");
    if(fw == NULL) return;
    fwrite(magic_, 1, mlen, fw);
    fwrite(&row_, sizeof(uint32_t), 1, fw);
    fwrite(&col_, sizeof(uint32_t), 1, fw);
    for(uint32_t row = 0; row < row_; ++row) {
        fwrite(data_[row].data(), sizeof(double), col_, fw);
    }
    fwrite(magic_, 1, mlen, fw);
    fclose(fw);
}

void Matrix::Load(const std::string &fname) {
    char magic[mlen];
    FILE *fr = fopen(fname.c_str(), "rb");
    if(fr == NULL) {
        std::cerr << "[" << GetCurTime() << "] Could not open " << fname << " for reading\n";
        exit(EXIT_FAILURE);
    }
    fread(magic, 1, mlen, fr);
    if(strncmp(magic, magic_, mlen) != 0) {
        std::cerr << "[" << GetCurTime() << "] File is not correct\n";
        exit(EXIT_FAILURE);
    }
    fread(&row_, sizeof(uint32_t), 1, fr);
    fread(&col_, sizeof(uint32_t), 1, fr);
    data_.resize(row_);
    std::vector<double> column(col_);
    for(uint32_t row = 0; row < row_; ++row) {
        fread(&column[0], sizeof(double), col_, fr);
        SetCol(row, column);
    }
    char magic_end[mlen];
    fread(&magic_end, 1, mlen, fr);
    if(strncmp(magic_end, magic_, mlen) != 0) {
        std::cerr << "[" << GetCurTime() << "] File is not correct\n";
        exit(EXIT_FAILURE);
    }
    fclose(fr);
}

int Matrix::AddLink(uint32_t row, uint32_t col) {
    if(row >= row_ || col >= col_) return 0;
    data_[row][col] += 1;
    return 1;
}

void Matrix::Print() const {
    std::cerr << "[" << GetCurTime() << "] Matrix size: " << row_ << " X " << col_ << "\n";
    for(auto row: data_) {
        for(auto e: row) {
            std::cout << e << " ";
        }
        std::cout << "\n";
    }
}

void Matrix::Print(uint32_t row, uint32_t col) const {
    std::cerr << "[" << GetCurTime() << "] " << row << " X " << col << "\n";
    for(uint32_t r = 0; r < row; ++r) {
        for(uint32_t c = 0; c < col; ++c) {
            std::cout << data_[r][c] << " ";
        }
        std::cout << "\n";
    }
}