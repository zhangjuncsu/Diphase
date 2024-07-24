#pragma once

#include <string>
#include <vector>
#include <memory>
#include <cassert>
#include <algorithm>

#include <stdio.h>

class Seq {
public:
    static std::string& ReverseComplement(const std::string &seq);
    static void ReverseComplementInPlace(std::string &seq);
};

class SerialDNATable {
public:
    SerialDNATable() {
        std::fill(table_, table_ + 256, -1);
        table_['a'] = table_['A'] = 0;
        table_['c'] = table_['C'] = 1;
        table_['g'] = table_['G'] = 2;
        table_['t'] = table_['T'] = 3;

        table_['n'] = table_['N'] = 4;
        table_['-'] = 5;
    }

    uint8_t operator[](char c) const { return table_[(uint8_t)c]; }

private:
    uint8_t table_[256];
};

class ComplementDNATable {
public:
    ComplementDNATable() {
        std::fill(table_, table_ + 256, -1);
        table_['a'] = table_['A'] = 0;
        table_['c'] = table_['C'] = 1;
        table_['g'] = table_['G'] = 2;
        table_['t'] = table_['T'] = 3;

        table_['n'] = table_['N'] = 4;
        table_['-'] = 5;

        std::swap(table_['a'], table_['t']);
        std::swap(table_['c'], table_['g']);
        // std::swap(table_['g'], table_['c']);
        // std::swap(table_['t'], table_['a']);
        std::swap(table_['A'], table_['T']);
        std::swap(table_['C'], table_['G']);
        // std::swap(table_['G'], table_['C']);
        // std::swap(table_['T'], table_['A']);
    }

    char operator[](char c) const { return table_[(uint8_t)c]; }

private:
    uint8_t table_[256];
};

class DNASeq {
public:
    DNASeq(): length_(0) {}
    DNASeq(const std::string &seq): length_(seq.size()), data_((seq.size() + 3) / 4) { Reset(seq); }

    ~DNASeq() = default;

    void Reset(const std::string &seq) {
        length_ = seq.size();
        data_.assign((seq.size() + 3) / 4, 0);
        for(std::size_t i = 0; i < seq.size(); ++i) {
            data_[i >> 2] |= (serial_dna_table_[seq[i]] << ((i & 3) << 1));
        }
    }

    std::size_t Size() const { return length_; }
    
    uint8_t operator[](std::size_t i) const {
        assert(i < length_);
        return (data_[i >> 2] >> ((i & 3) << 1)) & 0x3;
    }

    DNASeq operator+(const DNASeq &seq) const {
        return DNASeq(*ToString() + *seq.ToString());
    }

    bool operator<(const DNASeq &seq) const {
        for(std::size_t i = 0; i < length_; ++i) {
            if(i == seq.Size()) return true;
            else if(this->operator[](i) != seq[i]) return this->operator[](i) < seq[i];
        }
        return false;
    }

    std::shared_ptr<std::string> ToString(bool upper = true) const {
        std::shared_ptr<std::string> str(new std::string(length_, 'N'));
        const char *base = upper? "ACGT": "acgt";
        for(std::size_t i = 0; i < length_; ++i) {
            str->operator[](i) = base[(*this)[i]];
        }
        return str;
    }

    std::vector<uint8_t> &ToUint8(int s = 0, int e = -1, bool rc = false) const {
        assert(s >= 0 && (std::size_t)s <= length_);
        assert(e < 0 || (e >= 0 && e > s && (std::size_t)e <= length_));
        if(e < 0) e = length_;

        std::vector<uint8_t> seq(e - s, -1);
        if(rc) {
            for(std::size_t i = 0; i < seq.size(); ++i) {
                seq[i] = 3 - this->operator[](s + seq.size() - 1 - i);
            }
        } else {
            for(std::size_t i = 0; i < seq.size(); ++i) {
                seq[i] = this->operator[](s + i);
            }
        }
        return seq;
    }

    DNASeq& ReverseComplement() const {
        DNASeq seq;
        seq.Reset(Seq::ReverseComplement(*(this->ToString())));
        return seq;
    }

    const DNASeq& SubSeq(std::size_t pos, std::size_t size, bool rc = false) const {
        assert(pos <= length_ && pos >= 0);
        if(size > length_) size = length_;
        if(pos + size > length_) size = length_ - pos;
        DNASeq seq;
        if(rc) seq.Reset(Seq::ReverseComplement((*(this->ToString())).substr(pos, size)));
        else seq.Reset((*(this->ToString())).substr(pos, size));
        return seq;
    }

    const DNASeq SubSeq(std::size_t pos, bool rc = false) const {
        return SubSeq(pos, length_, rc);
    }

    bool operator==(const DNASeq &seq) const {
        if(length_ != seq.length_) return false;
        for(std::size_t i = 0; i < length_; ++i) {
            if(this->operator[](i) != seq[i]) return false;
        }
        return true;
    }

    static uint8_t Serial(char c) { return serial_dna_table_[c]; }
    static char Complement(char c) { return complement_dna_table_[c]; }
    static bool Cheak(const std::string &seq, bool n = true) {
        const static SerialDNATable table;
        if(n) {
            return std::find_if(seq.begin(), seq.end(), [](char c) -> bool {
                return table[c] == -1;
            }) == seq.end();
        } else {
            return std::find_if(seq.begin(), seq.end(), [](char c) -> bool {
                return table[c] < 4;
            }) == seq.end();
        }
    }

    static SerialDNATable serial_dna_table_;
    static ComplementDNATable complement_dna_table_;

private:
    std::size_t length_;
    std::vector<uint8_t> data_;
};