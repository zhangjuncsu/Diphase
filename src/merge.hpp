#pragma once

#include "bam_reader.hpp"
#include "read_store.hpp"

struct Options {
    std::string pfname;
    std::string bfname { "-" };
    std::string prefix { "phasing" };

    int threads { 1 };
    int threshold { 5 };

    int Check();
};

class Merge {
public:
    Merge(Options &opt): opt_(opt) {}

    enum { 
        START = 1, 
        END
    };

    std::vector<std::array<int64_t, 2>> GenerateSubctg(std::vector<int64_t> &poss, std::vector<char> &marks);
    std::unordered_map<std::string, std::string> MergeCtg(std::vector<bam1_t*> &records);
    void MergeCtgs();
    void Run();

private:
    Options opt_;
    ReadStore rs_;
    BamReader reader_;
};