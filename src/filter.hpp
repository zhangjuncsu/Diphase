#pragma once

#include "bam_reader.hpp"

#include <string>
#include <unordered_map>

struct Options {
    std::string varfname;
    std::string infname { "-" };
    std::string outfname { "-" };

    int threads { 1 };
    int mapq { 0 };
    int clip { 1000 };

    bool sameid { false };
    
    int Check();
};

class Filter {
public:
    Filter(Options &opt): opt_(opt) {}
    bool FilterVariants(bam1_t *record);
    int FilterVariants(bam1_t *a, bam1_t *b);
    void LoadVariants();
    void Run();
    void Run1();
    void Run2();

private:
    Options opt_;
    BamReader reader_;

    std::unordered_map<std::string, std::vector<int64_t>> vars_;
};