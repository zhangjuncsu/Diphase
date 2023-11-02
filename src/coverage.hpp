#pragma once

#include "bam_reader.hpp"

#include <unordered_map>

struct Options {
    std::string bfname { "-" };
    std::string rfname;
    std::string ofname;
    std::string cname;
    std::string vfname;
    
    int threads { 1 };

    int Check();
};

class Coverage {
public:
    Coverage(Options &opt): opt_(opt) {}
    void Init();
    void GetCoverage();
    void LoadVariants();
    int FilterHiCPairWithSNP(std::array<bam1_t*, 2> &pair);
    void Run();

private:
    Options opt_;
    BamReader reader_;

    std::unordered_map<std::string, std::vector<int>> cov_;
    std::unordered_map<std::string, std::vector<long>> vars_;
};