#pragma once

#include <array>
#include <string>
#include <vector>
#include <unordered_map>

struct Option {
    std::string b1fname;
    std::string b2fname;
    std::string hfname;
    std::string prefix { "phasing" };

    int threads { 1 };

    int Check();
};

class Bam2Bed {
public:
    Bam2Bed(Option &opt): opt_(opt) {}
    void LoadBed();
    void ToBed();
    void Run();

private:
    Option opt_;
    std::unordered_map<std::string, std::array<std::string, 2>> map_;
};