#pragma once

#include "matrix.hpp"
#include "read_store.hpp"

#include <array>
#include <unordered_map>

struct Options {
    std::string pfname;
    std::string mfname;
    std::string cfname;
    std::string prefix { "phasing" };
 
    int threads { 1 };
    int iteration { 10000 };
    int seed { 10 };

    int Check();
};

class Phasing {
public:
    Phasing(Options &opt): opt_(opt) {}

    enum STATUS {
        INIT, WEAK, STRONG
    };

    struct Overlap {
        int first, second;
        int label;
        int status;
    };

    struct Group {
        int num_ctg;
        std::string group_name;
        std::vector<Overlap> ovlps;
    };

    std::vector<std::array<std::string, 2>> &LoadPair();
    void GroupCluster(std::vector<std::array<std::string, 2>> &pairs);
    void RandomPhasing();
    void Run();

private:
    std::vector<Group> group_;
    Matrix matrix_;
    ReadStore rs;
    Options opt_;
};