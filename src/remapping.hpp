#pragma once

#include "read_store.hpp"

struct Options {
    std::string paffname;
    std::string prifname;
    std::string altfname;
    std::string breakfname;
    std::string prefix { "phasing" };

    int threads { 1 };
    int block { 10000 };
    
    bool breaks { true };
    
    int Check();
};

class Remapping {
public:
    Remapping(const Options &opt): opt_(opt) {}

    std::unordered_map<std::string, std::vector<std::size_t>> LoadBps();

    void Remap();
    void Cut();
    void DumpPaf();
    void LoadPaf();
    void Run();

    void Load() {
        pstore_.Load(opt_.prifname);
        astore_.Load(opt_.altfname);
    }

    bool Empty() { return paf_.empty(); }

    bool Nested(int as, int ae, int bs, int be) {
        if(as > bs && ae < be) return true;
        return false;
    }

    bool Duplicated(int as, int ae, int bs, int be) {
        if(as == bs && ae == be) return true;
        return false;
    }

private:
    ReadStore pstore_;
    ReadStore astore_;
    Options opt_;

    std::unordered_map<std::string, std::vector<std::array<std::string, 13>>> paf_;
};