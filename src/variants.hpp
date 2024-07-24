#pragma once

#include "read_store.hpp"
#include "bam_reader.hpp"
#include "utility.hpp"
#include "logger.hpp"

#include <vector>
#include <string>
#include <unordered_map>

struct Table {
    void SetRef(uint8_t b) { counts[0] = b; }
    void IncreaseM(uint8_t b) { counts[b + 1]++; }
    void IncreaseCov() { counts[5]++; }
    bool Valid() const { return var[0] != var[1]; }
    bool AtM(uint8_t b) const { return b == var[0] || b == var[1]; }

    void Confirm();
    int Offset(uint8_t b) const {
        return var[0] == b? 0: (var[1] == b? 1: -1);
    }
    uint16_t counts[6] = { 0 };
    uint8_t var[2] = { 0 };
};

struct Options {
    int min_length { 3000 };
    int min_align_length { 2000 };
    int overhang { 2000 };
    int mapq { 5 };
    int snp_match_length { 5 };
    int threads { 1 };

    double align_rate { 0.6 };
    double overhang_rate { 0.2 };
    double identity { 0.7 };

    bool check { false };
    bool write_all { false };

    std::string ctgfname;
    std::string bamfname;
    std::string vfname;
    std::string pcfname;
    std::string acfname;
    std::string mapfname;

    std::string prefix { "phasing" };

    int Check() const;
};

class Variants {
public:
    Variants(const ReadStore &ctg, const std::vector<bam1_t*> &record, 
            std::string &ctg_name, std::vector<long> &pos, Options &opt);
    
    uint8_t Convert(uint8_t b) const {
        if(b == 1) return 0;
        else if(b == 2) return 1;
        else if(b == 4) return 2;
        else if(b == 8) return 3;
        else if(b == 15) {
            LOG(WARNING)("There is a 'N' in sequence\n");
            // std::cerr << "[" << GetCurTime() << "] There is a 'N' in sequence\n";
            return 4;
        } else {
            LOG(WARNING)("Unrecognized character in sequence\n");
            // std::cerr << "[" << GetCurTime() << "] Unrecognized character in sequence\n";
            return 5;
        }
    }

    bool Valid(bam1_t* record, Options &opt);

    void FindSNPInContig(Options &opt);
    void ConfirmSNPs();
    void FindSNPInReads(Options &opt);

    void DumpSNP(std::ostream &out, bool dump_all) const;
    void DumpReadInfo(std::ostream &out) const;

private:
    std::string ctg_name_;

    std::vector<bam1_t*> bam_record_;
    std::vector<Table> snps_;
    std::vector<long> pos_;
    const ReadStore ctg_store_;
    std::unordered_map<std::string, std::vector<std::array<int, 4>>> read_info_;
};

class SNP {
    public:
    SNP(Options &opt): opt_(opt) {}

    void Load();
    std::unordered_map<std::string, std::vector<long>> GetVariantsWithClair3() const;
    bool Filtered(bam1_t *record) const;
    void Run();

private:
    Options opt_;
    ReadStore ctg_store_;
    BamReader reader_;
};