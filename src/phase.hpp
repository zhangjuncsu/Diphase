#pragma once

#include "matrix.hpp"
#include "utility.hpp"
#include "read_store.hpp"
#include "bam_reader.hpp"

#include <unordered_map>
#include <unordered_set>

struct Options {
    std::string bfname; // alt2pri.sorted.bam
    std::string hfname; // hic2minced.bam
    std::string vfname; // var.txt
    std::string cfname; // mince.fasta
    std::string pfname; // pair.txt
    std::string mfname; // matrix.binmat
    std::string sfname; // clusters.by_name.txt
    std::string b1fname; // {prefix}.phased0.bed
    std::string b2fname; // {prefix}.phased1.bed
    std::string pcfname; // primary clair3
    std::string acfname; // alt clair3
    std::string prefix { "phasing" };

    int mapq { 10 };
    int threads { 1 };
    int iteration { 10000 };
    int seed { 10000 };

    bool dump_var { false };
    bool dump_filtered { false };
    bool dump_mat { false };
    bool filtered { false };
    bool print_mat { false };
    bool print_score { false };
    bool scaffold { false };
    bool porec { false };

    int Check();
};

class Phase {
public:
    Phase(Options& opt): opt_(opt) {}

    struct Overlap {
        uint32_t first, second;
        short label, label_tmp;
        int count1, count2;
    };

    struct Ctg {
        int num;
        std::string name;
        std::vector<Overlap> olps;
    };

    void GenerateVariantsWithClair3();
    void GenerateVariants();
    void LoadVariants();

    std::unordered_map<std::string, std::string> LoadBed() const;
    void GenerateMatrixScaffold(std::unordered_map<std::string, std::string> &bed);

    int FilterHiCPair(std::array<bam1_t*, 2> &pair);
    int FilterHiCPairWithSNP(std::array<bam1_t*, 2> &pair);
    void GenerateMatrix();
    void FilterAndGenerateMatrix();

    std::vector<std::size_t> FilterPoreC(std::vector<bam1_t*> &vec);
    std::vector<std::size_t> FilterPoreCWithSNP(std::vector<bam1_t*> &vec);
    void GenerateMatrixPoreC();
    void FilterAndGenerateMatrixPoreC();
    void IncreaseCov(std::string &name, bam1_t *a, bam1_t *b);

    // void IncreaseCov(std::string &name, std::array<bam1_t*, 2> &pair);
    std::unordered_map<std::string, std::size_t> DetectSwitchError();
    void FilterSwitchError(std::unordered_map<std::string, std::size_t> &se, std::vector<std::array<std::string, 2>> &pair);
    void FixPhase(std::string &bfname);

    void FixPhasePoreC(std::string &bfname);

    std::vector<std::array<std::string, 2>> LoadPair();
    void ClusterCtg(std::vector<std::array<std::string, 2>> &pair);
    void ClusterCtgScaffold(std::vector<std::array<std::string, 2>> &pair);

    int GetMaxAll(std::vector<Overlap> &olp) const;
    int GetMaxDiag(std::vector<Overlap> &olp) const;
    int GetMaxNoDiag(std::vector<Overlap> &olp, std::unordered_set<uint32_t> &filt) const;
    int GetMaxAllCount(std::vector<Overlap> &olp) const;

    void LocalMaximumAll(std::vector<Overlap> &olp);
    void LocalMaximumDiag(std::vector<Overlap> &olp);
    void LocalMaximumNoDiag(std::vector<Overlap> &olp, std::unordered_set<uint32_t> &filt);
    void LocalMaximumAllCount(std::vector<Overlap> &olp);

    void Print(std::vector<Overlap> &olp);
    void PrintScore(std::vector<Overlap> &olp, std::ofstream &out);

    void Phasing();
    void Run();
    void Run1();

private:
    Options opt_;
    ReadStore rs_;
    Matrix matrix_;
    BamReader reader_;
    std::vector<Ctg> contigs_;
    std::vector<std::string> names_;
    std::unordered_map<std::string, std::size_t> name2id_;
    std::unordered_map<std::string, std::vector<long>> vars_;
    std::unordered_map<std::string, std::vector<std::size_t>> cov_;
    std::unordered_map<std::string, std::size_t> switch_error_;
    std::unordered_map<std::string, char> consistency_;
};