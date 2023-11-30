#pragma once

#include "file_io.hpp"
#include "read_store.hpp"

#include <unordered_map>

struct Option {
    std::vector<std::string> vfname;
    std::vector<std::string> rfname;
    std::vector<std::string> h1fname;
    std::vector<std::string> h2fname;

    std::string prefix_ { "phasing" };

    int k { 31 };
    int cnt_thres { 50 };
    int threads { 1 };

    int Check() const;
};

class HiC {
public:
    struct Kmer {
        Kmer() = default;
        Kmer(uint64_t hash, uint64_t pos): hash(hash), pos(pos) {}
        ~Kmer() = default;

        uint32_t Id() const { return static_cast<uint32_t>( pos >> 32 ); }
        uint32_t Pos() const { return static_cast<uint32_t>( pos >> 1 ); }
        bool Strand() const { return pos & 1; }
        static uint64_t SortByHash(const Kmer &kmer) { return kmer.hash; }

        uint64_t hash;
        uint64_t pos;
    };

    struct Index {
        Index() = default;
        ~Index() = default;
        uint32_t Find(uint64_t key, const uint64_t **dst) const;
        struct Hash {
            size_t operator()(uint64_t key) const {
                return std::hash<uint64_t>()(key >> 1);
            }
        };
        struct KeyEqual {
            bool operator()(uint64_t lhs, uint64_t rhs) const {
                return (lhs >> 1) == (rhs >> 1);
            }
        };
        std::vector<uint64_t> offset_;
        std::unordered_map<uint64_t, uint64_t, Hash, KeyEqual> index_;
    };

    struct Candidate {
        Candidate() = default;
        Candidate(uint64_t ref, uint64_t offset): ref(ref), offset(offset) {}
        ~Candidate() = default;

        uint64_t ref;
        uint64_t offset;
    };

    struct PeHit { 
        uint64_t s, e, id, len; 
    };

    HiC(const Option &opt): opt_(opt), index_(1ULL << std::min(14, opt.k * 2)), index_all_(1ULL << std::min(14, opt.k * 2)) {};
    
    void LoadVariants();

    uint64_t Hash(uint64_t key, uint64_t mask);
    std::vector<Kmer> Sketch(const ReadStore::Unit &unit);
    void BuildIndex();
    void BuildIndexSingle();
    void StatIndex();

    std::vector<Kmer> GetKmer(const ReadStore::Unit &unit);
    void BuildIndexAll();

    std::vector<Candidate> CompressMappedPos(std::vector<Candidate> &cand, uint64_t thres);
    template <typename RandomAccessIterator>
    uint64_t Collect(RandomAccessIterator first, RandomAccessIterator last);
    void InterprePos(Candidate &cand, uint64_t &rev, uint64_t &id, uint64_t &ref_p, uint64_t &self_p, uint64_t &exact_len, uint64_t *total_len);
    int ExactMatch(const std::string &query, std::size_t q_beg, const std::string &target, std::size_t t_beg, uint64_t rev, uint64_t dir);
    void GetLongestHit(const std::string &seq, std::size_t i, int strand, std::vector<Candidate> &cand, const uint64_t *list, int cnt, uint64_t &suffix);
    std::vector<Candidate> GetAlignment(const std::string &seq, std::size_t id);

    void Get53List(std::vector<Candidate> &cand1, std::vector<Candidate> &l5, std::vector<Candidate> &l3);
    inline void SetPos(std::vector<Candidate> &cand1, std::vector<Candidate> &cand2, uint64_t id, PeHit &ph);
    std::vector<PeHit> DedupHits(std::vector<PeHit> &hits);
    void GenerateContact();

    void Run();

    template <typename RandomAccessIterator, typename KeySort>
    static void RadixSort(RandomAccessIterator first, RandomAccessIterator last, uint8_t max_bits, KeySort ks);

private:
    Option opt_;

    std::vector<Index> index_;
    std::vector<Index> index_all_;
    ReadStore rs_ref_;

    std::vector<std::string> hic_name_;
    std::unordered_map<std::string, std::vector<std::size_t>> variants_;
};