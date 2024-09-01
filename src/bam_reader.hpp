#pragma once

#include "htslib/hts.h"
#include "htslib/sam.h"

#include <array>
#include <string>
#include <vector>
#include <memory>

struct BamAllocator {
    using value_type = bam1_t*;

    template <typename T>
    struct rebind { using other = BamAllocator; };

    bam1_t** allocate(std::size_t n) {
        return static_cast<bam1_t**>(malloc(n * sizeof(bam1_t*)));
    }

    void deallocate(bam1_t** p, std::size_t n) {
        for (std::size_t i = 0; i < n; ++i) {
            bam_destroy1(p[i]);
        }
        free(p);
    }
};

struct BamDeleter {
    void operator()(bam1_t *record) {
        bam_destroy1(record);
    }
};

using BamVecs = std::vector<std::unique_ptr<bam1_t, BamDeleter>>;

class BamReader {
public:
    BamReader() = default;
    ~BamReader();
    void Initialize(const std::string &fname);

    const bam_hdr_t *Header() const { return header_; };
    const std::string GetCtgName(int id) const { return header_->target_name[id]; }

    std::vector<std::string> LoadRefName() const;
    std::vector<bam1_t*> LoadAll();
    int LoadPair(std::array<bam1_t*, 2> &pair);
    int LoadOneRecord(bam1_t* record);
    std::vector<bam1_t*> Load(const std::string &ctg);
    // std::size_t LoadBatchPair(std::vector<std::vector<bam1_t*>> &vp);
    std::size_t LoadBatchPair(std::vector<BamVecs> &vp);
    std::size_t LoadBatchReads(std::vector<std::vector<bam1_t*>> &vec);

private:
    std::string fname_;
    bool used_ { true };
    bool eof_ { false };
    bam1_t *record_;
    // samFile *reader_{ NULL };
    htsFile *hfile_{ NULL };
    bam_hdr_t *header_{ NULL };
};