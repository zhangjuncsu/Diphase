#pragma once

#include "htslib/hts.h"
#include "htslib/sam.h"

#include <array>
#include <string>
#include <vector>

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
    std::size_t LoadBatchPair(std::vector<std::vector<bam1_t*>> &vp);
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