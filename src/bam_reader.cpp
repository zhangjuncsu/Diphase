#include "bam_reader.hpp"
#include "utility.hpp"
#include "logger.hpp"

#include <iostream>

BamReader::~BamReader() {
    bam_destroy1(record_);
    bam_hdr_destroy(header_);
    // sam_close(reader_);
    hts_close(hfile_);
}

void BamReader::Initialize(const std::string &fname) {
    fname_ = fname;
    hfile_ = hts_open(fname.c_str(), "r");
    const htsFormat *fmt = hts_get_format(hfile_);
    LOG(INFO)("Detect %s file type\n", hts_format_file_extension(fmt));
    // std::cerr << "[" << GetCurTime() << "] Detect " << hts_format_file_extension(fmt) << " file type\n";

    // if((reader_ = sam_open_format(fname.c_str(), "r", fmt)) == 0) {
    //     std::cerr << "[" << GetCurTime() << "] Could not open " << fname << " for reading\n";
    //     exit(EXIT_FAILURE);
    // }
    if((header_ = sam_hdr_read(hfile_)) == 0) {
        LOG(ERROR)("Failed to read the header from %s\n", fname.c_str());
        // std::cerr << "[" << GetCurTime() << "] Failed to read the header from " << fname << "\n";
        // exit(EXIT_FAILURE);
    }
    record_ = bam_init1();
    if(record_ == NULL) {
        LOG(ERROR)("Failed to init bam1_t structure\n");
        // std::cerr << "[" << GetCurTime() << "] Failed to init bam1_t structure\n";
        // exit(EXIT_FAILURE);
    }
}

std::vector<bam1_t*> BamReader::LoadAll() {
    std::vector<bam1_t*> rec;
    while(1) {
        int ret = sam_read1(hfile_, header_, record_);
        if(ret == -1) break;
        else if(ret < -1) {
            LOG(ERROR)("An error occur when reading %s\n", fname_.c_str());
            // std::cerr << "[" << GetCurTime() << "] An error occur when reading " << fname_ << "\n";
            // exit(EXIT_FAILURE);
        }
        rec.emplace_back(bam_dup1(record_));
    }
    return rec;
}

std::vector<bam1_t*> BamReader::Load(const std::string &ctg) {
    std::vector<bam1_t*> rec;
    htsFile *hfile = hts_open(fname_.c_str(), "r");
    hts_idx_t *bam_idx = sam_index_load(hfile, fname_.c_str());
    if(bam_idx == NULL) {
        LOG(ERROR)("Could not load .bai index file for %s\nPlease run \"samtools index %s\" before\n", fname_.c_str(), fname_.c_str());
        // std::cerr << "[" << GetCurTime() << "] Could not load .bai index file for " << fname_ << "\n";
        // std::cerr << "Please run \"samtools index " << fname_ << "\" brfore\n";
        // exit(EXIT_FAILURE);
    }
    int contig_id = bam_name2id(header_, ctg.c_str());
    bam1_t *record = bam_init1();
    hts_itr_t *iter = sam_itr_queryi(bam_idx, contig_id, 0, header_->target_len[contig_id]);
    int ret;
    while((ret = sam_itr_next(hfile, iter, record)) >= 0) {
        rec.emplace_back(bam_dup1(record));
    }
    sam_itr_destroy(iter);
    hts_idx_destroy(bam_idx);
    bam_destroy1(record);
    hts_close(hfile);
    return rec;
}

// 0 on success, -1 on EOF
int BamReader::LoadPair(std::array<bam1_t*, 2> &pair) {
    int ret1 = sam_read1(hfile_, header_, pair[0]);
    if(ret1 == -1) return -1;
    else if(ret1 < -1) {
        LOG(ERROR)("An error occur when reading %s\n", fname_.c_str());
        // std::cerr << "[" << GetCurTime() << "] An error occur when reading " << fname_ << "\n";
        // exit(EXIT_FAILURE);
    }
    int ret2 = sam_read1(hfile_, header_, pair[1]);
    if(ret2 < -1) {
        LOG(ERROR)("An error occur when reading %s\n", fname_.c_str());
        // std::cerr << "[" << GetCurTime() << "] An error occur when reading " << fname_ << "\n";
        // exit(EXIT_FAILURE);
    }
    if(ret1 > -1 && ret2 > -1) return 0;
    else return -1;
}

int BamReader::LoadOneRecord(bam1_t *record) {
    return sam_read1(hfile_, header_, record);
}

std::size_t BamReader::LoadBatchReads(std::vector<std::vector<bam1_t*>> &vec) {
    std::size_t index = 0;
    if(eof_) return 0;
    while(1) {
        if(!used_) {
            vec[index].emplace_back(bam_dup1(record_));
            used_ = true;
        }
        int r = sam_read1(hfile_, header_, record_);
        if(r < 0) {
            if(r == -1) eof_ = true;
            index++;
            break;
        }
        if(!vec[index].empty() && strcmp(bam_get_qname(vec[index].back()), bam_get_qname(record_)) != 0) index++;
        if(index == vec.size()) {
            used_ = false;
            break;
        }
        vec[index].emplace_back(bam_dup1(record_));
    }
    return index;
}

std::size_t BamReader::LoadBatchPair(std::vector<std::vector<bam1_t*>> &vp) {
    std::size_t index = 0;
    if(eof_) return 0;
    while(1) {
        if(!used_) {
            vp[index].emplace_back(bam_dup1(record_));
            used_ = true;
        }
        int r = sam_read1(hfile_, header_, record_);
        if(r < 0) {
            if(r == -1) eof_ = true;
            index++;
            break;
        }
        if(!vp[index].empty() && strcmp(bam_get_qname(vp[index].back()), bam_get_qname(record_)) != 0) index++;
        if(index == vp.size()) {
            used_ = false;
            break;
        }
        vp[index].emplace_back(bam_dup1(record_));
    }
    return index;
}

std::vector<std::string> BamReader::LoadRefName() const {
    return std::vector<std::string>(header_->target_name, header_->target_name + header_->n_targets);
}