#include "variants.hpp"
#include "utility.hpp"
#include "file_io.hpp"

#include "getopt.h"

#include <regex>
#include <mutex>
#include <atomic>
#include <iostream>

void USAGE(Options &opt) {
    std::cerr << "Usage:\tSNPCall [options] -a alignment -c contigs\n\n";
    std::cerr << "\tOptions:\n";
    std::cerr << "\t\t-a | --align             <FILE>   mapped file name of reads against alt contig, sam or bam format\n";
    std::cerr << "\t\t-c | --contig            <FILE>   alt contig file name\n";
    std::cerr << "\t\t-v | --vars              <FILE>   variants file name\n";
    std::cerr << "\t\t-i | --identity          <FLOAT>  minimal mapping identity to filter [" << opt.identity << "]\n";
    std::cerr << "\t\t-l | --min_length        <INT>    minimal read length to filter [" << opt.min_length << "]\n";
    std::cerr << "\t\t-m | --min_align_length  <INT>    minimal mapping length to filter [" << opt.min_align_length << "]\n";
    std::cerr << "\t\t-o | --overhang          <INT>    max overhang to filter [" << opt.overhang << "]\n";
    std::cerr << "\t\t-p | --prefix            <STR>    prefix of output files [" << opt.overhang << "]\n";
    std::cerr << "\t\t-e | --overhang_rate     <FLOAT>  max overhang rate to filter [" << opt.overhang_rate << "]\n";
    std::cerr << "\t\t-q | --mapq              <INT>    minimal map quality to filter [" << opt.mapq << "]\n";
    std::cerr << "\t\t-f | --align_rate        <FLOAT>  minimal mapping rate to filter [" << opt.align_rate << "]\n";
    std::cerr << "\t\t-s | --snp_match         <INT>    number of matched base around snp (--check) [" << opt.snp_match_length << "]\n";
    std::cerr << "\t\t-t | --threads           <INT>    number of threads to use " << opt.threads << "]\n";
    std::cerr << "\t\t-h | --help                       display this message\n";
    std::cerr << "\t\t     --mapping           <FILE>   mapping file of alt to primary in paf\n";
    std::cerr << "\t\t     --priclair          <FILE>   primary clair file name\n";
    std::cerr << "\t\t     --altclair          <FILE>   alt clair file name\n";
    std::cerr << "\t\t     --check                      check the mapping around snp\n";
    std::cerr << "\t\t     --write_all                  write all base infomation\n";
}

int ParseArgument(int argc, char **argv, Options &opt) {
    const char *short_opt = "a:e:c:f:i:l:m:o:p:q:s:t:v:h";
    struct option long_opt[] = {
        {"align",               required_argument,  NULL,   'a'}, 
        {"contig",              required_argument,  NULL,   'c'}, 
        {"vars",                optional_argument,  NULL,   'v'}, 
        {"align_rate",          optional_argument,  NULL,   'f'}, 
        {"identity",            optional_argument,  NULL,   'i'}, 
        {"min_length",          optional_argument,  NULL,   'l'}, 
        {"min_align_length",    optional_argument,  NULL,   'm'}, 
        {"overhang",            optional_argument,  NULL,   'o'}, 
        {"overhang_rate",       optional_argument,  NULL,   'e'}, 
        {"prefix",              optional_argument,  NULL,   'p'}, 
        {"mapq",                optional_argument,  NULL,   'q'}, 
        {"snp_match",           optional_argument,  NULL,   's'}, 
        {"threads",             optional_argument,  NULL,   't'}, 
        {"help",                no_argument,        NULL,   'h'}, 
        {"check_snp",           no_argument,        NULL,   300}, 
        {"write_all",           no_argument,        NULL,   301}, 
        {"priclair",            required_argument,  NULL,   302}, 
        {"altclair",            required_argument,  NULL,   303}, 
        {"mapping",             required_argument,  NULL,   304}, 
        {0, 0,  0,  0}
    };
    int c;
    while((c = getopt_long(argc, argv, short_opt, long_opt, NULL)) != -1) {
        if(c == 'f') opt.align_rate = atof(optarg);
        else if(c == 'i') opt.identity = atof(optarg);
        else if(c == 'l') opt.min_length = atoi(optarg);
        else if(c == 'm') opt.min_align_length = atoi(optarg);
        else if(c == 'o') opt.overhang = atoi(optarg);
        else if(c == 'e') opt.overhang_rate = atof(optarg);
        else if(c == 'p') opt.prefix = optarg;
        else if(c == 'q') opt.mapq = atoi(optarg);
        else if(c == 's') opt.snp_match_length = atoi(optarg);
        else if(c == 't') opt.threads = atoi(optarg);
        else if(c == 'a') opt.bamfname = optarg;
        else if(c == 'c') opt.ctgfname = optarg;
        else if(c == 'v') opt.vfname = optarg;
        else if(c == 'h') {USAGE(opt); exit(EXIT_SUCCESS); }
        else if(c == 300) opt.check = true;
        else if(c == 301) opt.write_all = true;
        else if(c == 302) opt.pcfname = optarg;
        else if(c == 303) opt.acfname = optarg;
        else if(c == 304) opt.mapfname = optarg;
        else if(c == ':') {
            LOG(WARNING)("missing option argument in %c\n", optopt);
            // std::cerr << "[" << GetCurTime() << "] ERROR missing option argument in " << optopt << "\n";
            return 0;
        } else if(c == '?') {
            LOG(WARNING)("unknown option in %c\n", optopt);
            // std::cerr << "[" << GetCurTime() << "] ERROR unknown option in " << optopt << "\n";
            return 0;
        }
    }
    return opt.Check();
}

void Table::Confirm() {
    int total = 0;
    for(int i = 1; i < 5; ++i) total += counts[i];
    std::vector<int> ok;
    for(int i = 1; i < 5; ++i) {
        if(counts[i] > total * 0.8) break;
        if(counts[i] > total * 0.2) ok.emplace_back(i);
    }
    if(ok.size() >= 2) {
        std::sort(ok.begin(), ok.end(), [this](int a, int b) { return counts[a] > counts[b]; });
        var[0] = ok[0] - 1;
        var[1] = ok[1] - 1;
    }
}

int Options::Check() const {
    if(ctgfname.empty()) {
        LOG(WARNING)("Please specify the alt contig file name [-c | --contig]\n");
        // std::cerr << "[" << GetCurTime() << "] Please specify the alt contig file name [-c | --contig]\n";
        return 0;
    }
    if(bamfname.empty()) {
        LOG(WARNING)("Please specify the mapped file name of reads against alt contig [-a | --align]\n");
        // std::cerr << "[" << GetCurTime() << "] Please specify the mapped file name of reads against alt contig [-a | --align]\n";
        return 0;
    }
    if(threads < 1) {
        LOG(WARNING)("Threads must be >= 1 [-t | --threads]\n");
        // std::cerr << "[" << GetCurTime() << "] Threads must be >= 1 [-t | --threads]\n";
        return 0;
    }
    if(min_length < 0) {
        LOG(WARNING)("min_length must be >= 0 [-l | --min_length]\n");
        // std::cerr << "[" << GetCurTime() << "] min_length must be >= 0 [-l | --min_length]\n";
        return 0;
    }
    if(min_align_length < 0) {
        LOG(WARNING)("min_align_length should be >= 0 [-m | --min_align_length]\n");
        // std::cerr << "[" << GetCurTime() << "] min_align_length should be >= 0 [-m | --min_align_length]\n";
        return 0;
    }
    if(overhang < 0) {
        LOG(WARNING)("overhang should be >= 0 [-o | --overhang]\n");
        // std::cerr << "[" << GetCurTime() << "] overhang should be >= 0 [-o | --overhang]\n";
        return 0;
    }
    if(mapq < 0 || mapq > 60) {
        LOG(WARNING)("mapq should be >= 0 && <= 60 [-q | --mapq]\n");
        // std::cerr << "[" << GetCurTime() << "] mapq should be >= 0 && <= 60 [-q | --mapq]\n";
        return 0;
    }
    if(align_rate < 0 || align_rate > 1) {
        LOG(WARNING)("align_rate should be >= 0.0 && <= 1.0 [-f | --align_rate]\n");
        // std::cerr << "[" << GetCurTime() << "] align_rate should be >= 0.0 && <= 1.0 [-f | --align_rate]\n";
        return 0;
    }
    if(overhang_rate < 0 || overhang_rate > 1) {
        LOG(WARNING)("overhang_rate should be >= 0.0 && <= 1.0 [-p | --overhang_rate]\n");
        // std::cerr << "[" << GetCurTime() << "] overhang_rate should be >= 0.0 && <= 1.0 [-p | --overhang_rate]\n";
        return 0;
    }
    if(identity < 0 || identity > 1) {
        LOG(WARNING)("identity should be >= 0.0 && <= 1.0 [-i | --identity]\n");
        // std::cerr << "[" << GetCurTime() << "] identity should be >= 0.0 && <= 1.0 [-i | --identity]\n";
        return 0;
    }
    return 1;
}

bool Variants::Valid(bam1_t *record, Options &opt) {
    if(record->core.qual < opt.mapq) return false;
    if((record->core.flag & 0x2316) > 0) return false;
    return true;
    int length = 0, aligned_length = 0, overhang = 0;
    uint32_t *cigar = bam_get_cigar(record);
    for(int i = 0; i < record->core.n_cigar; ++i) {
        int op = cigar[i] & 0xf;
        int len = cigar[i] >> 4;
        if(op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            length += len;
            aligned_length += len;
        } else if(op == BAM_CINS) {
            length += len;
        } else if(op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) {
            length += len;
            overhang = std::max(len, overhang);
        }
    }
    if(length < opt.min_length) return false;
    if(aligned_length < opt.min_align_length || aligned_length < opt.align_rate * length) return false;
    if(overhang > opt.overhang || overhang > opt.overhang_rate * length) return false;
    if(aligned_length * 1.0 / length < opt.identity) return false;
    return true;
}

uint8_t TABLE[256] = {
    255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255, 
    255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255, 
    255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255, 
    255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255, 
    255, 0,   255, 1,    255, 255, 255, 2,    255, 255, 255, 255,  255, 255, 255, 255, 
    255, 255, 255, 255,  3,   255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255, 
    255, 0,   255, 1,    255, 255, 255, 2,    255, 255, 255, 255,  255, 255, 255, 255, 
    255, 255, 255, 255,  3,   255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255, 
    255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255, 
    255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255, 
    255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255, 
    255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255, 
    255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255, 
    255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255, 
    255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255, 
    255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255,  255, 255, 255, 255, 
};

Variants::Variants(const ReadStore &rs, const std::vector<bam1_t*> &record, std::string &ctg, std::vector<long> &p, Options &opt): 
        ctg_store_(rs), bam_record_(record), ctg_name_(ctg), pos_(p) {
    snps_.assign(ctg_store_.GetSeqLength(ctg_name_), Table());
    if(ctg == "000003F_001") {
        for(auto i: pos_) {
            std::cerr << i << " ";
        }
        std::cerr << std::endl;
    }
    if(!opt.vfname.empty()) {
        GzFileReader reader(opt.vfname);
        if(reader.Valid()) {
            bool found = false;
            std::string line = reader.GetNoEmptyLine();
            while(!line.empty()) {
                if(line[0] != '#' && line.substr(0, ctg.size()) == ctg) {
                    found = true;
                    auto items = SplitString(line, '\t');
                    if(items[3].size() == 1 && items[4].size() == 1) {
                        long pos = atol(items[1].c_str()) - 1;
                        // auto iter = std::lower_bound(pos_.begin(), pos_.end(), pos);
                        // if(iter == pos_.end() || *iter != pos) continue;
                        // if(std::find(pos_.begin(), pos_.end(), pos) == pos_.end()) continue;
                        snps_[pos].SetRef(TABLE[items[3][0]]);
                        snps_[pos].IncreaseM(TABLE[items[3][0]]);
                        snps_[pos].IncreaseM(TABLE[items[4][0]]);
                        snps_[pos].var[0] = TABLE[items[3][0]];
                        snps_[pos].var[1] = TABLE[items[4][0]];
                    }
                } else if(line[0] != '#' && line.substr(0, ctg.size()) != ctg) {
                    if(found) break;
                }
                line = reader.GetNoEmptyLine();
            }
        }
    } else {
        FindSNPInContig(opt);
        ConfirmSNPs();
    }
    FindSNPInReads(opt);
}

void Variants::FindSNPInContig(Options &opt) {
    if(bam_record_.empty()) {
        LOG(WARNING)("There is no record of %s in map file\n", ctg_name_.c_str());
        // std::cerr << "[" << GetCurTime() << "] There is no record of " << ctg_name_ << " in map file\n";
        return;
    }
    const DNASeq &ctg_seq = ctg_store_.GetSeq(ctg_name_);
    for(std::size_t i = 0; i < ctg_seq.Size(); ++i) {
        snps_[i].SetRef(ctg_seq[i]);
    }
    for(auto &record: bam_record_) {
        if(!Valid(record, opt)) continue;
        uint8_t *rd_seq = bam_get_seq(record); // 1:A  2:C  4:G  8:T
        std::size_t ctg_pos = record->core.pos;
        std::size_t rd_pos = 0;
        uint32_t *cigar = bam_get_cigar(record);
        for(uint32_t i = 0; i < record->core.n_cigar; ++i) {
            int len = bam_cigar_oplen(cigar[i]);
            int op = bam_cigar_op(cigar[i]);
            // int len = cigar[i] >> 4;
            // int op = cigar[i] & 0xf;
            if(op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                if(opt.check && len >= opt.snp_match_length) {
                    for(int j = opt.snp_match_length / 2; j < len - opt.snp_match_length / 2; ++j) {
                        std::size_t ctg_i = ctg_pos + j;
                        uint8_t rd_base = Convert(bam_seqi(rd_seq, rd_pos + j));
                        bool s = true;
                        for(int jj = j - opt.snp_match_length / 2; jj < j + opt.snp_match_length / 2 + 1; ++jj) {
                            if(jj != j) {
                                std::size_t ctg_i = ctg_pos + jj;
                                uint8_t rd_base = Convert(bam_seqi(rd_seq, rd_pos + jj));
                                if(ctg_seq[ctg_i] != rd_base) {
                                    s = false;
                                    break;
                                }
                            }
                        }
                        if(s) snps_[ctg_i].IncreaseM(rd_base);
                    }
                } else {
                    for(int j = 0; j < len; ++j) {
                        std::size_t ctg_i = ctg_pos + j;
                        uint8_t rd_base = Convert(bam_seqi(rd_seq, rd_pos + j));
                        snps_[ctg_i].IncreaseM(rd_base);
                    }
                }
                ctg_pos += len;
                rd_pos += len;
            } else if(op == BAM_CINS) {
                rd_pos += len;
            } else if(op == BAM_CDEL) {
                ctg_pos += len;
            } else if(op == BAM_CSOFT_CLIP) {
                rd_pos += len;
            }
        }
        for(std::size_t i = record->core.pos; i < ctg_pos; ++i) snps_[i].IncreaseCov();
    }
}

void Variants::ConfirmSNPs() {
    for(std::size_t i = 0; i < snps_.size(); ++i) {
        snps_[i].Confirm();
    }
}

void Variants::FindSNPInReads(Options &opt) {
    for(auto &record: bam_record_) {
        if(!Valid(record, opt)) {
            bam_destroy1(record);
            continue;
        }
        std::unordered_map<int, std::array<int, 4>> cand_snps;
        uint8_t *rd_seq = bam_get_seq(record);
        std::size_t rd_pos = 0;
        std::size_t ctg_pos = record->core.pos;
        uint32_t *cigar = bam_get_cigar(record);
        for(int i = 0; i < record->core.n_cigar; ++i) {
            int len = bam_cigar_oplen(cigar[i]);
            int op = bam_cigar_op(cigar[i]);
            // int len = cigar[i] >> 4;
            // int op = cigar[i] & 0xf;
            if(op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                if(opt.check && len >= opt.snp_match_length) {
                    for(int j = opt.snp_match_length / 2; j < len - opt.snp_match_length / 2; ++j) {
                        std::size_t ctg_i = ctg_pos + j;
                        uint8_t rd_base = Convert(bam_seqi(rd_seq, rd_pos + j));
                        if(snps_[ctg_i].Valid()) {
                            if(snps_[ctg_i].AtM(rd_base)) {
                                read_info_[bam_get_qname(record)].push_back( {(int)ctg_i, (int)rd_pos + j, rd_base, snps_[ctg_i].Offset(rd_base)} );
                            }
                        }
                    }
                } else {
                    for(int j = 0; j < len; ++j) {
                        std::size_t ctg_i = ctg_pos + j;
                        uint8_t rd_base = Convert(bam_seqi(rd_seq, rd_pos + j));
                        if(snps_[ctg_i].Valid()) {
                            if(snps_[ctg_i].AtM(rd_base)) {
                                read_info_[bam_get_qname(record)].push_back( {(int)ctg_i, (int)rd_pos + j, rd_base, snps_[ctg_i].Offset(rd_base)} );
                            }
                        }
                    }
                }
                rd_pos += len;
                ctg_pos += len;
            } else if(op == BAM_CINS) {
                rd_pos += len;
            } else if(op == BAM_CDEL) {
                ctg_pos += len;
            } else if(op == BAM_CSOFT_CLIP) {
                rd_pos += len;
            }
        }
        bam_destroy1(record);
    }
}

void Variants::DumpSNP(std::ostream &out, bool dump_all) const {
    for(std::size_t i = 0; i < snps_.size(); ++i) {
        if(!dump_all && !snps_[i].Valid()) continue;
        out << ctg_name_ << "\t" << i;
        for(std::size_t j = 0; j < 6; ++j) {
            out << "\t" << snps_[i].counts[j];
        }
        out << "\t" << snps_[i].Valid() << "\n";
    }
}

void Variants::DumpReadInfo(std::ostream &out) const {
    for(auto &ri: read_info_) {
        if(ri.second.size() > 0) {
            out << ctg_name_ << "\t" << ri.first << "\t" << ri.second[0][0] - ri.second[0][1];
            for(auto &v: ri.second) {
                int choice = v[3];
                if(choice == 0 || choice == 1) {
                    out << "\t" << v[0] << "|" << int(snps_[v[0]].var[choice]) << "|" << v[1] << "|" << v[2];
                } else {
                    out << "\t" << v[0] << "|" << -1 << "|" << v[1] << "|" << v[2];
                }
            }
            out << "\n";
        }
    }
}

void SNP::Load() {
    ctg_store_.Load(opt_.ctgfname);
    reader_.Initialize(opt_.bamfname);
}

std::unordered_map<std::string, std::vector<long>> SNP::GetVariantsWithClair3() const {
    std::unordered_map<std::string, std::vector<long>> vars;
    std::unordered_map<std::string, std::unordered_map<long, std::string>> variants;
    // load variants in primary
    GzFileReader preader(opt_.pcfname);
    if(preader.Valid()) {
        std::string line = preader.GetNoEmptyLine();
        while(!line.empty()) {
            if(line[0] != '#') {
                auto items = SplitString(line, '\t');
                if(items[3].size() == 1 && items[4].size() == 1) {
                    long pos = atol(items[1].c_str()) - 1;
                    variants[items[0]][pos] = items[3] + items[4];
                }
            }
            line = preader.GetNoEmptyLine();
        }
    } else {
        LOG(ERROR)("Failed to open %s for reading", opt_.pcfname.c_str());
        // std::cerr << "[" << GetCurTime() << "] Failed to open " << opt_.pcfname << " for reading\n";
        // exit(EXIT_FAILURE);
    }
    // load variants in alt
    GzFileReader areader(opt_.acfname);
    if(areader.Valid()) {
        std::string line = areader.GetNoEmptyLine();
        while(!line.empty()) {
            if(line[0] != '#') {
                auto items = SplitString(line, '\t');
                if(items[3].size() == 1 && items[4].size() == 1) {
                    long pos = atol(items[1].c_str()) - 1;
                    variants[items[0]][pos] = items[4] + items[3];
                }
            }
            line = areader.GetNoEmptyLine();
        }
    } else {
        LOG(ERROR)("Failed to open %s for reading", opt_.acfname.c_str());
        // std::cerr << "[" << GetCurTime() << "] Failed to open " << opt_.acfname << " for reading\n";
        // exit(EXIT_FAILURE);
    }
    // std::cerr << "variants[000003F_001] size " << variants["000003F_001"].size() << std::endl;
    // load paf
    std::unordered_map<std::string, std::vector<std::array<std::string, 13>>> paf;
    std::ifstream in(opt_.mapfname, std::ios::in);
    if(!in.is_open()) {
        LOG(ERROR)("Could not open %s for reading", opt_.mapfname.c_str());
        // std::cerr << "[" << GetCurTime() << "] Error could not open " << opt_.mapfname << " for reading\n";
        // exit(EXIT_FAILURE);
    }
    std::string line;
    while(std::getline(in, line)) {
        auto items = SplitString(line, '\t');
        std::array<std::string, 13> tmp;
        for(std::size_t i = 0; i < items.size(); ++i) {
            if(i < 12) tmp[i] = items[i];
            else if(items[i].substr(0, 5) == "cg:Z:") {
                tmp[12] = items[i];
                break;
            }
        }
        paf[items[5]].emplace_back(tmp);
    }
    in.close();
    // sort by position on target
    std::vector<std::string> keys;
    for(auto &ctg: paf) {
        keys.emplace_back(ctg.first);
        std::sort(ctg.second.begin(), ctg.second.end(), [](std::array<std::string, 13> &a, std::array<std::string, 13> &b) {
            return atol(a[7].c_str()) < atol(b[7].c_str());
        });
    }
    std::sort(keys.begin(), keys.end());

    auto find_diff = [&](std::vector<std::array<std::string, 13>> &records, std::string &ctg) -> std::unordered_map<std::string, std::vector<uint64_t>> {
        std::unordered_map<std::string, std::vector<uint64_t>> diff;
        std::regex pattern("(\\d+)([MIDX=]?)");
        for(auto &record: records) {
            long ref_pos = atol(record[7].c_str());
            std::string query = record[0];
            long query_pos = atol(record[2].c_str());
            if(record[4] == "-") query_pos = atol(record[3].c_str()) - 1;
            std::string cigar = record[12].substr(5);
            for(std::sregex_iterator iter(cigar.begin(), cigar.end(), pattern), end; iter != end; ++iter) {
                int len =  atoi(iter->str(1).c_str());
                std::string op = iter->str(2);
                if(op == "M" || op == "=") {
                    ref_pos += len;
                    if(record[4] == "+") query_pos += len;
                    else query_pos -= len;
                } else if(op == "X") {
                    for(int i = 0; i < len; ++i) {
                        // if(query == "000003F_001") {
                        //     std::cerr << "variants[" << ctg << "][" << ref_pos + i << "] " << variants[ctg][ref_pos + i] << " ";
                        //     std::cerr << "variamts[" << query << "][" << query_pos + i << "] " << variants[query][query_pos + i] << std::endl;
                        // }
                        if(record[4] == "+") {
                            if(variants.find(ctg) != variants.end() && 
                                variants.find(query) != variants.end() && 
                                variants[ctg].find(ref_pos + i) != variants[ctg].end() && 
                                variants[query].find(query_pos + i) != variants[query].end() && 
                                variants[ctg][ref_pos + i] == variants[query][query_pos + i]) {
                                    diff[ctg].emplace_back(ref_pos + i);
                                    diff[query].emplace_back(query_pos + i);
                            }
                        } else {
                            if(variants.find(ctg) != variants.end() && 
                                variants.find(query) != variants.end() && 
                                variants[ctg].find(ref_pos + i) != variants[ctg].end() && 
                                variants[query].find(query_pos - i) != variants[query].end() && 
                                variants[ctg][ref_pos + i] == variants[query][query_pos - i]) {
                                    diff[ctg].emplace_back(ref_pos + i);
                                    diff[query].emplace_back(query_pos - i);
                            }
                        }
                    }
                    ref_pos += len;
                    if(record[4] == "+") query_pos += len;
                    else query_pos -= len;
                } else if(op == "D") {
                } else if(op == "I") {
                    if(record[4] == "+") query_pos += len;
                    else query_pos -= len;
                }
            }
            if(record[4] == "-") {
                std::sort(diff[query].begin(), diff[query].end());
            }
        }
        return diff;
    };

    std::mutex mtx;
    std::atomic<std::size_t> index { 0 };
    auto func = [&](std::size_t tid) {
        for(auto cur = index.fetch_add(1); cur < keys.size(); cur = index.fetch_add(1)) {
            std::string ctg = keys[cur];
            auto diff = find_diff(paf[ctg], ctg);
            {
                std::lock_guard<std::mutex> lock(mtx);
                // if(ctg == "000003F") std::cerr << "diff[000003F_001] size " << diff["000003F_001"].size() << std::endl;
                for(auto dif: diff) {
                    for(auto pos: dif.second) {
                        vars[dif.first].emplace_back(pos);
                    }
                }
            }
        }
    };
    MultiThreads(opt_.threads, func);
    // std::cerr << "vars[000003F_001] size " << vars["000003F_001"].size() << std::endl;
    return vars;
}

void SNP::Run() {
    Load();
    auto variants = GetVariantsWithClair3();
    std::mutex mtx, mtx_err, mtx_read;
    std::ofstream out_var(opt_.prefix + ".variants");
    std::ofstream out_ri(opt_.prefix + ".readinfo");

    auto ctg_names = reader_.LoadRefName();

    auto read = [&](const std::string &ctg) -> std::vector<bam1_t*> {
        std::lock_guard<std::mutex> lock(mtx_read);
        return reader_.Load(ctg);
    };

    auto dump = [&](std::ostringstream &oss_snp, std::ostringstream &oss_ri) {
        std::lock_guard<std::mutex> lock(mtx);
        out_var << oss_snp.str();
        oss_snp.str("");
        out_ri << oss_ri.str();
        oss_ri.str("");
    };

    std::atomic<std::size_t> index{ 0 };
    auto task = [&](std::size_t tid) {
        std::ostringstream oss_snp;
        std::ostringstream oss_ri;
        for(auto cur = index.fetch_add(1); cur < ctg_names.size(); cur = index.fetch_add(1)) {
            std::string ctg_name = ctg_names[cur];
            std::vector<bam1_t*> record = read(ctg_name);
            if(record.empty()) {
                std::lock_guard<std::mutex> lock(mtx_err);
                std::cout << ctg_name << "\n";
                continue;
            }
            Variants var(ctg_store_, record, ctg_name, variants[ctg_name], opt_);
            var.DumpSNP(oss_snp, opt_.write_all);
            var.DumpReadInfo(oss_ri);
            dump(oss_snp, oss_ri);
        }
    };
    MultiThreads(opt_.threads, task);
}

int main(int argc, char **argv) {
    Options opt;
    if(ParseArgument(argc, argv, opt) == 0) exit(EXIT_FAILURE);
    SNP snp(opt);
    snp.Run();
    return 0;
}