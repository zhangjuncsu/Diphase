#include "filter.hpp"
#include "utility.hpp"

#include "getopt.h"

#include <mutex>
#include <regex>
#include <atomic>
#include <cassert>
#include <fstream>
#include <iostream>

int Options::Check() {
    // if(varfname.empty()) {
    //     std::cerr << "[" << GetCurTime() << "] Please specify a variant file name [-v | --var]\n";
    //     return 1;
    // }
    if(threads <= 0) {
        std::cerr << "[" << GetCurTime() << "] threads must be >= 1 [-t | --threads]\n";
        return 1;
    }
    return 0;
    // if(mapq < 0) {
    //     std::cerr << "[" << GetCurTime() << "] mapping quality must be >= 0 [-q | --mapq]\n";
    //     return 1;
    // }
}

void USAGE(Options &opt) {
    std::cerr << "USAGE:";
    std::cerr << "\tfilter [options] in.bam\n\n";
    std::cerr << "\tOptions:\n";
    std::cerr << "\t\t-o | --output  <FILE> output file name [standard output]\n";
    std::cerr << "\t\t-v | --var     <FILE> variants file name\n";
    std::cerr << "\t\t-t | --threads <INT> number of threads [" << opt.threads << "]\n";
    std::cerr << "\t\t-q | --mapq    <INT> have mapping quality >= INT [" << opt.mapq << "]\n";
    std::cerr << "\t\t-c | --clip    <INT> overhang to be filtered [" << opt.clip << "]\n";
    std::cerr << "\t\t-h | --help          display this message\n";
}

int ParseArgument(int argc, char **argv, Options &opt) {
    struct option long_opt[] = {
        {"output",  required_argument,  NULL,   'o'}, 
        {"var",     optional_argument,  NULL,   'v'}, 
        {"threads", optional_argument,  NULL,   't'}, 
        {"mapq",    optional_argument,  NULL,   'q'}, 
        {"clip",    optional_argument,  NULL,   'c'}, 
        {"help",    no_argument,        NULL,   'h'}, 
        {0, 0,  0,  0},
    };
    const char *short_opt = "v:t:o:q:c:h";
    int c;
    while((c = getopt_long(argc, argv, short_opt, long_opt, NULL)) != -1) {
        if(c == 'v') opt.varfname = optarg;
        else if(c == 'o') opt.outfname = optarg;
        else if(c == 't') opt.threads = atoi(optarg);
        else if(c == 'q') opt.mapq = atoi(optarg);
        else if(c == 'c') opt.clip = atoi(optarg);
        else if(c == 'h') {
            USAGE(opt);
            exit(EXIT_SUCCESS);
        } else if(c == ':') {
            std::cerr << "[" << GetCurTime() << "] ERROR missing option argument in " << optopt << "\n";
            return 1;
        } else if(c == '?') {
            std::cerr << "[" << GetCurTime() << "] ERROR unknown option in " << optopt << "\n";
            return 1;
        }
    }
    if(optind != argc) opt.infname = argv[optind];
    return opt.Check();
}

void Filter::LoadVariants() {
    std::string line;
    std::ifstream in(opt_.varfname);
    if(!in.is_open()) {
        std::cerr << "[" << GetCurTime() << "] Could not open " << opt_.varfname << " for reading\n";
        exit(EXIT_FAILURE);
    }
    while(std::getline(in, line)) {
        auto items = SplitString(line, '\t');
        vars_[items[0]].emplace_back(atol(items[1].c_str()));
    }
}

bool Filter::FilterVariants(bam1_t *record) {
    // bool ret = true;
    std::regex pattern("(\\w+):(\\d+)-(\\d+)");
    std::string ref_name = reader_.Header()->target_name[record->core.tid];
    int beg = 0;
    std::smatch result;
    if(std::regex_search(ref_name, result, pattern)) {
        ref_name = result[1].str();
        beg = atoi(result[2].str().c_str());
    }
    if(vars_.find(ref_name) == vars_.end()) return false;
    int64_t ref_pos = record->core.pos;
    int ref_len = bam_cigar2rlen(record->core.n_cigar, bam_get_cigar(record));
    auto low_bound = std::lower_bound(vars_[ref_name].begin(), vars_[ref_name].end(), ref_pos + beg);
    if(!(low_bound != vars_[ref_name].end() && *low_bound <= ref_pos + beg + ref_len)) return true;
    // if(ret) return ret;

    uint8_t *md_str = bam_aux_get(record, "MD");
    if(md_str != NULL) {
        std::string md = bam_aux2Z(md_str);
        auto s = md.begin();
        while(s < md.end()) {
            auto e = std::find_if(s, md.end(), [](char c) { return !::isdigit(c); });
            int num = atoi(std::string(s, e).c_str());
            ref_pos += num;
            s = std::find_if(e, md.end(), [](char c) { return ::isdigit(c); });
            if(s == md.end()) break;
            int distance = std::distance(e, s);
            if(e[0] == '^') ref_pos = ref_pos + distance - 1;
            else {
                for(int i = 0; i < distance; ++i) {
                    if(std::binary_search(vars_[ref_name].begin(), vars_[ref_name].end(), ref_pos + i + beg)) {
                        // ret = true;
                        // break;
                        return true;
                    }
                }
                ref_pos += distance;
                // if(ret) break;
            }
        }
    }
    // return ret;
    return false;
}

int Filter::FilterVariants(bam1_t *a, bam1_t *b) {
    uint8_t *md_a = bam_aux_get(a, "MD");
    uint8_t *md_b = bam_aux_get(b, "MD");
    if(md_a != NULL && md_b != NULL) {
        std::string md_str_a = bam_aux2Z(md_a);
        std::string md_str_b = bam_aux2Z(md_b);
        if(md_str_a == md_str_b) return 0;
        int64_t ref_a = a->core.pos;
        int64_t ref_b = a->core.pos;
        int ref_len_a = bam_cigar2rlen(a->core.n_cigar, bam_get_cigar(a));
        int ref_len_b = bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b));
        std::string name_a = bam_get_qname(a);
        std::string name_b = bam_get_qname(b);
        std::regex pattern("(\\w+):(\\d+)-(\\d+)");
        std::smatch result;
        int beg_a = 0, beg_b = 0;
        if(std::regex_search(name_a, result, pattern)) {
            name_a = result[1].str();
            beg_a = atoi(result[2].str().c_str());
        }
        if(std::regex_search(name_b, result, pattern)) {
            name_b = result[1].str();
            beg_b = atoi(result[2].str().c_str());
        }
        bool span_a = false, span_b = false;
        auto lb_a = std::lower_bound(vars_[name_a].begin(), vars_[name_a].end(), ref_a + beg_a);
        if(lb_a != vars_[name_a].end() && *lb_a <= ref_a + beg_a + ref_len_a) {
            span_a = true;
        }
        auto lb_b = std::lower_bound(vars_[name_b].begin(), vars_[name_b].end(), ref_b + beg_b);
        if(lb_b != vars_[name_b].end() && *lb_b <= ref_b + beg_b + ref_len_b) {
            span_b = true;
        }
        if(!span_a || !span_b) return 0;
        bool fa = false, fb = false;
        if(span_a) {
            auto s = md_str_a.begin();
            while(s < md_str_a.end()) {
                auto e = std::find_if(s, md_str_a.end(), [](char c) { return !::isdigit(c); });
                int num = atoi(std::string(s, e).c_str());
                ref_a += num;
                s = std::find_if(e, md_str_a.end(), [](char c) { return ::isdigit(c); });
                if(s == md_str_a.end()) break;
                int distance = std::distance(e, s);
                if(e[0] == '^') ref_a += distance - 1;
                else {
                    if(std::binary_search(vars_[name_a].begin(), vars_[name_a].end(), ref_a + beg_a)) {
                        fa = true;
                    }
                    ref_a += 1;
                }
                if(fa) break;
            }
        }
        if(span_b) {
            auto s = md_str_b.begin();
            while(s < md_str_b.end()) {
                auto e = std::find_if(s, md_str_b.end(), [](char c) { return !::isdigit(c); });
                int num = atoi(std::string(s, e).c_str());
                ref_b += num;
                s = std::find_if(e, md_str_b.end(), [](char c) { return ::isdigit(c); });
                if(s == md_str_b.end()) break;
                int distance = std::distance(e, s);
                if(e[0] == '^') ref_b += distance - 1;
                else {
                    if(std::binary_search(vars_[name_b].begin(), vars_[name_b].end(), ref_b + beg_b)) {
                        fb = true;
                    }
                    ref_b += 1;
                }
                if(fb) break;
            }
        }
        if(!fa && fb) return 1;
        else if(fa && !fb) return 2;
        else return 0;
    }
    return 0;
}

void Filter::Run1() {
    std::ios::sync_with_stdio(false);
    reader_.Initialize(opt_.infname);
    LoadVariants();
    samFile *out = sam_open(opt_.outfname.c_str(), "wb");
    if(sam_hdr_write(out, reader_.Header()) != 0) {
        std::cerr << "[" << GetCurTime() << "] Failed to write header to " << opt_.outfname << "\n";
        exit(EXIT_FAILURE);
    }

    std::size_t batch_size = 1000;
    std::vector<std::vector<bam1_t*>> reads(1000);
    while(1) {
        std::size_t size = reader_.LoadBatchPair(reads);
        if(size == 0) break;
        std::mutex mtx;
        std::atomic<std::size_t> index { 0 };
        std::vector<std::vector<int8_t>> rm(size);
        auto func = [&](std::size_t tid) {
            for(auto cur = index.fetch_add(1); cur < size; cur = index.fetch_add(1)) {
                auto records = reads[cur];
                rm[cur].resize(records.size(), 0);
                if(records.size() < 2 || (records.front()->core.flag & 0x80) || (records.back()->core.flag & 0x40)) {
                    for(std::size_t i = 0; i < records.size(); ++i) rm[cur][i] = 1;
                    continue;
                }
                bool ff = false, sf = false;
                for(std::size_t i = 0; i < records.size(); ++i) {
                    if((ff && (records[i]->core.flag & 0x40)) || (sf && (records[i]->core.flag & 0x80))) {
                        rm[cur][i] = 1;
                        continue;
                    }
                    if(records[i]->core.qual < opt_.mapq || FilterVariants(records[i])) {
                        rm[cur][i] = 1;
                    } else {
                        if(records[i]->core.flag & 0x40) ff = true;
                        else sf = true;
                    }
                }
                int num = std::count(rm[cur].begin(), rm[cur].end(), 0);
                if(num != 2) {
                    for(std::size_t i = 0; i < rm[cur].size(); ++i) rm[cur][i] = 1;
                    if(num > 2) {
                        std::lock_guard<std::mutex> lock(mtx);
                        std::cerr << bam_get_qname(records[0]) << "\n";
                    }
                }
            }
        };
        MultiThreads(opt_.threads, func);
        for(std::size_t i = 0; i < size; ++i) {
            for(std::size_t j = 0; j < reads[i].size(); ++j) {
                if(rm[i][j] != 1) sam_write1(out, reader_.Header(), reads[i][j]);
                bam_destroy1(reads[i][j]);
            }
            reads[i].clear();
        }
    }
    sam_close(out);
}

struct bamCmp {
    bool operator()(const bam1_t *a, const bam1_t *b) {
        int am = 0, bm = 0;
        uint32_t *cigar_a = bam_get_cigar(a);
        uint32_t *cigar_b = bam_get_cigar(b);
        for(int i = 0; i < a->core.n_cigar; ++i) {
            int size = cigar_a[i] >> 4;
            int op = cigar_a[i] & 0xf;
            if(op == BAM_CMATCH) am += size;
        }
        for(int i = 0; i < b->core.n_cigar; ++i) {
            int size = cigar_b[i] >> 4;
            int op = cigar_b[i] & 0xf;
            if(op == BAM_CMATCH) bm += size;
        }
        if(am != bm) return am > bm;
        else return strcmp(bam_get_qname(a), bam_get_qname(b));
    }
};

void Filter::Run2() {
    std::ios::sync_with_stdio(false);
    // std::cin.tie(NULL);
    reader_.Initialize(opt_.infname);
    LoadVariants();
    htsFile *out = sam_open(opt_.outfname.c_str(), "wb");
    if(sam_hdr_write(out, reader_.Header()) != 0) {
        std::cerr << "[" << GetCurTime() << "] Failed to write header to " << opt_.outfname << "\n";
        exit(EXIT_FAILURE);
    }

    const std::size_t batch_size = 1000;
    std::vector<std::vector<bam1_t*>> reads(batch_size);
    int count = 0;
    int remaining = 0;
    while(1) {
        std::size_t size = reader_.LoadBatchPair(reads);
        if(size == 0) break;
        std::atomic<std::size_t> index { 0 };
        std::mutex mtx;
        std::vector<std::vector<bam1_t*>> result(size);

        auto cmp = [&](bam1_t *a, bam1_t *b) {
            int am = 0, bm = 0;
            uint32_t *cigar_a = bam_get_cigar(a);
            uint32_t *cigar_b = bam_get_cigar(b);
            for(int i = 0; i < a->core.n_cigar; ++i) {
                int size = cigar_a[i] >> 4;
                int op = cigar_a[i] & 0xf;
                if(op == BAM_CMATCH) am += size;
            }
            for(int i = 0; i < b->core.n_cigar; ++i) {
                int size = cigar_b[i] >> 4;
                int op = cigar_b[i] & 0xf;
                if(op == BAM_CMATCH) bm += size;
            }
            if(am != bm) return am > bm;
            else {
                std::regex pattern("(\\w+):(\\d+)-(\\d+)");
                std::smatch result;
                std::string name_a = reader_.Header()->target_name[a->core.tid];
                std::string name_b = reader_.Header()->target_name[b->core.tid];
                if(std::regex_search(name_a, result, pattern)) name_a = result[1].str();
                if(std::regex_search(name_b, result, pattern)) name_b = result[1].str();
                return name_a < name_b;
            }
        };

        auto filt = [&](std::vector<bam1_t*> &records, std::vector<bam1_t*> &res) {
            std::size_t mid = 0;
            for(std::size_t i = 0; i < records.size(); ++i) {
                if(records[i]->core.flag & 0x80) break;
                mid++;
            }
            if(mid == 0) return;
            std::regex pattern("(\\w+):(\\d+)-(\\d+)");
            // get the first segment mapping record
            bool fe = true;
            if(mid == 1) {
                res[0] = bam_dup1(records[0]);
                fe = true;
            } else {
                std::sort(records.begin(), records.begin() + mid, cmp);
                for(std::size_t i = 0; i < mid - 1; ++i) {
                    std::string name_a = reader_.Header()->target_name[records[i]->core.tid];
                    std::string name_b = reader_.Header()->target_name[records[i + 1]->core.tid];
                    std::smatch result;
                    if(std::regex_search(name_a, result, pattern)) name_a = result[1].str();
                    if(std::regex_search(name_b, result, pattern)) name_b = result[1].str();
                    std::size_t cmp_len = std::min(name_a.size(), name_b.size());
                    if(name_a.substr(0, cmp_len) == name_b.substr(0, cmp_len)) {
                        bool ff = FilterVariants(records[i]);
                        bool fs = FilterVariants(records[i + 1]);
                        if(ff && !fs) {
                            res[0] = bam_dup1(records[i + 1]);
                            fe = false;
                            break;
                        } else if(!ff && fs) {
                            res[0] = bam_dup1(records[i]);
                            fe = false;
                            break;
                        } else {
                            i += 1;
                            continue;
                        }
                    } else {
                        res[0] = bam_dup1(records[i]);
                        fe = false;
                        break;
                    }
                }
            }
            if(fe) {
                res.clear();
                return;
            }
            // get the second segment mapping record
            bool se = true;
            if(mid == records.size() - 1) {
                res[0] = bam_dup1(records[mid]);
                se = false;
            } else {
                std::sort(records.begin() + mid, records.end(), cmp);
                for(std::size_t i = mid; i < records.size() - 1; ++i) {
                    std::string name_a = reader_.Header()->target_name[records[i]->core.tid];
                    std::string name_b = reader_.Header()->target_name[records[i + 1]->core.tid];
                    std::smatch result;
                    if(std::regex_search(name_a, result, pattern)) name_a = result[1].str();
                    if(std::regex_search(name_b, result, pattern)) name_b = result[1].str();
                    std::size_t cmp_len = std::min(name_a.size(), name_b.size());
                    if(name_a.substr(0, cmp_len) == name_b.substr(0, cmp_len)) {
                        bool ff = FilterVariants(records[i]);
                        bool fs = FilterVariants(records[i + 1]);
                        if(ff && !fs) {
                            res[1] = bam_dup1(records[i + 1]);
                            se = false;
                            break;
                        } else if(!ff && fs) {
                            res[1] = bam_dup1(records[i]);
                            se = false;
                            break;
                        } else {
                            i += 1;
                            continue;
                        }
                    } else {
                        res[1] = bam_dup1(records[i]);
                        se = false;
                        break;
                    }
                }
            }
            if(se) {
                bam_destroy1(res[0]);
                res.clear();
            }
        };

        auto filt1 = [&](std::vector<bam1_t*> &records, std::vector<bam1_t*> &res) {
            std::size_t mid = 0;
            for(std::size_t i = 0; i < records.size(); ++i) {
                if(records[i]->core.flag & 0x80) break;
                mid++;
            }
            if(mid == 0) return;
            for(std::size_t i = 0; i < mid; ++i) {
                if(!FilterVariants(records[i])) {
                    res.emplace_back(bam_dup1(records[i]));
                    break;
                }
            }
            if(res.empty()) return;
            for(std::size_t i = mid; i < records.size(); ++i) {
                if(!FilterVariants(records[i])) {
                    res.emplace_back(bam_dup1(records[i]));
                    break;
                }
            }
            if(res.size() != 2) {
                for(auto &b: res) {
                    bam_destroy1(b);
                }
                res.clear();
                return;
            }
            assert(res.size() == 2);
        };

        auto func = [&](std::size_t tid) {
            for(auto cur = index.fetch_add(1); cur < size; cur = index.fetch_add(1)) {
                auto& records = reads[cur];
                std::stable_sort(records.begin(), records.end(), [](bam1_t *a, bam1_t *b) { return (a->core.flag & 0x80) < (b->core.flag & 0x40); });
                if((records.front()->core.flag) & 0x80 || (records.back()->core.flag & 0x40)) {
                    continue;
                } 
                // else if(records.size() == 2) {
                //     result[cur].emplace_back(bam_dup1(records.front()));
                //     result[cur].emplace_back(bam_dup1(records.back()));
                //     continue;
                // }
                result[cur].resize(2);
                filt1(records, result[cur]);
            }
        };
        MultiThreads(opt_.threads, func);
        for(auto &records: reads) {
            for(auto &b: records) bam_destroy1(b);
            records.clear();
        }
        for(auto &record: result) {
            if(record.empty()) continue;
            sam_write1(out, reader_.Header(), record[0]);
            sam_write1(out, reader_.Header(), record[1]);
            bam_destroy1(record[0]);
            bam_destroy1(record[1]);
            remaining += 1;
        }
        count += size;
        if((count % 1000000) == 0) {
            std::cerr << "[" << GetCurTime() << "] Parsed " << count << " read pairs\n";
        }
    }
    hts_close(out);
    std::cerr << "[" << GetCurTime() << "] " << remaining << " " << (remaining / (double)count) * 100 << "% mate pair passed filtering\n";
}
/*
void Filter::Run() {
    std::ios::sync_with_stdio(false);
    // std::cin.tie(NULL);
    reader_.Initialize(opt_.infname);
    // LoadVariants();
    htsFile *out = sam_open(opt_.outfname.c_str(), "wb");
    if(sam_hdr_write(out, reader_.Header()) != 0) {
        std::cerr << "[" << GetCurTime() << "] Failed to write header to " << opt_.outfname << "\n";
        exit(EXIT_FAILURE);
    }
    std::size_t SIZE = 100000;
    std::vector<bam1_t*> reads(SIZE);
    for(std::size_t i = 0; i < SIZE; ++i) {
        reads[i] = bam_init1();
    }
    while(1) {
        std::size_t size = 0;
        for(std::size_t i = 0; i < SIZE; ++i, ++size) {
            int ret = reader_.LoadOneRecord(reads[i]);
            if(ret < 0) break;
        }
        if(size == 0) break;

        std::vector<short> flag(size, 0);
        std::atomic<std::size_t> index { 0 };
        auto func = [&] (std::size_t tid) {
            for(auto cur = index.fetch_add(1); cur < size; cur = index.fetch_add(1)) {
                bam1_t *record = reads[cur];
                uint32_t *cigar = bam_get_cigar(record);
                int clip = 0;
                for(int i = 0; i < record->core.n_cigar; ++i) {
                    int op = cigar[i] & 0xf;
                    int len = cigar[i] >> 4;
                    if(op == BAM_CHARD_CLIP || op == BAM_CSOFT_CLIP) clip += len;
                    // else if(op == BAM_CDEL || op == BAM_CINS) {
                    //     if(len >= 100) {
                    //         flag[cur] = 1;
                    //         break;
                    //     }
                    // }
                }
                if(clip > opt_.clip) flag[cur] = 1;
            }
        };
        MultiThreads(opt_.threads, func);
        for(std::size_t i = 0; i < size; ++i) {
            if(flag[i] == 0)
                sam_write1(out, reader_.Header(), reads[i]);
        }
    }
    for(std::size_t i = 0; i < SIZE; ++i) bam_destroy1(reads[i]);
    hts_close(out);
}
*/
void Filter::Run() {
    std::ios::sync_with_stdio(false);
    reader_.Initialize(opt_.infname);
    std::array<bam1_t*, 2> pair;
    pair[0] = bam_init1();
    pair[1] = bam_init1();
    htsFile *out = sam_open(opt_.outfname.c_str(), "wb");
    if(sam_hdr_write(out, reader_.Header()) != 0) {
        std::cerr << "[" << GetCurTime() << "] Failed to write header to " << opt_.outfname << "\n";
        exit(EXIT_FAILURE);
    }
    std::vector<std::string> names { "000023F:35606921-49533002", "000007F:19439856-23777517 2336700-2655068", 
                                    "000023F_009:0-13865647 4604447-5006217", "000007F_007:0-4334208 2333570-2670989", 
                                    "000013F:29622135-39727347 7567077-7944962", "000013F_012:0-10028763 7464712-7868276" };
    std::unordered_map<std::string, std::array<int, 2>> regions;
    regions["000007F:19439856-23777517"] = std::array<int, 2> { 2336700, 2655068 };
    regions["000007F_007:0-4334208"] = std::array<int, 2> { 2333570, 2670989 };
    regions["000013F:29622135-39727347"] = std::array<int, 2> { 7567077, 7944962 };
    regions["000013F_012:0-10028763"] = std::array<int, 2> { 7464712, 7868276 };
    regions["000023F:35606921-49533002"] = std::array<int, 2> { 4578191, 5035620 };
    regions["000023F_009:0-13865647"] = std::array<int, 2> { 4604447, 5006217 };
    while(1) {
        int r = reader_.LoadPair(pair);
        if(r < 0) break;
        if(pair[0]->core.tid != pair[1]->core.tid) continue;
        std::string name = reader_.Header()->target_name[pair[0]->core.tid];
        if(std::find(names.begin(), names.end(), name) == names.end()) continue;
        sam_write1(out, reader_.Header(), pair[0]);
        sam_write1(out, reader_.Header(), pair[1]);
    }
    std::cerr << "1\n";
    bam_destroy1(pair[0]);
    std::cerr << "2\n";
    bam_destroy1(pair[1]);
    std::cerr << "3\n";
    hts_close(out);
    std::cerr << "4\n";
}

int main(int argc, char **argv) {
    Options opt;
    if(ParseArgument(argc, argv, opt) == 1) exit(EXIT_FAILURE);
    Filter filt(opt);
    filt.Run();
}