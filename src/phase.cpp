#include "phase.hpp"

#include <regex>
#include <mutex>
#include <atomic>
#include <numeric>

#include "getopt.h"

// void Phase::GenerateVariants() {
//     // --eqx
//     BamReader reader;
//     reader.Initialize(opt_.bfname);
    
//     std::ofstream out;
//     if(opt_.dump_var) {
//         out.open(opt_.prefix + ".var.txt");
//         if(!out.is_open()) {
//             std::cerr << "[" << GetCurTime() << "] Warning: could not open " << opt_.prefix << ".var.txt for writing\n";
//         }
//     }
//     std::vector<std::string> ref = reader.LoadRefName();
//     std::atomic<std::size_t> index { 0 };
//     std::mutex mtx;
    
//     auto find_diff = [](std::vector<bam1_t*> &records, std::string &ctg) -> std::unordered_map<std::string, std::vector<uint64_t>> {
//         std::unordered_map<std::string, std::vector<uint64_t>> diff;
//         uint64_t ref_pos = 0;
//         for(auto &record: records) {
//             uint64_t query_pos = 0;
//             ref_pos = record->core.pos;
//             std::string query = bam_get_qname(record);
//             uint32_t *cigar = bam_get_cigar(record);
//             for(uint32_t i = 0; i < record->core.n_cigar; ++i) {
//                 int op = bam_cigar_op(cigar[i]);
//                 int len = bam_cigar_oplen(cigar[i]);
//                 if(op == BAM_CMATCH || op == BAM_CEQUAL) {
//                     ref_pos += len;
//                     query_pos += len;
//                 } else if(op == BAM_CDIFF) {
//                     for(int j = 0; j < len; ++j) {
//                         diff[ctg].emplace_back((ref_pos + j) << 2 | 0);
//                         diff[query].emplace_back((query_pos + j) << 2 | 0);
//                     }
//                     ref_pos += len;
//                     query_pos += len;
//                 } else if(op == BAM_CDEL) {
//                     for(int j = 0; j < len; ++j) {
//                         diff[ctg].emplace_back((ref_pos + j) << 2 | 2);
//                     }
//                     diff[query].emplace_back((query_pos + 1) << 2 | 2);
//                     ref_pos += len;
//                 } else if(op == BAM_CINS) {
//                     for(int j = 0; j < len; ++j) {
//                         diff[query].emplace_back((query_pos + j) << 2 | 1);
//                     }
//                     diff[ctg].emplace_back((ref_pos + 1) << 2 | 1);
//                     query_pos += len;
//                 } else if(op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) {
//                     query_pos += len;
//                 }
//             }
//             bam_destroy1(record);
//         }
//         return diff;
//     };

//     auto func = [&] (std::size_t tid) {
//         for(auto cur = index.fetch_add(1); cur < ref.size(); cur = index.fetch_add(1)) {
//             std::string ctg = ref[cur];
//             auto records = reader.Load(ctg);
//             auto diff = find_diff(records, ctg);
//             {
//                 std::lock_guard<std::mutex> lock(mtx);
//                 for(auto dif: diff) {
//                     for(auto pos: dif.second) {
//                         vars_[dif.first].emplace_back(pos >> 2);
//                         if(opt_.dump_var) {
//                             out << dif.first << "\t" << (pos >> 2) << "\t" << "MID"[pos & 0x3] << "\n";
//                         }
//                     }
//                 }
//             }
//         }
//     };
//     MultiThreads(opt_.threads, func);
//     if(opt_.dump_var) out.close();
// }

void Phase::GenerateVariantsWithClair3() {
    // std::unordered_map<std::string, std::vector<long>> variants;
    std::unordered_map<std::string, std::unordered_map<long, std::string>> variants;
    // load variants in primary
    GzFileReader preader(opt_.pcfname);
    if(preader.Valid()) {
        std::string line = preader.GetNoEmptyLine();
        while(!line.empty()) {
            if(line[0] != '#') {
                auto items = SplitString(line, '\t');
                if(items[3].size() == 1 && items[4].size() == 1) {
                    // uint8_t ref = 0x03 & ((items[3][0] >> 2) ^ (items[3][0] >> 1));
                    // uint8_t alt = 0x03 & ((items[4][0] >> 2) ^ (items[4][0] >> 1));
                    long pos = atol(items[1].c_str()) - 1;
                    // variants[items[0]].emplace_back((pos << 4 | ref << 2 | alt));
                    variants[items[0]][pos] = items[3] + items[4];
                }
            }
            line = preader.GetNoEmptyLine();
        }
    } else {
        std::cerr << "[" << GetCurTime() << "] Failed to open " << opt_.pcfname << " for reading\n";
        exit(EXIT_FAILURE);
    }
    // load variants in alt
    GzFileReader areader(opt_.acfname);
    if(areader.Valid()) {
        std::string line = areader.GetNoEmptyLine();
        while(!line.empty()) {
            if(line[0] != '#') {
                auto items = SplitString(line, '\t');
                if(items[3].size() == 1 && items[4].size() == 1) {
                    // uint8_t ref = 0x03 & ((items[3][0] >> 2) ^ (items[3][0] >> 1));
                    // uint8_t alt = 0x03 & ((items[4][0] >> 2) ^ (items[4][0] >> 1));
                    long pos = atol(items[1].c_str()) - 1;
                    // variants[items[0]].emplace_back((pos << 4 | ref << 2 | alt));
                    variants[items[0]][pos] = items[4] + items[3];
                }
            }
            line = areader.GetNoEmptyLine();
        }
    } else {
        std::cerr << "[" << GetCurTime() << "] Failed to open " << opt_.acfname << " for reading\n";
        exit(EXIT_FAILURE);
    }
    // load paf
    std::unordered_map<std::string, std::vector<std::array<std::string, 13>>> paf;
    std::ifstream in(opt_.bfname, std::ios::in);
    if(!in.is_open()) {
        std::cerr << "[" << GetCurTime() << "] Error could not open " << opt_.bfname << " for reading\n";
        exit(EXIT_FAILURE);
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
                        if(record[4] == "+") {
                            // if(std::find(variants[ctg].begin(), variants[ctg].end(), ref_pos + i) != variants[ctg].end() && 
                            //     std::find(variants[query].begin(), variants[query].end(), query_pos + i) != variants[query].end()) {
                            //         diff[ctg].emplace_back(ref_pos + i);
                            //         diff[query].emplace_back(query_pos + i);
                            //     }
                            if(variants.find(ctg) != variants.end() && 
                                variants.find(query) != variants.end() && 
                                variants[ctg].find(ref_pos + i) != variants[ctg].end() && 
                                variants[query].find(query_pos + i) != variants[query].end() && 
                                variants[ctg][ref_pos + i] == variants[query][query_pos + i]) {
                                    diff[ctg].emplace_back(ref_pos + i);
                                    diff[query].emplace_back(query_pos + i);
                            }
                        } else {
                            // if(std::find(variants[ctg].begin(), variants[ctg].end(), ref_pos + i) != variants[ctg].end() && 
                            //     std::find(variants[query].begin(), variants[query].end(), query_pos - i) != variants[query].end()) {
                            //         diff[ctg].emplace_back(ref_pos + i);
                            //         diff[query].emplace_back(query_pos - i);
                            //     }
                            if(variants.find(ctg) != variants.end() && 
                                variants.find(query) != variants.end() && 
                                variants[ctg].find(ref_pos + i) != variants[ctg].end() && 
                                variants[query].find(query_pos - i) != variants[query].end() && 
                                variants[ctg][ref_pos + i] == variants[query][query_pos - i]) {
                                    diff[ctg].emplace_back(ref_pos + i);
                                    diff[query].emplace_back(query_pos - i);
                            }
                        }
                        // diff[ctg].emplace_back((ref_pos + i));
                        // if(record[4] == "+") diff[query].emplace_back((query_pos + i));
                        // else diff[query].emplace_back((query_pos - i));
                    }
                    ref_pos += len;
                    if(record[4] == "+") query_pos += len;
                    else query_pos -= len;
                } else if(op == "D") {
                    // for(int i = 0; i < len; ++i) {
                    //     diff[ctg].emplace_back((ref_pos + i));
                    // }
                    // if(record[4] == "+") diff[query].emplace_back((query_pos + 1));
                    // else diff[query].emplace_back((query_pos - 1));
                    ref_pos += len;
                } else if(op == "I") {
                    // for(int i = 0; i < len; ++i) {
                    //     if(record[4] == "+") diff[query].emplace_back((query_pos + i));
                    //     else diff[query].emplace_back((query_pos - i));
                    // }
                    // diff[ctg].emplace_back((ref_pos + 1));
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
                for(auto dif: diff) {
                    for(auto pos: dif.second) {
                        vars_[dif.first].emplace_back(pos);
                    }
                }
            }
        }
    };
    MultiThreads(opt_.threads, func);

    if(opt_.dump_var) {
        std::ofstream out(opt_.prefix + ".var.txt");
        if(!out.is_open()) {
            std::cerr << "[" << GetCurTime() << "] Warning: could not open " << opt_.prefix << ".var.txt for writing\n";
        }
        keys.clear();
        for(auto ctg: vars_) keys.emplace_back(ctg.first);
        std::sort(keys.begin(), keys.end());
        for(auto ctg: keys) {
            for(auto &var: vars_[ctg]) {
                out << ctg << "\t" << var << "\n";
            }
        }
        out.close();
    }
}

void Phase::GenerateVariants() {
    // load paf
    std::unordered_map<std::string, std::vector<std::array<std::string, 13>>> paf;
    std::ifstream in(opt_.bfname, std::ios::in);
    if(!in.is_open()) {
        std::cerr << "[" << GetCurTime() << "] Error could not open " << opt_.bfname << " for reading\n";
        exit(EXIT_FAILURE);
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

    auto find_diff = [](std::vector<std::array<std::string, 13>> &records, std::string &ctg) -> std::unordered_map<std::string, std::vector<uint64_t>> {
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
                        diff[ctg].emplace_back((ref_pos + i) << 2 | 0);
                        if(record[4] == "+") diff[query].emplace_back((query_pos + i) << 2 | 0);
                        else diff[query].emplace_back((query_pos - i) << 2 | 0);
                    }
                    ref_pos += len;
                    if(record[4] == "+") query_pos += len;
                    else query_pos -= len;
                } else if(op == "D") {
                    for(int i = 0; i < len; ++i) {
                        diff[ctg].emplace_back((ref_pos + i) << 2 | 2);
                    }
                    if(record[4] == "+") diff[query].emplace_back((query_pos + 1) << 2 | 1);
                    else diff[query].emplace_back((query_pos - 1) << 2 | 1);
                    ref_pos += len;
                } else if(op == "I") {
                    for(int i = 0; i < len; ++i) {
                        if(record[4] == "+") diff[query].emplace_back((query_pos + i) << 2 | 1);
                        else diff[query].emplace_back((query_pos - i) << 2 | 1);
                    }
                    diff[ctg].emplace_back((ref_pos + 1) << 2 | 1);
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
                for(auto dif: diff) {
                    for(auto pos: dif.second) {
                        vars_[dif.first].emplace_back(pos >> 2);
                        // if(opt_.dump_var) {
                        //     out << dif.first << "\t" << (pos >> 2) << "\t" << "MID"[pos & 0x3] << "\n";
                        // }
                    }
                }
            }
        }
    };
    MultiThreads(opt_.threads, func);
    
    if(opt_.dump_var) {
        std::ofstream out(opt_.prefix + ".var.txt");
        if(!out.is_open()) {
            std::cerr << "[" << GetCurTime() << "] Warning: could not open " << opt_.prefix << ".var.txt for writing\n";
        }
        keys.clear();
        for(auto ctg: vars_) keys.emplace_back(ctg.first);
        std::sort(keys.begin(), keys.end());
        for(auto ctg: keys) {
            for(auto &var: vars_[ctg]) {
                // out << ctg << "\t" << (var >> 2) << "\t" << "MID"[var & 0x3] << "\n";
                out << ctg << "\t" << var << "\n";
            }
        }
        out.close();
    }

    // for(auto &ctg: vars_) {
    //     for(auto &pos: ctg.second) {
    //         pos = (pos >> 2);
    //     }
    // }
}

void Phase::LoadVariants() {
    std::string line;
    std::ifstream in(opt_.vfname);
    if(!in.is_open()) {
        std::cerr << "[" << GetCurTime() << "] Error could not open " << opt_.vfname << " for reading\n";
        exit(EXIT_FAILURE);
    }
    while(std::getline(in, line)) {
        auto items = SplitString(line, '\t');
        vars_[items[0]].emplace_back(atol(items[1].c_str()));
    }
    in.close();
}

std::unordered_map<std::string, std::string> Phase::LoadBed() const {
    std::unordered_map<std::string, std::string> bed;
    std::ifstream in1(opt_.b1fname);
    if(!in1.is_open()) {
        std::cerr << "[" << GetCurTime() << "] Couldn't open file " << opt_.b1fname << " for reading\n";
        exit(EXIT_FAILURE);
    }
    std::string line;
    while(std::getline(in1, line)) {
        auto items = SplitString(line, '\t');
        assert(items.size() == 4);
        std::string key = items[0] + ':' + items[1] + '-' + items[2];
        bed[key] = items[3];
    }
    std::ifstream in2(opt_.b2fname);
    if(!in2.is_open()) {
        std::cerr << "[" << GetCurTime() << "] Couldn't open file " << opt_.b2fname << " for reading\n";
        exit(EXIT_FAILURE);
    }
    while(std::getline(in2, line)) {
        auto items = SplitString(line, '\t');
        assert(items.size() == 4);
        std::string key = items[0] + ':' + items[1] + '-' + items[2];
        bed[key] = items[3];
    }
    return bed;
}

void Phase::GenerateMatrixScaffold(std::unordered_map<std::string, std::string> &bed) {
    reader_.Initialize(opt_.hfname);
    matrix_.Resize(name2id_.size());
    std::array<bam1_t*, 2> pair;
    pair[0] = bam_init1();
    pair[1] = bam_init1();
    long count = 0;
    while(1) {
        int r = reader_.LoadPair(pair);
        if(r < 0) break;
        std::string name1 = bam_get_qname(pair[0]);
        std::string name2 = bam_get_qname(pair[1]);
        // if(pair[0]->core.qual < 10 || pair[1]->core.qual < 10) continue;
        // else {
        //     uint8_t *nm1 = bam_aux_get(pair[0], "NM");
        //     if(nm1 != NULL && bam_aux2i(nm1) > 5) continue;
        //     uint8_t *nm2 = bam_aux_get(pair[1], "NM");
        //     if(nm2 != NULL && bam_aux2i(nm2) > 5) continue;
        // }
        // if(name1 != name2) {
        //     std::cerr << "[" << GetCurTime() << "] Read pair " << name1 << " != " << name2 << "\n";
        //     exit(EXIT_FAILURE);
        // }
        std::size_t id0 = name2id_[bed[reader_.Header()->target_name[pair[0]->core.tid]]];
        std::size_t id1 = name2id_[bed[reader_.Header()->target_name[pair[1]->core.tid]]];
        if(matrix_.AddLink(id0, id1) < 1) {
            std::cerr << "[" << GetCurTime() << "] Error in AddLink\n";
            exit(EXIT_FAILURE);
        }
        if(matrix_.AddLink(id1, id0) < 1) {
            std::cerr << "[" << GetCurTime() << "] Error in AddLink\n";
            exit(EXIT_FAILURE);
        }
        count += 1;
        if((count % 1000000) == 0) {
            std::cerr << "[" << GetCurTime() << "] Parsed " << count << " read pairs\n";
        }
    }
    bam_destroy1(pair[0]);
    bam_destroy1(pair[1]);
}

// void Phase::IncreaseCov(std::string &name, std::array<bam1_t*, 2> &pair) {
//     std::size_t pos1 = pair[0]->core.pos;
//     std::size_t pos2 = pair[1]->core.pos;
//     if(pos1 < pos2) {
//         pos2 += bam_cigar2rlen(pair[1]->core.n_cigar, bam_get_cigar(pair[1]));
//         cov_[name][pos1] += 1;
//         cov_[name][pos2] -= 1;
//     } else {
//         pos1 += bam_cigar2rlen(pair[0]->core.n_cigar, bam_get_cigar(pair[0]));
//         cov_[name][pos2] += 1;
//         cov_[name][pos1] -= 1;
//     }
// }

void Phase::IncreaseCov(std::string &name, bam1_t *a, bam1_t *b) {
    std::size_t pos1 = a->core.pos;
    std::size_t pos2 = b->core.pos;
    if(pos1 < pos2) {
        pos2 = bam_endpos(b);
        cov_[name][pos1] += 1;
        cov_[name][pos2] -= 1;
    } else {
        pos1 = bam_endpos(a);
        cov_[name][pos2] += 1;
        cov_[name][pos1] -= 1;
    }
}

std::unordered_map<std::string, std::size_t> Phase::DetectSwitchError() {
    std::unordered_map<std::string, std::size_t> positions;
    std::unordered_map<std::string, std::array<size_t, 2>> regions;
    std::unordered_map<std::string, std::vector<std::array<int, 2>>> stats;
    std::regex pattern("(\\w+):(\\d+)-(\\d+)");
    std::smatch result;
    for(auto ctg: cov_) {
        auto cov = ctg.second;
        for(std::size_t i = 1; i < cov.size(); ++i) {
            cov[i] = cov[i] + cov[i - 1];
        }
        if(*std::max_element(cov.begin(), cov.end()) < 30) continue;
        double average = std::accumulate(cov.begin(), cov.end(), 0.0) / cov.size();
        std::vector<std::size_t> pos;
        for(int div = 5; div <= 15; ++div) {
            double cutoff = average / div;
            std::vector<int> delta;
            for(std::size_t i = 0; i < cov.size(); ++i) {
                if(cov[i] < cutoff) {
                    delta.emplace_back(10);
                } else {
                    delta.emplace_back(-10);
                }
            }
            std::size_t sz = delta.size();
            int max_so_far = -100000, max_ending_here = 0;
            std::size_t start = sz / 1000, end = sz / 1000, s = sz / 1000;
            for(std::size_t i = sz / 1000; i < 99 * sz / 100; ++i) {
                max_ending_here += delta.at(i);
                if(max_so_far < max_ending_here) {
                    max_so_far = max_ending_here;
                    start = s;
                    end = i;
                }
                if(max_ending_here < 0) {
                    max_ending_here = 0;
                    s = i + 1;
                }
            }
            if(start >= sz / 1000 + 50000 && end <= 99 * sz / 100 - 50000) {
                pos.emplace_back((start + end) / 2);
            }
        }
        if(pos.size() > 0) {
            std::size_t consensus_pos = std::accumulate(pos.begin(), pos.end(), 0) / pos.size();
            positions[ctg.first] = consensus_pos;
            std::size_t st, en;
            for(st = consensus_pos; st >= 0; --st) {
                if(cov[st] > cov[consensus_pos] + 10) break;
            }
            for(en = consensus_pos; en < cov.size(); ++en) {
                if(cov[en] > cov[consensus_pos] + 10) break;
            }
            std::string name = ctg.first;
            int offset = 0;
            if(std::regex_search(ctg.first, result, pattern)) {
                name = result[1].str();
                offset = atoi(result[2].str().c_str());
            }
            if(vars_.find(name) != vars_.end()) {
                auto left = std::lower_bound(vars_[name].begin(), vars_[name].end(), st + offset);
                auto right = std::upper_bound(vars_[name].begin(), vars_[name].end(), en + offset);
                if(right - left > 1) {
                    regions[ctg.first] = std::array<std::size_t, 2> { st, en };
                    stats[ctg.first].resize(right - left);
                }
            }
        }
    }
    return positions;
/*   
    BamReader reader;
    reader.Initialize(opt_.prefix + ".hic.filtered.bam");
    std::array<bam1_t*, 2> pair;
    pair[0] = bam_init1();
    pair[1] = bam_init1();
    while(1) {
        int ret = reader.LoadPair(pair);
        if(ret < 0) break;
        if(pair[0]->core.tid != pair[1]->core.tid) continue;
        if(FilterHiCPairWithSNP(pair) != 0) continue;
        std::string rname = reader.Header()->target_name[pair[0]->core.tid];
        if(regions.find(rname) == regions.end()) break;
        std::string name = rname;
        int offset = 0;
        if(std::regex_search(name, result, pattern)) {
            name = result[1].str();
            offset = atoi(result[2].str().c_str());
        }
        std::size_t st1 = pair[0]->core.pos;
        std::size_t en1 = st1 + bam_cigar2rlen(pair[0]->core.n_cigar, bam_get_cigar(pair[0]));
        std::size_t st2 = pair[1]->core.pos;
        std::size_t en2 = st2 + bam_cigar2rlen(pair[1]->core.n_cigar, bam_get_cigar(pair[1]));
        if(st1 >= regions[rname][0] && en1 <= regions[rname][1]) {
            if(st2 < regions[rname][0] || en2 > regions[rname][1]) {
                auto lb = std::lower_bound(vars_[name].begin(), vars_[name].end(), regions[rname][0] + offset);
                auto left = std::lower_bound(vars_[name].begin(), vars_[name].end(), st1 + offset) - lb;
                auto right = std::upper_bound(vars_[name].begin(), vars_[name].end(), en1 + offset) - lb;
                for(auto i = left; i < right; ++i) {
                    if(st2 < regions[rname][0]) stats[rname][i][0] += 1;
                    else stats[rname][i][1] += 1;
                }
            }
        } else if(st2 >= regions[rname][0] && en2 <= regions[rname][1]) {
            if(st1 < regions[rname][0] || en1 > regions[rname][1]) {
                auto lb = std::lower_bound(vars_[name].begin(), vars_[name].end(), regions[rname][0] + offset);
                auto left = std::lower_bound(vars_[name].begin(), vars_[name].end(), st2 + offset) - lb;
                auto right = std::upper_bound(vars_[name].begin(), vars_[name].end(), en2 + offset) - lb;
                for(auto i = left; i < right; ++i) {
                    if(st1 < regions[rname][0]) stats[rname][i][0] += 1;
                    else stats[rname][i][1] += 1;
                }
            }
        }
    }
    bam_destroy1(pair[0]);
    bam_destroy1(pair[1]);
    for(auto ctg: stats) {
        std::string name = ctg.first;
        int offset = 0;
        if(std::regex_search(name, result, pattern)) {
            name = result[1].str();
            offset = atoi(result[2].str().c_str());
        }
        auto left = std::lower_bound(vars_[name].begin(), vars_[name].end(), regions[ctg.first][0] + offset) - vars_[name].begin();
        auto right = std::upper_bound(vars_[name].begin(), vars_[name].end(), regions[ctg.first][1] + offset) - vars_[name].begin();
        for(auto j = left; j < right; ++j) {
            std::cerr << vars_[name][j] << " ";
        }
        std::cerr << "\n";
        for(auto c: ctg.second) {
            std::cerr << "[" << c[0] << ", " << c[1] << "] ";
        }
        std::cerr << "\n";
        std::size_t i = 0;
        for(; i < ctg.second.size(); ++i) {
            if(ctg.second[i][0] < ctg.second[i][1]) break;
        }
        auto offsets = std::lower_bound(vars_[name].begin(), vars_[name].end(), regions[ctg.first][0] + offset) - vars_[name].begin();
        std::cerr << "switch error " << ctg.first << " " << vars_[name][offsets + i] << "-" << vars_[name][offsets + i + 1] << "\n";
    }
    return positions;
*/ 
}

void Phase::FilterSwitchError(std::unordered_map<std::string, std::size_t> &se, std::vector<std::array<std::string, 2>> &pair) {
    for(auto p: pair) {
        auto iter1 = se.find(p[0]);
        auto iter2 = se.find(p[1]);
        if(iter1 != se.end() && iter2 != se.end()) {
            switch_error_[iter1->first] = iter1->second;
            switch_error_[iter2->first] = iter2->second;
        }
    }
}
/*
void Phase::FixPhase(std::string &fname) {
    BamReader reader;
    reader.Initialize(fname);
    std::unordered_map<std::string, std::array<std::size_t, 2>> count;
    auto find_func = [&](std::array<bam1_t*, 2> &pair) -> int {
        if(pair[0]->core.tid == pair[1]->core.tid) return 0;
        std::string rname1 = reader.Header()->target_name[pair[0]->core.tid];
        std::string rname2 = reader.Header()->target_name[pair[1]->core.tid];
        if(switch_error_.find(rname1) != switch_error_.end()) return 1;
        if(switch_error_.find(rname2) != switch_error_.end()) return 2;
        return 0;
    };
    std::size_t SIZE = 100000;
    std::vector<std::array<bam1_t*, 2>> reads(SIZE);
    for(std::size_t i = 0; i < SIZE; ++i) {
        reads[i][0] = bam_init1();
        reads[i][1] = bam_init1();
    }
    while(1) {
        std::size_t size = 0;
        for(std::size_t i = 0; i < SIZE; ++i) {
            int r = reader.LoadPair(reads[i]);
            if(r < 0) break;
            size += 1;
        }
        if(size == 0) break;
        // std::cerr << "size " << size << "\n";
        int flag[size];
        std::atomic<std::size_t> index { 0 };
        auto func = [&](std::size_t tid) {
            for(auto cur = index.fetch_add(1); cur < size; cur = index.fetch_add(1)) {
                flag[cur] = find_func(reads[cur]);
            }
        };
        MultiThreads(opt_.threads, func);
        for(std::size_t i = 0; i < size; ++i) {
            if(flag[i] == 1) {
                std::string rname = reader.Header()->target_name[reads[i][0]->core.tid];
                std::size_t pos = reads[i][0]->core.pos;
                if(pos < switch_error_[rname]) count[rname][0] += 1;
                else count[rname][1] += 1;
            } else if(flag[i] == 2) {
                std::string rname = reader.Header()->target_name[reads[i][1]->core.tid];
                std::size_t pos = reads[i][1]->core.pos;
                if(pos < switch_error_[rname]) count[rname][0] += 1;
                else count[rname][1] += 1;
            }
        }
    }
    for(std::size_t i = 0; i < SIZE; ++i) {
        bam_destroy1(reads[i][0]);
        bam_destroy1(reads[i][1]);
    }

    for(auto ctg: count) {
        std::cerr << "counts " << ctg.first << " " << ctg.second[0] << " " << ctg.second[1] << std::endl;
        if(ctg.second[0] < ctg.second[1]) consistency_[ctg.first] = 'F';
        else consistency_[ctg.first] = 'S';
    }
}
*/
void Phase::FixPhase(std::string &fname) {
    BamReader reader;
    reader.Initialize(fname);
    std::array<bam1_t*, 2> pair;
    pair[0] = bam_init1();
    pair[1] = bam_init1();
    while(1) {
        int r = reader.LoadPair(pair);
        if(r < 0) break;
        if(pair[0]->core.tid == pair[1]->core.tid) continue;
        if(FilterHiCPairWithSNP(pair) != 0) continue;
        std::string name1 = reader.Header()->target_name[pair[0]->core.tid];
        std::string name2 = reader.Header()->target_name[pair[1]->core.tid];
        auto iter1 = switch_error_.find(name1);
        auto iter2 = switch_error_.find(name2);
        if(iter1 != switch_error_.end()) {
            std::size_t index1 = name2id_[name1];
            std::size_t index2 = name2id_[name2];
            std::size_t pos = pair[0]->core.pos;
            if(pos < iter1->second) {
                // first segment
                matrix_.AddLink(index1, index2);
            } else {
                // last segment
                matrix_.AddLink(index2, index1);
            }
        } else if(iter2 != switch_error_.end()) {
            std::size_t index1 = name2id_[name1];
            std::size_t index2 = name2id_[name2];
            std::size_t pos = pair[1]->core.pos;
            if(pos < iter2->second) {
                // first segment
                matrix_.AddLink(index2, index1);
            } else {
                // last segment
                matrix_.AddLink(index1, index2);
            }
        }
    }
    bam_destroy1(pair[0]);
    bam_destroy1(pair[1]);
    std::unordered_set<std::string> cand_ctg;
    for(auto ctg: switch_error_) {
        cand_ctg.insert(ctg.first.substr(0,7));
    }
    for(auto ctg: contigs_) {
        if(cand_ctg.find(ctg.name) == cand_ctg.end()) continue;
        for(std::size_t i = 0; i < ctg.olps.size(); ++i) {
            std::string name = names_[ctg.olps[i].first];
            if(switch_error_.find(name) == switch_error_.end()) continue;
            Print(ctg.olps);
            std::array<int, 4> count = { 0, 0, 0, 0 };
            for(std::size_t j = 0; j < ctg.olps.size(); ++j) {
                if(j == i) continue;
                int c1 = (int)matrix_.Get(ctg.olps[i].first, ctg.olps[j].first);
                int c2 = (int)matrix_.Get(ctg.olps[i].first, ctg.olps[j].second);
                int c3 = (int)matrix_.Get(ctg.olps[i].second, ctg.olps[j].first);
                int c4 = (int)matrix_.Get(ctg.olps[i].second, ctg.olps[j].second);
                int c5 = (int)matrix_.Get(ctg.olps[j].first, ctg.olps[i].first);
                int c6 = (int)matrix_.Get(ctg.olps[j].first, ctg.olps[i].second);
                int c7 = (int)matrix_.Get(ctg.olps[j].second, ctg.olps[i].first);
                int c8 = (int)matrix_.Get(ctg.olps[j].second, ctg.olps[i].second);
                if(ctg.olps[i].label == ctg.olps[j].label) {
                    if(c1 + c4 > c2 + c3) count[0] += 1;
                    else if(c1 + c4 < c2 + c3) count[1] += 1;
                    if(c5 + c8 > c6 + c7) count[2] += 1;
                    else if(c5 + c8 < c6 + c7) count[3] += 1;
                } else {
                    if(c1 + c4 > c2 + c3) count[1] += 1;
                    else if(c1 + c4 < c2 + c3) count[0] += 1;
                    if(c5 + c8 > c6 + c7) count[3] += 1;
                    else if(c5 + c8 < c6 + c7) count[2] += 1;
                }
            }
            if(count[0] > count[1] && count[2] < count[3]) {
                consistency_[name] = 'F';
                consistency_[names_[ctg.olps[i].second]] = 'F';
            } else if(count[0] < count[1] && count[2] > count[3]) {
                consistency_[name] = 'S';
                consistency_[names_[ctg.olps[i].second]] = 'S';
            } else if(count[0] < count[1] && count[2] < count[3]) {
                consistency_[name] = 'C';
                consistency_[names_[ctg.olps[i].second]] = 'C';
                // std::cerr << name << " " << count[0] << " " << count[1] << " " << count[2] << " " << count[3] << "\n";
            }
        }
    }
    // for(auto c: consistency_) {
    //     std::cerr << c.first << " " << c.second << "\n";
    // }
}

std::vector<std::size_t> Phase::FilterPoreC(std::vector<bam1_t*> &reads) {
    std::vector<std::size_t> ret;
    std::regex pattern("(\\w+):(\\d+)-(\\d+)");
    std::smatch result;

    for(std::size_t i = 0; i < reads.size(); ++i) {
        std::string name = reader_.Header()->target_name[reads[i]->core.tid];
        long pos = reads[i]->core.pos;
        int beg = 0;
        if(std::regex_search(name, result, pattern)) {
            name = result[1].str();
            beg = atoi(result[2].str().c_str());
        }
        if(reads[i]->core.qual >= opt_.mapq) {
            uint8_t *nm = bam_aux_get(reads[i], "NM");
            if(nm == NULL || bam_aux2i(nm) < 5) {
                ret.emplace_back(i);
            }
        } else if(vars_.find(name) != vars_.end()) {
            int endpos = bam_endpos(reads[i]);
            auto lb = std::lower_bound(vars_[name].begin(), vars_[name].end(), pos + beg);
            if(lb != vars_[name].end() && *lb > pos + beg && *lb < endpos + beg) ret.emplace_back(i);
        }
    }
    
    return ret;
}

std::vector<std::size_t> Phase::FilterPoreCWithSNP(std::vector<bam1_t*> &reads) {
    std::vector<std::size_t> ret;
    std::regex pattern("(\\w+):(\\d+)-(\\d+)");
    std::smatch result;

    for(std::size_t i = 0; i < reads.size(); ++i) {
        std::string name = reader_.Header()->target_name[reads[i]->core.tid];
        long pos = reads[i]->core.pos;
        int beg = 0;
        if(std::regex_search(name, result, pattern)) {
            name = result[1].str();
            beg = atoi(result[2].str().c_str());
        }
        if(vars_.find(name) != vars_.end()) {
            int endpos = bam_endpos(reads[i]);
            auto lb = std::lower_bound(vars_[name].begin(), vars_[name].end(), pos + beg);
            if(lb != vars_[name].end() && *lb > pos + beg && *lb < endpos + beg) ret.emplace_back(i);
        }
    }

    return ret;
}

void Phase::FilterAndGenerateMatrixPoreC() {
    reader_.Initialize(opt_.hfname);
    for(int32_t i = 0; i < reader_.Header()->n_targets; ++i) {
        cov_[reader_.Header()->target_name[i]].resize(reader_.Header()->target_len[i] + 1, 0);
    }
    htsFile *out;
    if(opt_.dump_filtered) {
        std::string ofname = opt_.prefix + ".hic.filtered.bam";
        out = hts_open(ofname.c_str(), "wb");
        if(sam_hdr_write(out, reader_.Header()) != 0) {
            std::cerr << "[" << GetCurTime() << "] Error failed to write header to " << ofname << "\n";
            exit(EXIT_FAILURE);
        }
    }
    matrix_.Resize(reader_.Header()->n_targets);
    std::size_t SIZE = 100000;
    std::vector<std::vector<bam1_t*>> reads(SIZE);
    long count = 0, filtered = 0, snp = 0;
    while(1) {
        std::size_t size = reader_.LoadBatchReads(reads);
        if(size == 0) break;
        // std::cerr << "size " << size << "\n";
        std::vector<std::vector<std::size_t>> flag(size);
        std::atomic<std::size_t> index { 0 };
        auto func = [&](std::size_t tid) {
            for(auto cur = index.fetch_add(1); cur < size; cur = index.fetch_add(1)) {
                flag[cur] = FilterPoreC(reads[cur]);
            }
        };
        MultiThreads(opt_.threads, func);
        // std::cerr << "filter" << std::endl;
        for(std::size_t i = 0; i < size; ++i) {
            for(std::size_t j = 0; j < flag[i].size(); ++j) {
                std::string name1 = bam_get_qname(reads[i][j]);
                for(std::size_t k = j + 1; k < flag[i].size(); ++k) {
                    std::string name2 = bam_get_qname(reads[i][k]);
                    if(name1 != name2) {
                        std::cerr << "[" << GetCurTime() << "] PoreC loaded wrong " << name1 << " != " << name2 << "\n";
                        exit(EXIT_FAILURE);
                    }
                    if(matrix_.AddLink(reads[i][flag[i][j]]->core.tid, reads[i][flag[i][k]]->core.tid) < 1) {
                        std::cerr << "[" << GetCurTime() << "] Error in AddLink\n";
                        exit(EXIT_FAILURE);
                    }
                    if(matrix_.AddLink(reads[i][flag[i][k]]->core.tid, reads[i][flag[i][j]]->core.tid) < 1) {
                        std::cerr << "[" << GetCurTime() << "] Error in AddLink\n";
                        exit(EXIT_FAILURE);
                    }
                    if(reads[i][flag[i][j]]->core.tid == reads[i][flag[i][k]]->core.tid) {
                        std::string name = reader_.Header()->target_name[reads[i][flag[i][j]]->core.tid];
                        IncreaseCov(name, reads[i][flag[i][j]], reads[i][flag[i][k]]);
                    }
                }
            }
            if(opt_.dump_filtered) {
                for(auto f: flag[i]) {
                    sam_write1(out, reader_.Header(), reads[i][f]);
                }
            }
        }
        // std::cerr << "dumpped" << std::endl;
        count += size;
        if((count % 1000000) == 0) {
            std::cerr << "[" << GetCurTime() << "] Parsed " << count << " reads\n";
        }

        for(std::size_t i = 0; i < size; ++i) {
            for(auto &r: reads[i]) bam_destroy1(r);
            reads[i].clear();
        }
        // std::cerr << "destroyed" << std::endl;
    }
    std::cerr << "[" << GetCurTime() << "] Total " << count << " read pairs\n";
    if(opt_.dump_filtered) hts_close(out);
}

void Phase::GenerateMatrixPoreC() {
    reader_.Initialize(opt_.hfname);
    for(int32_t i = 0; i < reader_.Header()->n_targets; ++i) {
        cov_[reader_.Header()->target_name[i]].resize(reader_.Header()->target_len[i] + 1, 0);
    }
    matrix_.Resize(reader_.Header()->n_targets);
    std::size_t SIZE = 100000;
    std::vector<std::vector<bam1_t*>> reads(SIZE);
    long count = 0;
    while(1) {
        int size = reader_.LoadBatchReads(reads);
        if(size < 0) break;
        for(std::size_t i = 0; i < size; ++i) {
            for(std::size_t j = 0; j < reads[i].size(); ++j) {
                std::string name1 = bam_get_qname(reads[i][j]);
                for(std::size_t k = j + 1; k < reads[i].size(); ++k) {
                    std::string name2 = bam_get_qname(reads[i][k]);
                    if(name1 != name2) {
                        std::cerr << "[" << GetCurTime() << "] PoreC loaded wrong " << name1 << " != " << name2 << "\n";
                        exit(EXIT_FAILURE);
                    }
                    if(matrix_.AddLink(reads[i][j]->core.tid, reads[i][k]->core.tid) < 1) {
                        std::cerr << "[" << GetCurTime() << "] Error in AddLink\n";
                        exit(EXIT_FAILURE);
                    }
                    if(matrix_.AddLink(reads[i][k]->core.tid, reads[i][j]->core.tid) < 1) {
                        std::cerr << "[" << GetCurTime() << "] Error in AddLink\n";
                        exit(EXIT_FAILURE);
                    }
                    if(reads[i][j]->core.tid == reads[i][k]->core.tid) {
                        std::string name = reader_.Header()->target_name[reads[i][j]->core.tid];
                        IncreaseCov(name, reads[i][j], reads[i][k]);
                    }
                }
            }
        }
        count += size;
        if((count % 1000000) == 0) {
            std::cerr << "[" << GetCurTime() << "] Parsed " << count << " read pairs\n";
        }
        for(std::size_t i = 0; i < size; ++i) {
            for(auto r: reads[i]) bam_destroy1(r);
        }
    }
}

void Phase::FixPhasePoreC(std::string &fname) {
    BamReader reader;
    reader.Initialize(fname);
    std::size_t SIZE = 100000;
    std::vector<std::vector<bam1_t*>> reads(SIZE);
    while(1) {
        std::size_t size = reader.LoadBatchReads(reads);
        if(size == 0) break;
        for(auto read: reads) {
            auto reads_snp = FilterPoreCWithSNP(read);
            for(std::size_t i = 0; i < reads_snp.size(); ++i) {
                std::string name1 = bam_get_qname(read[reads_snp[i]]);
                for(std::size_t j = i + 1; j < reads_snp.size(); ++j) {
                    if(read[reads_snp[j]]->core.tid == read[reads_snp[i]]->core.tid) continue;
                    std::string name2 = bam_get_qname(read[reads_snp[j]]);
                    if(name1 != name2) {
                        std::cerr << "[" << GetCurTime() << "] PoreC loaded wrong " << name1 << " != " << name2 << "\n";
                        exit(EXIT_FAILURE);
                    }
                    std::string rname1 = reader.Header()->target_name[read[reads_snp[i]]->core.tid];
                    std::string rname2 = reader.Header()->target_name[read[reads_snp[j]]->core.tid];
                    if(switch_error_.find(rname1) != switch_error_.end()) {
                        std::size_t index1 = name2id_[rname1];
                        std::size_t index2 = name2id_[rname2];
                        std::size_t pos = read[reads_snp[i]]->core.pos;
                        if(pos < switch_error_[rname1]) {
                            // first segment
                            matrix_.AddLink(index1, index2);
                        } else {
                            // last segment
                            matrix_.AddLink(index2, index1);
                        }
                    } else if(switch_error_.find(rname2) != switch_error_.end()) {
                        std::size_t index1 = name2id_[rname1];
                        std::size_t index2 = name2id_[rname2];
                        std::size_t pos = read[reads_snp[j]]->core.pos;
                        if(pos < switch_error_[rname2]) {
                            // first segment
                            matrix_.AddLink(index2, index1);
                        } else {
                            // last segment
                            matrix_.AddLink(index1, index2);
                        }
                    }
                }
            }
        }
        for(std::size_t i = 0; i < size; ++i) {
            for(auto r: reads[i]) bam_destroy1(r);
        }
    }
    std::unordered_set<std::string> cand_ctg;
    for(auto ctg: switch_error_) {
        cand_ctg.insert(ctg.first.substr(0,7));
    }
    for(auto ctg: contigs_) {
        if(cand_ctg.find(ctg.name) == cand_ctg.end()) continue;
        for(std::size_t i = 0; i < ctg.olps.size(); ++i) {
            std::string name = names_[ctg.olps[i].first];
            if(switch_error_.find(name) == switch_error_.end()) continue;
            Print(ctg.olps);
            std::array<int, 4> count = { 0, 0, 0, 0 };
            for(std::size_t j = 0; j < ctg.olps.size(); ++j) {
                if(j == i) continue;
                int c1 = (int)matrix_.Get(ctg.olps[i].first, ctg.olps[j].first);
                int c2 = (int)matrix_.Get(ctg.olps[i].first, ctg.olps[j].second);
                int c3 = (int)matrix_.Get(ctg.olps[i].second, ctg.olps[j].first);
                int c4 = (int)matrix_.Get(ctg.olps[i].second, ctg.olps[j].second);
                int c5 = (int)matrix_.Get(ctg.olps[j].first, ctg.olps[i].first);
                int c6 = (int)matrix_.Get(ctg.olps[j].first, ctg.olps[i].second);
                int c7 = (int)matrix_.Get(ctg.olps[j].second, ctg.olps[i].first);
                int c8 = (int)matrix_.Get(ctg.olps[j].second, ctg.olps[i].second);
                if(ctg.olps[i].label == ctg.olps[j].label) {
                    if(c1 + c4 > c2 + c3) count[0] += 1;
                    else if(c1 + c4 < c2 + c3) count[1] += 1;
                    if(c5 + c8 > c6 + c7) count[2] += 1;
                    else if(c5 + c8 < c6 + c7) count[3] += 1;
                } else {
                    if(c1 + c4 > c2 + c3) count[1] += 1;
                    else if(c1 + c4 < c2 + c3) count[0] += 1;
                    if(c5 + c8 > c6 + c7) count[3] += 1;
                    else if(c5 + c8 < c6 + c7) count[2] += 1;
                }
            }
            if(count[0] > count[1] && count[2] < count[3]) {
                consistency_[name] = 'F';
                consistency_[names_[ctg.olps[i].second]] = 'F';
            } else if(count[0] < count[1] && count[2] > count[3]) {
                consistency_[name] = 'S';
                consistency_[names_[ctg.olps[i].second]] = 'S';
            } else if(count[0] < count[1] && count[2] < count[3]) {
                consistency_[name] = 'C';
                consistency_[names_[ctg.olps[i].second]] = 'C';
                // std::cerr << name << " " << count[0] << " " << count[1] << " " << count[2] << " " << count[3] << "\n";
            }
        }
    }
}

int Phase::FilterHiCPair(std::array<bam1_t*, 2> &pair) {
    int span = 0;
    std::regex pattern("(\\w+):(\\d+)-(\\d+)");
    std::smatch result;
    
    std::string name1 = reader_.Header()->target_name[pair[0]->core.tid];
    long ref_pos1 = pair[0]->core.pos;
    int beg1 = 0;
    if(std::regex_search(name1, result, pattern)) {
        name1 = result[1].str();
        beg1 = atoi(result[2].str().c_str());
    }
    if(vars_.find(name1) != vars_.end()) {
        int len = bam_cigar2rlen(pair[0]->core.n_cigar, bam_get_cigar(pair[0]));
        auto lb = std::lower_bound(vars_[name1].begin(), vars_[name1].end(), ref_pos1 + beg1);
        if(lb != vars_[name1].end() && *lb > ref_pos1 + beg1 && *lb < ref_pos1 + beg1 + len) span = 2;
    }

    if(span == 0 && pair[0]->core.qual > opt_.mapq) {
        uint8_t *nm = bam_aux_get(pair[0], "NM");
        if(nm == NULL || bam_aux2i(nm) < 5) {
            span = 1;
        }
    } 
    if(span == 0) return 0;
    
    std::string name2 = reader_.Header()->target_name[pair[1]->core.tid];
    long ref_pos2 = pair[1]->core.pos;
    int beg2 = 0;
    if(std::regex_search(name2, result, pattern)) {
        name2 = result[1].str();
        beg2 = atoi(result[2].str().c_str());
    }
    if(vars_.find(name2) != vars_.end()) {
        int len = bam_cigar2rlen(pair[1]->core.n_cigar, bam_get_cigar(pair[1]));
        auto lb = std::lower_bound(vars_[name2].begin(), vars_[name2].end(), ref_pos2 + beg2);
        if(lb != vars_[name2].end() && *lb > ref_pos2 + beg2 && *lb < ref_pos2 + beg2 + len) {
            if(span == 1) return 1;
            else return 2;
        }
    }
  
    if(pair[1]->core.qual > opt_.mapq) {
        uint8_t *nm = bam_aux_get(pair[1], "NM");
        if(nm == NULL || bam_aux2i(nm) < 5) {
            return 1;
        }
    } 
    return 0;
    // std::string name1 = reader_.Header()->target_name[pair[0]->core.tid];
    // std::string name2 = reader_.Header()->target_name[pair[1]->core.tid];
    // std::regex pattern("(\\w+):(\\d+)-(\\d+)");
    // std::smatch result;
    // long ref_pos = pair[0]->core.pos;
    // int beg = 0;
    // if(std::regex_search(name1, result, pattern)) {
    //     name1 = result[1].str();
    //     beg = atoi(result[2].str().c_str());
    // }
    // bool span = false;
    // if(vars_.find(name1) != vars_.end()) {
    //     int len = bam_cigar2rlen(pair[0]->core.n_cigar, bam_get_cigar(pair[0]));
    //     auto lb = std::lower_bound(vars_[name1].begin(), vars_[name1].end(), ref_pos + beg);
    //     if(lb != vars_[name1].end() && *lb > ref_pos + beg && *lb < ref_pos + beg + len) span = true;
    // }
    // ref_pos = pair[1]->core.pos;
    // beg = 0;
    // if(std::regex_search(name2, result, pattern)) {
    //     name2 = result[1].str();
    //     beg = atoi(result[2].str().c_str());
    // }
    // if(span && vars_.find(name2) != vars_.end()) {
    //     int len = bam_cigar2rlen(pair[1]->core.n_cigar, bam_get_cigar(pair[1]));
    //     auto lb = std::lower_bound(vars_[name2].begin(), vars_[name2].end(), ref_pos + beg);
    //     if(lb != vars_[name2].end() && *lb > ref_pos + beg && *lb < ref_pos + beg + len) return 0;
    // }
    // return 1;
}

int Phase::FilterHiCPairWithSNP(std::array<bam1_t*, 2> &pair) {
    std::string name1 = reader_.Header()->target_name[pair[0]->core.tid];
    std::string name2 = reader_.Header()->target_name[pair[1]->core.tid];
    std::regex pattern("(\\w+):(\\d+)-(\\d+)");
    std::smatch result;
    long ref_pos = pair[0]->core.pos;
    int beg = 0;
    if(std::regex_search(name1, result, pattern)) {
        name1 = result[1].str();
        beg = atoi(result[2].str().c_str());
    }
    bool span = false;
    if(vars_.find(name1) != vars_.end()) {
        int len = bam_cigar2rlen(pair[0]->core.n_cigar, bam_get_cigar(pair[0]));
        auto lb = std::lower_bound(vars_[name1].begin(), vars_[name1].end(), ref_pos + beg);
        if(lb != vars_[name1].end() && *lb > ref_pos + beg && *lb < ref_pos + beg + len) span = true;
    }
    ref_pos = pair[1]->core.pos;
    beg = 0;
    if(std::regex_search(name2, result, pattern)) {
        name2 = result[1].str();
        beg = atoi(result[2].str().c_str());
    }
    if(span && vars_.find(name2) != vars_.end()) {
        int len = bam_cigar2rlen(pair[1]->core.n_cigar, bam_get_cigar(pair[1]));
        auto lb = std::lower_bound(vars_[name2].begin(), vars_[name2].end(), ref_pos + beg);
        if(lb != vars_[name2].end() && *lb > ref_pos + beg && *lb < ref_pos + beg + len) return 0;
    }
    return 1;
}

void Phase::FilterAndGenerateMatrix() {
    reader_.Initialize(opt_.hfname);
    for(int32_t i = 0; i < reader_.Header()->n_targets; ++i) {
        cov_[reader_.Header()->target_name[i]].resize(reader_.Header()->target_len[i] + 1, 0);
    }
    htsFile *out;
    if(opt_.dump_filtered) {
        std::string ofname = opt_.prefix + ".hic.filtered.bam";
        out = hts_open(ofname.c_str(), "wb");
        if(sam_hdr_write(out, reader_.Header()) != 0) {
            std::cerr << "[" << GetCurTime() << "] Error failed to write header to " << ofname << "\n";
            exit(EXIT_FAILURE);
        }
    }
    matrix_.Resize(reader_.Header()->n_targets);
    std::size_t SIZE = 100000;
    std::vector<std::array<bam1_t*, 2>> reads(SIZE);
    for(std::size_t i = 0; i < SIZE; ++i) {
        reads[i][0] = bam_init1();
        reads[i][1] = bam_init1();
    }
    long count = 0, filtered = 0, snp = 0;
    while(1) {
        std::size_t size = 0;
        for(std::size_t i = 0; i < SIZE; ++i) {
            int r = reader_.LoadPair(reads[i]);
            if(r < 0) break;
            size += 1;
        }
        if(size == 0) break;
        // std::cerr << "size " << size << "\n";
        int flag[size];
        std::atomic<std::size_t> index { 0 };
        auto func = [&](std::size_t tid) {
            for(auto cur = index.fetch_add(1); cur < size; cur = index.fetch_add(1)) {
                flag[cur] = FilterHiCPair(reads[cur]);
            }
        };
        MultiThreads(opt_.threads, func);
        for(std::size_t i = 0; i < size; ++i) {
            if(flag[i] != 0) {
                std::string name1 = bam_get_qname(reads[i][0]);
                std::string name2 = bam_get_qname(reads[i][1]);
                if(name1 != name2) {
                    std::cerr << "[" << GetCurTime() << "] Read pair " << name1 << " != " << name2 << "\n";
                    exit(EXIT_FAILURE);
                }
                filtered += 1;
                if(matrix_.AddLink(reads[i][0]->core.tid, reads[i][1]->core.tid) < 1) {
                    std::cerr << "[" << GetCurTime() << "] Error in AddLink\n";
                    exit(EXIT_FAILURE);
                }
                if(matrix_.AddLink(reads[i][1]->core.tid, reads[i][0]->core.tid) < 1) {
                    std::cerr << "[" << GetCurTime() << "] Error in AddLink\n";
                    exit(EXIT_FAILURE);
                }
                if(opt_.dump_filtered) {
                    sam_write1(out, reader_.Header(), reads[i][0]);
                    sam_write1(out, reader_.Header(), reads[i][1]);
                }
                if(flag[i] == 2 && reads[i][0]->core.tid == reads[i][1]->core.tid) {
                    snp += 1;
                    std::string ref_name = reader_.Header()->target_name[reads[i][0]->core.tid];
                    IncreaseCov(ref_name, reads[i][0], reads[i][1]);
                }
            }
        }
        count += size;
        if((count % 1000000) == 0) {
            std::cerr << "[" << GetCurTime() << "] Parsed " << count << " read pairs\n";
        }
    }
    // std::cerr << "[" << GetCurTime() << "] Total " << count << " read pairs\n";
    std::cerr << "[" << GetCurTime() << "] " << filtered << "/" << count << " pass filtered " << snp << "/" << count << " SNP\n";
    for(std::size_t i = 0; i < SIZE; ++i) {
        bam_destroy1(reads[i][0]);
        bam_destroy1(reads[i][1]);
    }
    if(opt_.dump_filtered) hts_close(out);
}

void Phase::GenerateMatrix() {
    reader_.Initialize(opt_.hfname);
    for(int32_t i = 0; i < reader_.Header()->n_targets; ++i) {
        cov_[reader_.Header()->target_name[i]].resize(reader_.Header()->target_len[i] + 1, 0);
    }
    matrix_.Resize(reader_.Header()->n_targets);
    std::array<bam1_t*, 2> pair;
    pair[0] = bam_init1();
    pair[1] = bam_init1();
    long count = 0;
    while(1) {
        int r = reader_.LoadPair(pair);
        if(r < 0) break;
        std::string name1 = bam_get_qname(pair[0]);
        std::string name2 = bam_get_qname(pair[1]);
        if(name1 != name2) {
            std::cerr << "[" << GetCurTime() << "] Read pair " << name1 << " != " << name2 << "\n";
            exit(EXIT_FAILURE);
        }
        if(matrix_.AddLink(pair[0]->core.tid, pair[1]->core.tid) < 1) {
            std::cerr << "[" << GetCurTime() << "] Error in AddLink\n";
            exit(EXIT_FAILURE);
        }
        if(matrix_.AddLink(pair[1]->core.tid, pair[0]->core.tid) < 1) {
            std::cerr << "[" << GetCurTime() << "] Error in AddLink\n";
            exit(EXIT_FAILURE);
        }
        if(pair[0]->core.tid == pair[1]->core.tid) {
            std::string ref_name = reader_.Header()->target_name[pair[0]->core.tid];
            IncreaseCov(ref_name, pair[0], pair[1]);
        }
        count += 1;
        if((count % 1000000) == 0) {
            std::cerr << "[" << GetCurTime() << "] Parsed " << count << " read pairs\n";
        }
    }
    bam_destroy1(pair[0]);
    bam_destroy1(pair[1]);
}

std::vector<std::array<std::string, 2>> Phase::LoadPair() {
    std::ifstream pin(opt_.pfname);
    if(!pin.is_open()) {
        std::cerr << "[" << GetCurTime() << "] Error could not open " << opt_.pfname << " for reading\n";
        exit(EXIT_FAILURE);
    }
    std::vector<std::array<std::string, 2>> pairs;
    std::string line;
    while(std::getline(pin, line)) {
        auto items = SplitString(line, '\t');
        pairs.emplace_back(std::array<std::string, 2> { items[0], items[1] });
    }
    return pairs;
}

void Phase::ClusterCtg(std::vector<std::array<std::string, 2>> &pairs) {
    std::string prev_pri;
    Ctg group;
    for(auto p: pairs) {
        auto items = SplitString(p[0], ':');
        if(prev_pri.empty() || items[0] != prev_pri) {
            if(!prev_pri.empty()) contigs_.push_back(group);
            prev_pri = items[0];
            group.num = 0;
            group.name = items[0];
            group.olps.clear();
        }
        group.num += 1;
        Overlap o;
        o.first = name2id_[p[0]];
        o.second = name2id_[p[1]];
        // o.first = rs_.GetIdByName(p[0]);
        // o.second = rs_.GetIdByName(p[1]);
        int init = rand() & 1;
        o.label = init;
        o.label_tmp = init ^ 1;
        o.count1 = 0;
        o.count2 = 0;
        group.olps.emplace_back(o);
    }
    contigs_.push_back(group);
}

void Phase::ClusterCtgScaffold(std::vector<std::array<std::string, 2>> &pair) {
    std::unordered_map<std::string, std::string> pair_map;
    for(auto p: pair) {
        pair_map[p[0]] = p[1];
        pair_map[p[1]] = p[0];
    }
    GzFileReader reader(opt_.sfname);
    if(reader.Valid()) {
        std::string line = reader.GetNoEmptyLine();
        std::size_t index = 0;
        while(!line.empty()) {
            if(line[0] != '#') {
                Ctg group;
                auto items = SplitString(line, '\t');
                std::string group_name = "cluster" + std::to_string(index);
                group.name = group_name;
                group.num = 0;
                group.olps.clear();
                for(auto item: items) {
                    Overlap o;
                    o.first = name2id_[item];
                    o.second = name2id_[pair_map[item]];
                    int init = rand() & 1;
                    o.label = init;
                    o.label_tmp = init ^ 1;
                    o.count1 = 0;
                    o.count2 = 0;
                    group.olps.emplace_back(o);
                }
                contigs_.emplace_back(group);
                index++;
            }
            line = reader.GetNoEmptyLine();
        }
    } else {
        std::cerr << "[" << GetCurTime() << "] Failed to open " << opt_.sfname << " for reading\n";
        exit(EXIT_FAILURE);
    }
}

int Phase::GetMaxDiag(std::vector<Overlap> &olp) const {
    int max = 0;
    for(std::size_t i = 0; i < olp.size(); ++i) {
        for(std::size_t j = 0; j < olp.size(); ++j) {
            if(j == i) continue;
            int c1 = (int)matrix_.Get(olp[i].first, olp[j].first);
            int c2 = (int)matrix_.Get(olp[i].first, olp[j].second);
            int c3 = (int)matrix_.Get(olp[i].second, olp[j].first);
            int c4 = (int)matrix_.Get(olp[i].second, olp[j].second);
            if(olp[i].label == olp[j].label) {
                if(c1 * c4 > 0 && c2 + c3 == 0) max += 1;
            } else {
                if(c2 * c3 > 0 && c1 + c4 == 0) max += 1;
            }
        }
    }
    return max;
}

int Phase::GetMaxNoDiag(std::vector<Overlap> &olp, std::unordered_set<uint32_t> &flit) const {
    int max = 0;
    for(std::size_t i = 0; i < olp.size(); ++i) {
        for(std::size_t j = 0; j < olp.size(); ++j) {
            if(j == i) continue;
            int c1 = (int)matrix_.Get(olp[i].first, olp[j].first);
            int c2 = (int)matrix_.Get(olp[i].first, olp[j].second);
            int c3 = (int)matrix_.Get(olp[i].second, olp[j].first);
            int c4 = (int)matrix_.Get(olp[i].second, olp[j].second);
            if(olp[i].label == olp[j].label) {
                if(c1 + c4 > c2 + c3) max += 1;
            } else {
                if(c2 + c3 > c1 + c4) max += 1;
            }
        }
    }
    return max;
}

int Phase::GetMaxAll(std::vector<Overlap> &olp) const {
    int max = 0;
    for(std::size_t i = 0; i < olp.size(); ++i) {
        for(std::size_t j = 0; j < olp.size(); ++j) {
            if(j == i) continue;
            int c1 = (int)matrix_.Get(olp[i].first, olp[j].first);
            int c2 = (int)matrix_.Get(olp[i].first, olp[j].second);
            int c3 = (int)matrix_.Get(olp[i].second, olp[j].first);
            int c4 = (int)matrix_.Get(olp[i].second, olp[j].second);
            if(olp[i].label == olp[j].label) {
                if(c1 + c4 > c2 + c3) max += 1;
            } else {
                if(c2 + c3 > c1 + c4) max += 1;
            }
        }
    }
    return max;
}

int Phase::GetMaxAllCount(std::vector<Overlap> &olp) const {
    int max = 0;
    for(std::size_t i = 0; i < olp.size(); ++i) {
        for(std::size_t j = 0; j < olp.size(); ++j) {
            if(j == i) continue;
            int c1 = (int)matrix_.Get(olp[i].first, olp[j].first);
            int c2 = (int)matrix_.Get(olp[i].first, olp[j].second);
            int c3 = (int)matrix_.Get(olp[i].second, olp[j].first);
            int c4 = (int)matrix_.Get(olp[i].second, olp[j].second);
            if(olp[i].label == olp[j].label) {
                max = max + (c1 + c4) - (c2 + c3);
            } else {
                max = max + (c2 + c3) - (c1 + c4);
            }
        }
    }
    return max;
}

void Phase::LocalMaximumDiag(std::vector<Overlap> &olp) {
    for(std::size_t i = 0; i < olp.size(); ++i) {
        int count = 0;
        for(std::size_t j = 0; j < olp.size(); ++j) {
            if(j == i) continue;
            int c1 = (int)matrix_.Get(olp[i].first, olp[j].first);
            int c2 = (int)matrix_.Get(olp[i].first, olp[j].second);
            int c3 = (int)matrix_.Get(olp[i].second, olp[j].first);
            int c4 = (int)matrix_.Get(olp[i].second, olp[j].second);
            if(olp[i].label == olp[j].label) {
                if(c1 * c4 > 0 && c2 + c3 == 0) count -= 1;
                else if(c2 * c3 > 0 && c1 + c4 == 0) count += 1;
            } else {
                if(c1 * c4 > 0 && c2 + c3 == 0) count += 1;
                else if(c2 * c3 > 0 && c1 + c4 == 0) count -= 1;
            }
        }
        if(count > 0) olp[i].label ^= 1;
    }
}

void Phase::LocalMaximumNoDiag(std::vector<Overlap> &olp, std::unordered_set<uint32_t> &filt) {
    for(std::size_t i = 0; i < olp.size(); ++i) {
        if(filt.find(i) != filt.end()) continue;
        int count = 0;
        for(std::size_t j = 0; j < olp.size(); ++j) {
            if(j == i) continue;
            int c1 = (int)matrix_.Get(olp[i].first, olp[j].first);
            int c2 = (int)matrix_.Get(olp[i].first, olp[j].second);
            int c3 = (int)matrix_.Get(olp[i].second, olp[j].first);
            int c4 = (int)matrix_.Get(olp[i].second, olp[j].second);
            if(olp[i].label == olp[j].label) {
                if(c1 + c4 > c2 + c3) count -= 1;
                else if(c2 + c3 > c1 + c4) count += 1;
            } else {
                if(c1 + c4 > c2 + c3) count += 1;
                else if(c2 + c3 > c1 + c4) count -= 1;
            }
        }
        if(count > 0) olp[i].label ^= 1;
    }
}

void Phase::LocalMaximumAll(std::vector<Overlap> &olp) {
    for(std::size_t i = 0; i < olp.size(); ++i) {
        int count = 0;
        for(std::size_t j = 0; j < olp.size(); ++j) {
            if(j == i) continue;
            int c1 = (int)matrix_.Get(olp[i].first, olp[j].first);
            int c2 = (int)matrix_.Get(olp[i].first, olp[j].second);
            int c3 = (int)matrix_.Get(olp[i].second, olp[j].first);
            int c4 = (int)matrix_.Get(olp[i].second, olp[j].second);
            if(olp[i].label == olp[j].label) {
                if(c1 + c4 > c2 + c3) count -= 1;
                else if(c2 + c3 > c1 + c4) count += 1;
            } else {
                if(c1 + c4 > c2 + c3) count += 1;
                else if(c2 + c3 > c1 + c4) count -= 1;
            }
        }
        if(count > 0) olp[i].label ^= 1;
    }
}

void Phase::LocalMaximumAllCount(std::vector<Overlap> &olp) {
    for(std::size_t i = 0; i < olp.size(); ++i) {
        int count = 0;
        for(std::size_t j = 0; j < olp.size(); ++j) {
            if(j == i) continue;
            int c1 = (int)matrix_.Get(olp[i].first, olp[j].first);
            int c2 = (int)matrix_.Get(olp[i].first, olp[j].second);
            int c3 = (int)matrix_.Get(olp[i].second, olp[j].first);
            int c4 = (int)matrix_.Get(olp[i].second, olp[j].second);
            if(olp[i].label == olp[j].label) {
                count = count + (c2 + c3) - (c1 + c4);
            } else {
                count = count + (c1 + c4) - (c2 + c3);
            }
        }
        if(count > 0) olp[i].label ^= 1;
    }
}

void Phase::Print(std::vector<Overlap> &olp) {
    for(std::size_t i = 0; i < olp.size(); ++i) {
        std::cout << names_[olp[i].first] << " " << names_[olp[i].second] << " ";
    }
    std::cout << "\n";
    for(std::size_t i = 0; i < olp.size(); ++i) {
        for(std::size_t j = 0; j < olp.size(); ++j) {
            std::cout << (int)matrix_.Get(olp[i].first, olp[j].first) << " " << (int)matrix_.Get(olp[i].first, olp[j].second) << " ";
        }
        std::cout << "\n";
        for(std::size_t j = 0; j < olp.size(); ++j) {
            std::cout << (int)matrix_.Get(olp[i].second, olp[j].first) << " " << (int)matrix_.Get(olp[i].second, olp[j].second) << " ";
        }
        std::cout << "\n";
    }
}

void Phase::PrintScore(std::vector<Overlap> &olp, std::ofstream &out) {
    // int consistent = 0, inconsistent = 0;
    for(std::size_t i = 0; i < olp.size(); ++i) {
        std::array<int, 2> consistent { 0, 0 };
        std::array<int, 2> inconsistent { 0, 0 };
        int equal = 0;
        for(std::size_t j = 0; j < olp.size(); ++j) {
            if(j == i) continue;
            int c1 = matrix_.Get(olp[i].first, olp[j].first);
            int c2 = matrix_.Get(olp[i].first, olp[j].second);
            int c3 = matrix_.Get(olp[i].second, olp[j].first);
            int c4 = matrix_.Get(olp[i].second, olp[j].second);
            if(olp[i].label == olp[j].label) {
                consistent[0] = consistent[0] + c1 + c4;
                inconsistent[0] = inconsistent[0] + c2 + c3;
                if(c1 + c4 > c2 + c3) consistent[1] += 1;
                else if(c1 + c4 < c2 + c3) inconsistent[1] += 1;
                else if(c1 + c4 == c2 + c3 && c1 + c4 != 0) equal += 1;
            } else {
                consistent[0] = consistent[0] + c2 + c3;
                inconsistent[0] = inconsistent[0] + c1 + c4;
                if(c2 + c3 > c1 + c4) consistent[1] += 1;
                else if(c2 + c3 < c1 + c4) inconsistent[1] += 1;
                else if(c1 + c4 == c2 + c3 && c1 + c4 != 0) equal += 1;
            }
        }
        out << names_[olp[i].first] << " " << names_[olp[i].second] << " " << consistent[0] << " " << inconsistent[0] << " " << consistent[1] << " " << inconsistent[1] << " " << equal << "\n";
    }
}

void Phase::Phasing() {
    std::string ofname = opt_.prefix + ".result.txt";
    std::ofstream out(ofname);
    if(!out.is_open()) {
        std::cerr << "[" << GetCurTime() << "] Error could not open " << ofname << " for writing\n";
        exit(EXIT_FAILURE);
    }
    std::ofstream sout;
    if(opt_.print_score) {
        std::string sfname = opt_.prefix + ".score";
        sout.open(sfname);
        if(!sout.is_open()) {
            std::cerr << "[" << GetCurTime() << "] Error could not open " << sfname << " for writing\n";
            exit(EXIT_FAILURE);
        }
    }

    auto mask_switch_error_func = [&](std::vector<Overlap> &olps) {
        for(std::size_t i = 0; i < olps.size(); ++i) {
            std::string name = names_[olps[i].first];
            if(switch_error_.find(name) == switch_error_.end()) continue;
            for(std::size_t j = 0; j < olps.size(); ++j) {
                matrix_.Set(olps[i].first, olps[j].first, 0);
                matrix_.Set(olps[i].first, olps[j].second, 0);
                matrix_.Set(olps[i].second, olps[j].first, 0);
                matrix_.Set(olps[i].second, olps[j].second, 0);
                matrix_.Set(olps[j].first, olps[i].first, 0);
                matrix_.Set(olps[j].first, olps[i].second, 0);
                matrix_.Set(olps[j].second, olps[i].first, 0);
                matrix_.Set(olps[j].second, olps[i].second, 0);
            }
        }
    };

    std::mutex mtx;
    std::atomic<std::size_t> index { 0 };
    auto func = [&](std::size_t tid) {
        for(auto cur = index.fetch_add(1); cur < contigs_.size(); cur = index.fetch_add(1)) {
            Ctg &c = contigs_[cur];
            if(c.olps.empty()) continue;
            {
                std::lock_guard<std::mutex> lock(mtx);
                std::cerr << "[" << GetCurTime() << "] Thread " << tid << " working on group " << c.name << "\n";
                if(opt_.print_mat) Print(c.olps);
            }
            if(!opt_.scaffold)
                mask_switch_error_func(c.olps);
// /*
            std::array<int, 2> all {0, 0};
            // if(opt_.scaffold) all[0] = GetMaxAllCount(c.olps);
            // else all[0] = GetMaxAll(c.olps);
            all[0] = GetMaxAll(c.olps);
            int best = 0;
            for(int iter = 0; iter < opt_.iteration; ++iter) {
                while(1) {
                    // if(opt_.scaffold) LocalMaximumAllCount(c.olps);
                    // else LocalMaximumAll(c.olps);
                    LocalMaximumAll(c.olps);
                    // if(opt_.scaffold) all[1] = GetMaxAllCount(c.olps);
                    // else all[1] = GetMaxAll(c.olps);
                    all[1] = GetMaxAll(c.olps);
                    if(all[1] > all[0]) all[0] = all[1];
                    else break;
                }
                if(all[0] > best) {
                    for(std::size_t i = 0; i < c.olps.size(); ++i) {
                        c.olps[i].label_tmp = c.olps[i].label;
                    }
                    best = all[0];
                } 
                // else {
                //     for(std::size_t i = 0; i < c.olps.size(); ++i) {
                //         c.olps[i].label = c.olps[i].label_tmp;
                //     }
                // }
                int random = rand() % c.olps.size();
                c.olps[random].label ^= 1;
                for(std::size_t i = 0; i < c.olps.size(); ++i) {
                    if(i == random) continue;
                    int c1 = (int)matrix_.Get(c.olps[random].first, c.olps[i].first);
                    int c2 = (int)matrix_.Get(c.olps[random].first, c.olps[i].second);
                    int c3 = (int)matrix_.Get(c.olps[random].second, c.olps[i].first);
                    int c4 = (int)matrix_.Get(c.olps[random].second, c.olps[i].second);
                    if(c1 + c2 + c3 + c4 > 0) c.olps[i].label ^= 1;
                }
                // if(opt_.scaffold) all[0] = GetMaxAllCount(c.olps);
                // else all[0] = GetMaxAll(c.olps);
                all[0] = GetMaxAll(c.olps);
                all[1] = 0;
            }
            // if(all[0] < best) {
                for(std::size_t i = 0; i < c.olps.size(); ++i) {
                    c.olps[i].label = c.olps[i].label_tmp;
                }
            // }
// */
/*
            std::unordered_set<uint32_t> stable;
            std::array<int, 2> max_diag { 0, 0 }; // last, cur
            for(std::size_t i = 0; i < c.olps.size(); ++i) {
                int count = 0;
                for(std::size_t j = 0; j < c.olps.size(); ++j) {
                    if(j == i) continue;
                    int c1 = (int)matrix_.Get(c.olps[i].first, c.olps[j].first);
                    int c2 = (int)matrix_.Get(c.olps[i].first, c.olps[j].second);
                    int c3 = (int)matrix_.Get(c.olps[i].second, c.olps[j].first);
                    int c4 = (int)matrix_.Get(c.olps[i].second, c.olps[j].second);
                    if((c1 * c4 > 0 && c2 + c3 == 0) || (c2 * c3 > 0 && c1 + c4 == 0)) {
                        stable.insert(i);
                        count += 1;
                    }
                }
                c.olps[i].count1 = count;
            }
            max_diag[0] = GetMaxDiag(c.olps);
            int best_diag = 0;
            for(int iter = 0; iter < opt_.iteration; ++iter) {
                while(1) {
                    LocalMaximumDiag(c.olps);
                    max_diag[1] = GetMaxDiag(c.olps);
                    if(max_diag[1] > max_diag[0]) max_diag[0] = max_diag[1];
                    else break;
                }
                if(max_diag[0] > best_diag) {
                    for(std::size_t i = 0; i < c.olps.size(); ++i) {
                        c.olps[i].label_tmp = c.olps[i].label;
                    }
                    best_diag = max_diag[0];
                } 
                // else {
                //     for(std::size_t i = 0; i < c.olps.size(); ++i) {
                //         c.olps[i].label = c.olps[i].label_tmp;
                //     }
                // }
                int random = rand() % c.olps.size();
                c.olps[random].label ^= 1;
                for(std::size_t i = 0; i < c.olps.size(); ++i) {
                    int c1 = (int)matrix_.Get(c.olps[random].first, c.olps[i].first);
                    int c2 = (int)matrix_.Get(c.olps[random].first, c.olps[i].second);
                    int c3 = (int)matrix_.Get(c.olps[random].second, c.olps[i].first);
                    int c4 = (int)matrix_.Get(c.olps[random].second, c.olps[i].second);
                    if((c1 * c4 > 0 && c2 + c3 == 0) || (c2 * c3 > 0 && c1 + c4 == 0)) c.olps[i].label ^= 1;
                }
                max_diag[0] = GetMaxDiag(c.olps);
                max_diag[1] = 0;
            }
            for(std::size_t i = 0; i < c.olps.size(); ++i) {
                c.olps[i].label = c.olps[i].label_tmp;
            }

            if(stable.size() != c.olps.size()) {
                std::array<int, 2> max_nodiag { 0, 0 };
                for(std::size_t i = 0; i < c.olps.size(); ++i) {
                    if(stable.find(i) != stable.end()) continue;
                    int count = 0;
                    for(std::size_t j = 0; j < c.olps.size(); ++j) {
                        if(j == i) continue;
                        int c1 = (int)matrix_.Get(c.olps[i].first, c.olps[j].first);
                        int c2 = (int)matrix_.Get(c.olps[i].first, c.olps[j].second);
                        int c3 = (int)matrix_.Get(c.olps[i].second, c.olps[j].first);
                        int c4 = (int)matrix_.Get(c.olps[i].second, c.olps[j].second);
                        if((c1 + c4 > c2 + c3) || (c2 + c3 > c1 + c4)) count += 1;
                    }
                    c.olps[i].count2 = count;
                }
                max_nodiag[0] = GetMaxNoDiag(c.olps, stable);
                int best_nodiag = 0;
                for(int iter = 0; iter < opt_.iteration; ++iter) {
                    while(1) {
                        LocalMaximumNoDiag(c.olps, stable);
                        max_nodiag[1] = GetMaxNoDiag(c.olps, stable);
                        if(max_nodiag[1] > max_nodiag[0]) max_nodiag[0] = max_nodiag[1];
                        else break;
                    }
                    if(max_nodiag[0] > best_nodiag) {
                        for(std::size_t i = 0; i < c.olps.size(); ++i) {
                            c.olps[i].label_tmp = c.olps[i].label;
                        }
                        best_nodiag = max_nodiag[0];
                    }
                    int random = rand() % c.olps.size();
                    while(stable.find(random) != stable.end()) random = rand() % c.olps.size();
                    c.olps[random].label ^= 1;
                    for(std::size_t i = 0; i < c.olps.size(); ++i) {
                        if(stable.find(i) != stable.end()) continue;
                        int c1 = (int)matrix_.Get(c.olps[random].first, c.olps[i].first);
                        int c2 = (int)matrix_.Get(c.olps[random].first, c.olps[i].second);
                        int c3 = (int)matrix_.Get(c.olps[random].second, c.olps[i].first);
                        int c4 = (int)matrix_.Get(c.olps[random].second, c.olps[i].second);
                        if((c1 + c4 > c2 + c3) || (c2 + c3 > c1 + c4)) c.olps[i].label ^= 1;
                    }
                    max_nodiag[0] = GetMaxNoDiag(c.olps, stable);
                    max_nodiag[1] = 0;
                }
                for(std::size_t i = 0; i < c.olps.size(); ++i) {
                    c.olps[i].label = c.olps[i].label_tmp;
                }
            }
*/           
            {
                std::lock_guard<std::mutex> lock(mtx);
                if(opt_.print_score) PrintScore(c.olps, sout);
                for(std::size_t i = 0; i < c.olps.size(); ++i) {
                    std::string s1 = names_[c.olps[i].first];
                    std::string s2 = names_[c.olps[i].second];
                    // std::string s1 = rs_.GetNameById(c.olps[i].first);
                    // std::string s2 = rs_.GetNameById(c.olps[i].second);
                    if(c.olps[i].label == 0) {
                        out << c.name << " " << s1 << " " << s2 << "\n";
                        // out << c.name << " " << s1 << " " << s2 << " 0.5 " << c.olps[i].count1 << " " << c.olps[i].count2 << " " << c.olps[i].first << " " << c.olps[i].second << "\n";
                    } else {
                        out << c.name << " " << s2 << " " << s1 << "\n";
                        // out << c.name << " " << s2 << " " << s1 << " 0.5 " << c.olps[i].count1 << " " << c.olps[i].count2 << " " << c.olps[i].second << " " << c.olps[i].first << "\n";
                    }
                }
            }
        }
    };
    MultiThreads(1, func);
    // MultiThreads(opt_.threads, func);
    out.close();
    if(opt_.print_score) sout.close();
}

int Options::Check() {
    // if(bfname.empty() && vfname.empty()) {
    //     std::cerr << "[" << GetCurTime() << "] Please specify a paf file name of alt contig against pri contig [-b | --bam] or a variant file name [-v | --var]\n";
    //     return 1;
    // }
    // if(cfname.empty() || cfname.substr(cfname.rfind('.')) != ".fai") {
    //     std::cerr << "[" << GetCurTime() << "] Please specify a contigs index file name (.fai) [-c | --ctg]\n";
    //     return 1;
    // }
    if(hfname.empty() && mfname.empty()) {
        std::cerr << "[" << GetCurTime() << "] Please specify a hic mapped file name [-i | --hic] or a matrix file name [-m | --mat]\n";
        return 1;
    }
    if(pfname.empty()) {
        std::cerr << "[" << GetCurTime() << "] Please specify a pair file name [-a | pair]\n";
        return 1;
    }
    if(threads < 1) {
        std::cerr << "[" << GetCurTime() << "] threads must be >= 1 [-t | --threads]\n";
        return 1;
    }
    if(iteration < 1) {
        std::cerr << "[" << GetCurTime() << "] iteration must be >= 1 [-n | --iter]\n";
        return 1;
    }
    return 0;
}

void USAGE(Options &opt) {
    std::cerr << "Usage: phase [options] -c ctg.fasta -a pair.txt\n";
    std::cerr << "\tOptions:\n";
    std::cerr << "\t\t-a | --pair       <FILE>  pair file name\n";
    std::cerr << "\t\t-b | --paf        <FILE>  mapped file name of alt against primary\n";
    std::cerr << "\t\t-c | --ctg        <FILE>  minced contig index file name (.fai)\n";
    std::cerr << "\t\t-i | --hic        <FILE>  mapped file of Hi-C read pair against minced contig\n";
    std::cerr << "\t\t-m | --mat        <FILE>  matrix file name of Hi-C mapping\n";
    std::cerr << "\t\t-v | --var        <FILE>  variant file name\n";
    std::cerr << "\t\t-p | --prefix     <STR>   prefix of output files [" << opt.prefix << "]\n";
    std::cerr << "\t\t-t | --threads    <INT>   number of threads [" << opt.threads << "]\n";
    std::cerr << "\t\t-n | --iter       <INT>   number of iteration of phasing [" << opt.iteration << "]\n";
    std::cerr << "\t\t-s | --seed       <INT>   seed for random function [" << opt.seed << "]\n";
    std::cerr << "\t\t-q | --mapq       <INT>   minimum mapping quality [" << opt.mapq << "]\n";
    std::cerr << "\t\t     --priclair   <FILE>  variants in primary contigs\n";
    std::cerr << "\t\t     --altclair   <FILE>  variants in alt contigs\n";
    std::cerr << "\t\t     --cluster    <FILE>  clusters.by_name.txt generated by lachesis\n";
    std::cerr << "\t\t     --bed1       <FILE>  {prefix}.phased0.bed\n";
    std::cerr << "\t\t     --bed2       <FILE>  {prefix}.phased1.bed\n";
    std::cerr << "\t\t     --dump-var           dump different site between primary contigs and alt contigs [" << opt.dump_var << "]\n";
    std::cerr << "\t\t     --dump-filt          dump filtered Hi-C mapping [" << opt.dump_filtered << "]\n";
    std::cerr << "\t\t     --dump-mat           dump matrix of Hi-C mapping [" << opt.dump_mat << "]\n";
    std::cerr << "\t\t     --filtered           Hi-C mapping are filtered [" << opt.filtered << "]\n";
    std::cerr << "\t\t     --print-mat          print matrix of Hi-C mapping [" << opt.print_mat << "]\n";
    std::cerr << "\t\t     --scaffold           phasing the result of scaffolding\n";
    std::cerr << "\t\t     --porec              phasing the result of pore-c\n";
    std::cerr << "\t\t-h | --help               display this message\n";
}

int ParseArgument(int argc, char **argv, Options &opt) {
    struct option long_opt[] = {
        {"pair",        required_argument,  NULL,   'a'}, 
        {"paf",         required_argument,  NULL,   'b'}, 
        {"ctg",         required_argument,  NULL,   'c'}, 
        {"var",         required_argument,  NULL,   'v'}, 
        {"prefix",      required_argument,  NULL,   'p'}, 
        {"hic",         optional_argument,  NULL,   'i'}, 
        {"mat",         optional_argument,  NULL,   'm'}, 
        {"threads",     optional_argument,  NULL,   't'}, 
        {"iter",        optional_argument,  NULL,   'n'}, 
        {"seed",        optional_argument,  NULL,   's'}, 
        {"mapq",        optional_argument,  NULL,   'q'}, 
        {"help",        no_argument,        NULL,   'h'}, 
        {"dump-var",    no_argument,        NULL,   300}, 
        {"dump-filt",   no_argument,        NULL,   301}, 
        {"dump-mat",    no_argument,        NULL,   302}, 
        {"filtered",    no_argument,        NULL,   303}, 
        {"print-mat",   no_argument,        NULL,   304}, 
        {"print-score", no_argument,        NULL,   305}, 
        {"priclair",    required_argument,  NULL,   306}, 
        {"altclair",    required_argument,  NULL,   307}, 
        {"cluster",     required_argument,  NULL,   308}, 
        {"bed1",        required_argument,  NULL,   309}, 
        {"bed2",        required_argument,  NULL,   310}, 
        {"scaffold",    no_argument,        NULL,   311}, 
        {"porec",       no_argument,        NULL,   312}, 
        {0, 0,  0,  0}
    };
    const char *short_opt = "a:b:c:i:v:p:t:n:s:m:e:f:q:h";
    int c;
    while((c = getopt_long(argc, argv, short_opt, long_opt, NULL)) != -1) {
        if(c == 'a') opt.pfname = optarg;
        else if(c == 'b') opt.bfname = optarg;
        else if(c == 'c') opt.cfname = optarg;
        else if(c == 'i') opt.hfname = optarg;
        else if(c == 'v') opt.vfname = optarg;
        else if(c == 'p') opt.prefix = optarg;
        else if(c == 'm') opt.mfname = optarg;
        else if(c == 't') opt.threads = atoi(optarg);
        else if(c == 'n') opt.iteration = atoi(optarg);
        else if(c == 's') opt.seed = atoi(optarg);
        else if(c == 'q') opt.mapq = atoi(optarg);
        else if(c == 300) opt.dump_var = true;
        else if(c == 301) opt.dump_filtered = true;
        else if(c == 302) opt.dump_mat = true;
        else if(c == 303) opt.filtered = true;
        else if(c == 304) opt.print_mat = true;
        else if(c == 305) opt.print_score = true;
        else if(c == 306) opt.pcfname = optarg;
        else if(c == 307) opt.acfname = optarg;
        else if(c == 308) opt.sfname = optarg;
        else if(c == 309) opt.b1fname = optarg;
        else if(c == 310) opt.b2fname = optarg;
        else if(c == 311) opt.scaffold = true;
        else if(c == 312) opt.porec = true;
        else if(c == 'h') {
            USAGE(opt);
            exit(EXIT_SUCCESS);
        } else if(c == ':') {
            std::cerr << "[" << GetCurTime() << "] ERROR missing argument in option " << optopt << "\n";
            return 1;
        } else if(c == '?') {
            std::cerr << "[" << GetCurTime() << "] ERROR unknown option " << optopt << "\n";
            return 1;
        }
    }
    return opt.Check();
}

void Phase::Run() {
    std::ios::sync_with_stdio(false);
    srand(opt_.seed);
    std::ifstream in(opt_.cfname);
    std::string line;
    // load names mapping
    while(std::getline(in, line)) {
        auto items = SplitString(line, '\t');
        name2id_[items[0]] = names_.size();
        names_.emplace_back(items[0]);
    }
    // load variants if specified, otherwise generate variants
    // if clair3 variants are specified, generate variants with clair3, 
    // otherwise generate variants with mapped alt contgis against primary contigs
    if(!opt_.vfname.empty()) {
        LoadVariants();
    } else {
        if(!opt_.pcfname.empty() && !opt_.acfname.empty()) {
            GenerateVariantsWithClair3();
        } else {
            GenerateVariants();
        }
    }
    // load matrix if specified, otherwise generate matrix
    // if bed files are specified, generate matrix with bed files
    // if the mapping is filtered, generate matrix with filtered mapping
    // otherwise generate matrix with all Hi-C or Pore-C mapping
    if(!opt_.mfname.empty()) {
        matrix_.Load(opt_.mfname);
    } else if(!opt_.b1fname.empty() && !opt_.b2fname.empty()) {
        auto bed = LoadBed();
        GenerateMatrixScaffold(bed);
    } else if(!opt_.hfname.empty()) {
        if(opt_.filtered) {
            GenerateMatrix();
        } else {
            if(opt_.porec) {
                FilterAndGenerateMatrixPoreC();
            } else {
                FilterAndGenerateMatrix();
            }
        }
    } else {
        std::cerr << "[" << GetCurTime() << "] Error please specify a matrix file name [-m | --mat] or a hic file name [-i | --hic]\n";
        exit(EXIT_FAILURE);
    }
    // dump matrix if specified
    if(opt_.dump_mat) {
        std::string mfname = opt_.prefix + ".binmat";
        matrix_.Dump(mfname);
    }
    // load pair file
    auto pair = LoadPair();
    // detect switch error if specified
    std::unordered_map<std::string, std::size_t> switch_error;
    if(!opt_.scaffold) {
        switch_error = DetectSwitchError();
        FilterSwitchError(switch_error, pair);
    }
    // cluster contigs
    if(opt_.sfname.empty()) {
        ClusterCtg(pair);
    } else {
        ClusterCtgScaffold(pair);
    }
    // phasing
    Phasing();
    // fix phase if specified
    if(!opt_.scaffold) {
        if(opt_.filtered) {
            FixPhase(opt_.hfname);
        } else {
            std::string fname = opt_.prefix + ".hic.filtered.bam";
            FixPhase(fname);
        }
    }
    // dump switch error if specified
    std::string sfname = opt_.prefix + ".switch";
    std::ofstream out(sfname);
    if(!out.is_open()) {
        std::cerr << "[" << GetCurTime() << "] Failed to open " << sfname << " for writing\n";
        exit(EXIT_FAILURE);
    }
    for(auto ctg: switch_error_) {
        out << ctg.first << "\t" << ctg.second << "\t" << consistency_[ctg.first] << "\n";
    }
    out.close();
}

void Phase::Run1() {
    std::ios::sync_with_stdio(false);
    srand(opt_.seed);
    // rs_.Load(opt_.cfname);
    std::ifstream in(opt_.cfname);
    std::string line;
    while(std::getline(in, line)) {
        auto items = SplitString(line, '\t');
        name2id_[items[0]] = names_.size();
        names_.emplace_back(items[0]);
    }
    if(!opt_.vfname.empty()) {
        LoadVariants();
    } else {
        // GenerateVariants();
        GenerateVariantsWithClair3();
    }
    if(!opt_.b1fname.empty() && !opt_.b2fname.empty()) {
        auto bed = LoadBed();
        GenerateMatrixScaffold(bed);
    } else if(!opt_.mfname.empty()) {
        matrix_.Load(opt_.mfname);
    } else {
        if(opt_.filtered) {
            GenerateMatrix();
        } else {
            FilterAndGenerateMatrix();
        }
    }
    std::unordered_map<std::string, std::size_t> switch_error;
    if(!opt_.scaffold) switch_error = DetectSwitchError();
    if(opt_.dump_mat) {
        std::string mfname = opt_.prefix + ".binmat";
        matrix_.Dump(mfname);
    }
    auto pair = LoadPair();
    if(!opt_.scaffold) FilterSwitchError(switch_error, pair);
    if(!opt_.sfname.empty()) {
        ClusterCtgScaffold(pair);
    } else {
        ClusterCtg(pair);
    }
    Phasing();
    if(opt_.filtered) {
        if(!opt_.scaffold) FixPhase(opt_.hfname);
    } else {
        std::string fname = opt_.prefix + ".hic.filtered.bam";
        if(!opt_.scaffold) FixPhase(fname);
    }
    std::string sfname = opt_.prefix + ".switch";
    std::ofstream out(sfname);
    if(!out.is_open()) {
        std::cerr << "[" << GetCurTime() << "] Failed to open " << sfname << " for writing\n";
        exit(EXIT_FAILURE);
    }
    for(auto ctg: switch_error_) {
        out << ctg.first << "\t" << ctg.second << "\t" << consistency_[ctg.first] << "\n";
    }
    out.close();
}

int main(int argc, char **argv) {
    Options opt;
    if(ParseArgument(argc, argv, opt) == 1) {
        USAGE(opt);
        exit(EXIT_FAILURE);
    }
    Phase p(opt);
    p.Run();
    return 0;
}