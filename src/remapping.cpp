#include "remapping.hpp"
#include "utility.hpp"

#include "minimap2-2.22/minimap.h"
#include "minimap2-2.22/kseq.h"

#include <zlib.h>

#include <getopt.h>

#include <regex>
#include <mutex>
#include <atomic>

KSEQ_INIT(gzFile, gzread)

void USAGE(Options &opt) {
    std::cerr << "Usage:\t";
    std::cerr << "remap [options] -A alt.fasta -P pri.fasta -a aln.paf\n\n";
    std::cerr << "\tOptions:\n";
    std::cerr << "\t\t-P | --pri        <FILE> primary contig file name\n";
    std::cerr << "\t\t-A | --alt        <FILE> alt contig file name\n";
    std::cerr << "\t\t-a | --paf        <FILE> mapping file of alt contigs against primary in paf format with cigar\n";
    std::cerr << "\t\t-B | --breakpoint <FILE> breakpoing file name\n";
    std::cerr << "\t\t-p | --prefix     <STR>  prefix of output file [" << opt.prefix << "]\n";
    std::cerr << "\t\t-b | --block      <INT>  block length [" << opt.block << "]\n";
    std::cerr << "\t\t-t | --threads    <INT>  number of threads [" << opt.threads << "]\n";
    std::cerr << "\t\t-h | --help              display this message\n";
    std::cerr << "\t\t     --break             break contigs [" << opt.breaks << "]\n";
}

int ParseArgument(int argc, char **argv, Options &opt) {
    struct option long_opt[] = {
        {"paf",         required_argument,  NULL,   'a'}, 
        {"pri",         optional_argument,  NULL,   'P'}, 
        {"alt",         optional_argument,  NULL,   'A'}, 
        {"breakpoint",  optional_argument,  NULL,   'B'}, 
        {"prefix",      optional_argument,  NULL,   'p'}, 
        {"block",       optional_argument,  NULL,   'b'}, 
        {"threads",     optional_argument,  NULL,   't'}, 
        {"help",        no_argument,        NULL,   'h'}, 
        {"break",       optional_argument,  NULL,   300}, 
        {0, 0,  0,  0}
    };
    const char *short_opt = "A:a:B:b:P:p:t:h";
    int c;
    while((c = getopt_long(argc, argv, short_opt, long_opt, NULL)) != -1) {
        if(c == 'P') opt.prifname = optarg;
        else if(c == 'A') opt.altfname = optarg;
        else if(c == 'a') opt.paffname = optarg;
        else if(c == 'B') opt.breakfname = optarg;
        else if(c == 'p') opt.prefix = optarg;
        else if(c == 'b') opt.block = atoi(optarg);
        else if(c == 't') opt.threads = atoi(optarg);
        else if(c == 300) opt.breaks = false;
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
    return opt.Check();
}

int Options::Check() {
    if(paffname.empty()) {
        std::cerr << "[" << GetCurTime() << "] Please specify a mapped file [-a | --paf]\n";
        return 1;
    }
    if(prifname.empty() && altfname.empty() && breakfname.empty()) {
        std::cerr << "[" << GetCurTime() << "] Please specify a primary contig file name [-P | --pri] and an alt contig file name [-A | --alt] or a breakpoint file name [-B | --breakpoing]\n";
        return 1;
    }
    if(!prifname.empty() && altfname.empty() && breakfname.empty()) {
        std::cerr << "[" << GetCurTime() << "] Please specify an alt contig file name [-A | --alt]\n";
        return 1;
    }
    if(!altfname.empty() && prifname.empty() && breakfname.empty()) {
        std::cerr << "[" << GetCurTime() << "] Please specify a primary contif file name [-P | --pri]\n";
        return 1;
    }
    if(threads < 1) {
        std::cerr << "[" << GetCurTime() << "] threads must be > 0 [-t | --threads]\n";
        return 1;
    }
    if(block < 0) {
        std::cerr << "[" << GetCurTime() << "] block must be > 0 [-b | --block]\n";
        return 1;
    }
    return 0;
}

void Remapping::LoadPaf() {
    std::ifstream in(opt_.paffname, std::ios::in);
    if(!in.is_open()) {
        std::cerr << "[" << GetCurTime() << "] Could not open " << opt_.paffname << " for reading\n";
        exit(EXIT_FAILURE);
    }
    std::string line;
    while(std::getline(in, line)) {
        if(line.empty()) continue;
        auto items = SplitString(line, '\t');
        std::array<std::string, 13> tmp;
        for(std::size_t i = 0; i < items.size(); ++i) {
            if(i < 12) tmp[i] = items[i];
            else if(items[i].substr(0, 5) == "cg:Z:") {
                tmp[12] = items[i];
                break;
            }
        }
        paf_[items[0]].emplace_back(tmp);
    }
    std::cerr << "[" << GetCurTime() << "] Load " << paf_.size() << " alignments from file " << opt_.paffname << "\n";
    in.close();
}

void Remapping::DumpPaf() {
    std::vector<std::string> keys;
    keys.reserve(paf_.size());
    for(auto k: paf_) keys.emplace_back(k.first);
    std::sort(keys.begin(), keys.end());
    std::string ofname = opt_.prefix + "." + opt_.paffname + ".new.paf";
    std::cerr << "[" << GetCurTime() << "] Dump paf record into " << ofname << "\n";
    std::ofstream out(ofname);
    if(!out.is_open()) {
        std::cerr << "[" << GetCurTime() << "] Could not open " << ofname << " for writing\n";
        exit(EXIT_FAILURE);
    }
    for(auto ctg: keys) {
        for(auto paf: paf_[ctg]) {
            for(std::size_t i = 0; i < paf.size(); ++i) {
                out << paf[i] << "\t";
            }
            out << "\n";
        }
    }
    out.close();
}

void Remapping::Remap() {
    // contigs need remapping
    std::vector<std::string> ctgs;
    int fuzzy = 100;
    for(auto &ctg: paf_) {
        ctgs.emplace_back(ctg.first);
    }

    // remapping
    std::mutex mtx;
    std::atomic<std::size_t> total_l { 0 };
    std::atomic<std::size_t> total { 0 };
    std::atomic<std::size_t> index { 0 };

    auto func_sort = [&](std::size_t tid) {
        for(auto cur = index.fetch_add(1); cur < ctgs.size(); cur = index.fetch_add(1)) {
            std::string ctg = ctgs[cur];
            std::sort(paf_[ctg].begin(), paf_[ctg].end(), [](std::array<std::string, 13> &a, std::array<std::string, 13> &b) {
                return (atoi(a[3].c_str()) - atoi(a[2].c_str()) > atoi(b[3].c_str()) - atoi(b[2].c_str()));
            });
            paf_[ctg].resize(1);
        }
    };
    MultiThreads(opt_.threads, func_sort);

    index.fetch_and(0);
    auto func = [&](std::size_t tid) {
        for(auto cur = index.fetch_add(1); cur < ctgs.size(); cur = index.fetch_add(1)) {
            std::array<std::string, 13> item;
            assert(paf_[ctgs[cur]].size() == 1);
            if(atoi(paf_[ctgs[cur]][0][2].c_str()) <= fuzzy && atoi(paf_[ctgs[cur]][0][1].c_str()) - atoi(paf_[ctgs[cur]][0][3].c_str()) <= fuzzy) continue;
            else item = paf_[ctgs[cur]][0];

            // set minimap2's parameters (preset asm20)
            mm_idxopt_t idxopt;
            mm_mapopt_t mapopt;
            mm_set_opt(0, &idxopt, &mapopt);
            mm_set_opt("asm20", &idxopt, &mapopt);
            mapopt.flag |= MM_F_CIGAR;
            mapopt.mid_occ = 0;
            mapopt.mid_occ_frac = 0.00001;
            mapopt.min_mid_occ = 0;
            mapopt.max_mid_occ = 2000;
            mapopt.end_bonus = 100;
            mapopt.flag |= MM_F_NO_PRINT_2ND;
            mm_verbose = 2;

            int shift_l = atoi(item[2].c_str());
            int beg = atoi(item[7].c_str()) - shift_l * 2;
            if(beg < 0) beg = 0;
            int size_ref = atoi(item[1].c_str());
            std::string ref = *pstore_.GetSeq(item[5]).ToString();
            std::string ctg = *astore_.GetSeq(item[0]).ToString();

            std::string ref_tmp_fname = item[5] + std::to_string(tid) + ".ref.tmp.fasta";
            std::string ctg_tmp_fname = item[5] + std::to_string(tid) + ".ctg.tmp.fasta";

            std::ofstream ref_tmp(ref_tmp_fname);
            if(!ref_tmp.is_open()) {
                std::lock_guard<std::mutex> lock(mtx);
                std::cerr << "[" << GetCurTime() << "] Could not open " << ref_tmp_fname << " for writing\n";
                return;
            }
            ref_tmp << ">" << item[5] << "\n";
            ref_tmp << ref;
            ref_tmp.close();

            std::ofstream ctg_tmp(ctg_tmp_fname);
            if(!ctg_tmp.is_open()) {
                std::lock_guard<std::mutex> lock(mtx);
                std::cerr << "[" << GetCurTime() << "] Could not open " << ctg_tmp_fname << " for writing\n";
                return;
            }
            ctg_tmp << ">" << item[0] << "\n";
            ctg_tmp << ctg;
            ctg_tmp.close();

            gzFile cf = gzopen(ctg_tmp_fname.c_str(), "r");
            assert(cf);
            kseq_t *ks = kseq_init(cf);
            mm_idx_reader_t *r = mm_idx_reader_open(ref_tmp_fname.c_str(), &idxopt, 0);
            mm_idx_t *mi;
            while((mi = mm_idx_reader_read(r, 1)) != 0) {
                mm_mapopt_update(&mapopt, mi);
                mm_tbuf_t *tbuf = mm_tbuf_init();
                gzrewind(cf);
                while(kseq_read(ks) >= 0) {
                    int n_reg = 0;
                    mm_reg1_t *reg;
                    reg = mm_map(mi, ks->seq.l, ks->seq.s, &n_reg, tbuf, &mapopt, 0);
                    // if(n_reg >= 1) paf_[ctgs[cur]].clear();
                    if(n_reg > 1) {
                        total.fetch_add(1);
                        total_l.fetch_add(ks->seq.l);
                    } else if(n_reg == 1) {
                        if(reg[0].qs > fuzzy && ks->seq.l - reg[0].qe > fuzzy) {
                            total.fetch_add(1);
                            total_l.fetch_add(ks->seq.l);
                        }
                    }
                    std::sort(reg, reg + n_reg, [](mm_reg1_t &a, mm_reg1_t &b) {
                        return (a.qe - a.qs > b.qe - b.qs);
                    });
                    assert(reg[0].p);
                    assert(!reg[0].sam_pri || reg[0].id == reg[0].parent);
                    if((mapopt.flag & MM_F_NO_PRINT_2ND) && reg[0].id != reg[0].parent) continue;
                    paf_[ctgs[cur]].clear();
                    item[2] = std::to_string(reg[0].qs);
                    item[3] = std::to_string(reg[0].qe);
                    item[4] = "+-"[reg[0].rev];
                    item[7] = std::to_string(reg[0].rs);
                    item[8] = std::to_string(reg[0].re);
                    item[9] = std::to_string(reg[0].mlen);
                    item[10] = std::to_string(reg[0].blen);
                    item[11] = std::to_string(reg[0].mapq);
                    std::ostringstream oss;
                    oss << "cg:Z:";
                    for(int j = 0; j < reg[0].p->n_cigar; ++j) {
                        oss << (reg[0].p->cigar[j] >> 4) << (MM_CIGAR_STR[reg[0].p->cigar[j] & 0xf]);
                    }
                    item[12] = oss.str();
                    paf_[ctgs[cur]].push_back(item);
                    // for(int i = 0; i < n_reg; ++i) {
                    //     mm_reg1_t *r = &reg[i];
                    //     assert(r->p);
                    //     assert(!r->sam_pri || r->id == r->parent);
                    //     if((mapopt.flag & MM_F_NO_PRINT_2ND) && r->id != r->parent) continue;
                    //     item[2] = std::to_string(r->qs);
                    //     item[3] = std::to_string(r->qe);
                    //     item[4] = "+-"[r->rev];
                    //     item[7] = std::to_string(r->rs);
                    //     item[8] = std::to_string(r->re);
                    //     item[9] = std::to_string(r->mlen);
                    //     item[10] = std::to_string(r->blen);
                    //     item[11] = std::to_string(r->mapq);

                    //     std::ostringstream oss;
                    //     oss << "cg:Z:";
                    //     for(int j = 0; j < r->p->n_cigar; ++j) {
                    //         oss << (r->p->cigar[j] >> 4) << (MM_CIGAR_STR[r->p->cigar[j] & 0xf]);
                    //     }
                    //     item[12] = oss.str();
                    //     paf_[ctgs[cur]].emplace_back(item);
                    // }
                    free(reg);
                }
                mm_tbuf_destroy(tbuf);
                mm_idx_destroy(mi);
            }
            mm_idx_reader_close(r);
            kseq_destroy(ks);
            gzclose(cf);
            if(access(ref_tmp_fname.c_str(), F_OK) == 0) remove(ref_tmp_fname.c_str());
            if(access(ctg_tmp_fname.c_str(), F_OK) == 0) remove(ctg_tmp_fname.c_str());
        }
    };
    MultiThreads(opt_.threads, func);
    std::cerr << "[" << GetCurTime() << "] After remapping " << total.load() << " contig(s) still need remapping. Total length: " << total_l.load() << "\n";
}

std::unordered_map<std::string, std::vector<std::size_t>> Remapping::LoadBps() {
    std::unordered_map<std::string, std::vector<std::size_t>> bps;
    if(!opt_.breakfname.empty()) {
        std::ifstream in(opt_.breakfname);
        if(!in.is_open()) {
            std::cerr << "[" << GetCurTime() << "] Could not open " << opt_.breakfname << " for reading\n";
            exit(EXIT_FAILURE);
        }
        std::string line;
        while(std::getline(in, line)) {
            auto items = SplitString(line, '\t');
            bps[items[0]].emplace_back(atoi(items[1].c_str()));
        }
        if(!opt_.breaks) {
            for(auto &ctg: bps) {
                if(ctg.second.size() > 2) {
                    ctg.second[1] = ctg.second.back();
                    ctg.second.resize(2);
                }
            }
        }
    } else {
        ReadStore rs;
        rs.Load(opt_.altfname);
        for(std::size_t i = 0; i < rs.Size(); ++i) {
            bps[rs.GetNameById(i)].emplace_back(0);
            bps[rs.GetNameById(i)].emplace_back(rs.GetSeqLength(i));
        }
    }
    return bps;
}

void Remapping::Cut() {
    // std::unordered_map<std::string, std::vector<std::size_t>> bps;
    // std::ifstream in(opt_.breakfname);
    // if(!in.is_open()) {
    //     std::cerr << "[" << GetCurTime() << "] Could not open " << opt_.breakfname << " for reading\n";
    //     exit(EXIT_FAILURE);
    // }
    // std::string line;
    // while(std::getline(in, line)) {
    //     auto items = SplitString(line, '\t');
    //     bps[items[0]].emplace_back(atoi(items[1].c_str()));
    // }
    auto bps = LoadBps();
    // ctg cl cs ce ref rl rs re strand
    std::vector<std::vector<std::array<std::string, 9>>> cutsites(opt_.threads);
    std::regex pattern("(\\d+)([MIDX=]?)");
    std::vector<std::string> ctgs;
    for(auto ctg: paf_) ctgs.emplace_back(ctg.first);
    std::sort(ctgs.begin(), ctgs.end());

    std::mutex mtx;
    std::atomic<std::size_t> index { 0 };
    int fuzzy = 1000;
    int num_skip = 0, num_clip = 0;
    auto func = [&](std::size_t tid) {
        for(auto cur = index.fetch_add(1); cur < ctgs.size(); cur = index.fetch_add(1)) {
            std::string ctg = ctgs[cur];
            if(paf_[ctg].size() != 1) {
                std::lock_guard<std::mutex> lock(mtx);
                std::cerr << ctg << " " << paf_[ctg].size() << "\n";
            }
            assert(paf_[ctg].size() == 1);
            if(paf_[ctg][0][11] != "60") continue;
            // std::sort(paf_[ctg].begin(), paf_[ctg].end(), 
            //     [](std::array<std::string, 13> &a, std::array<std::string, 13> &b) {
            //         return atoi(a[2].c_str()) < atoi(b[2].c_str());
            //     });
            // int prev_ctg = 0, prev_ref = atoi(paf_[ctg][0][7].c_str());
            std::string ref = paf_[ctg][0][5];
            std::string strand = paf_[ctg][0][4];
            // bool skip = false;
            // for(auto items: paf_[ctg]) {
            //     int s_ctg = atoi(items[2].c_str());
            //     int s_ref = atoi(items[7].c_str());
            //     if(items[5] != ref || items[4] != strand || abs(s_ctg - prev_ctg) > fuzzy || abs(s_ref - prev_ref) > fuzzy) {
            //         skip = true;
            //         num_skip += 1;
            //         break;
            //     }
            //     prev_ctg = atoi(items[3].c_str());
            //     prev_ref = atoi(items[8].c_str());
            // }
            // if(skip) continue;
            // if(abs(atoi(paf_[ctg][0][1].c_str()) - prev_ctg) > fuzzy) {
            //     num_skip += 1;
            //     continue;
            // }
            std::vector<int> map_pos(atoi(paf_[ctg][0][1].c_str()), -1);
            // prev_ctg = 0;
            int l_ctg = 0, l_ref = 0;
            for(auto item: paf_[ctg]) {
                l_ctg = atoi(item[1].c_str());
                int s_ctg = atoi(item[2].c_str());
                int e_ctg = atoi(item[3].c_str());
                l_ref = atoi(item[6].c_str());
                int s_ref = atoi(item[7].c_str());
                int e_ref = atoi(item[8].c_str());

                int cur_ref = s_ref;
                int cur_ctg = s_ctg;
                if(strand == "-") cur_ctg = e_ctg - 1 ;
                std::string cigar = item[12].substr(5);
                for(std::sregex_iterator iter(cigar.begin(), cigar.end(), pattern), end; iter != end; ++iter) {
                    int size =  atoi(iter->str(1).c_str());
                    std::string op = iter->str(2);
                    if(op == "M" || op == "X" || op == "=") {
                        if(strand == "+") {
                            for(int i = 0; i < size; ++i) {
                                map_pos[cur_ctg + i] = cur_ref + i;
                                // if(cur_ctg + i < prev_ctg) map_pos[cur_ctg + i] = -2;
                                // else map_pos[cur_ctg + i] = cur_ref + i;
                            }
                            cur_ctg += size;
                        } else {
                            for(int i = 0; i < size; ++i) {
                                map_pos[cur_ctg - i] = cur_ref + i;
                                // if(cur_ctg - i - 1 < prev_ctg) map_pos[cur_ctg - i] = -2;
                                // else map_pos[cur_ctg - i] = cur_ref + i;
                            }
                            cur_ctg -= size;
                        }
                        cur_ref += size;
                    } else if(op == "D") cur_ref += size;
                    else if(op == "I") {
                        if(strand == "+") cur_ctg += size;
                        else cur_ctg -= size;
                    }
                }
                // prev_ctg = e_ctg;
                assert(cur_ref = e_ref);
                assert(cur_ctg + 1 == s_ctg || cur_ctg == e_ctg);
            }
            const int BLOCK = opt_.block;
            auto cut = bps[ctg];
            int prev = cut.front();
            int tail = cut.back() - 1;
            if(map_pos[prev] < 0) {
                for(; prev < map_pos.size(); ++prev) {
                    if(map_pos[prev] >= 0) break;
                }
            }
            if(map_pos[tail] < 0) {
                for(; tail > prev; --tail) {
                    if(map_pos[tail] >= 0) break;
                }
            }
            if(tail < prev) continue;
            if(cut.size() == 2) {
                if(strand == "+")
                    cutsites[tid].emplace_back( std::array<std::string, 9>{
                                                ctg, 
                                                std::to_string(l_ctg), 
                                                std::to_string(prev), 
                                                std::to_string(tail + 1), 
                                                ref, 
                                                std::to_string(l_ref), 
                                                std::to_string(map_pos[prev]), 
                                                std::to_string(map_pos[tail] + 1), 
                                                strand} );
                else
                    cutsites[tid].emplace_back( std::array<std::string, 9>{
                                                ctg, 
                                                std::to_string(l_ctg), 
                                                std::to_string(prev), 
                                                std::to_string(tail + 1), 
                                                ref, 
                                                std::to_string(l_ref), 
                                                std::to_string(map_pos[tail]), 
                                                std::to_string(map_pos[prev] + 1), 
                                                strand} );
            } else {
                std::size_t last_i = cut.size() - 2;
                for(; last_i > 0; --last_i) {
                    if(tail - cut[last_i] > BLOCK) break;
                }
                int last_pos = cut[last_i];
                for(; last_pos > prev; --last_pos) if(map_pos[last_pos] != -1) break;
                for(std::size_t i = 1; i < last_i; ++i) {
                    if(cut[i] < prev) continue;
                    if(cut[i] - prev > BLOCK) {
                        int cur = cut[i];
                        for(; cur < tail; ++cur) {
                            if(map_pos[cur] >= 0) break;
                        }
                        if(tail - cur < BLOCK) {
                            cur = cut[i];
                            for(; cur >= prev; --cur) {
                                if(map_pos[cur] >= 0) break;
                            }
                            if(cur - prev < BLOCK) break;
                        }
                        if(strand == "+") 
                            cutsites[tid].emplace_back( std::array<std::string, 9>{
                                                        ctg, 
                                                        std::to_string(l_ctg), 
                                                        std::to_string(prev), 
                                                        std::to_string(cur), 
                                                        ref, 
                                                        std::to_string(l_ref), 
                                                        std::to_string(map_pos[prev]), 
                                                        std::to_string(map_pos[cur]), 
                                                        strand} );
                        else 
                            cutsites[tid].emplace_back( std::array<std::string, 9>{
                                                        ctg, 
                                                        std::to_string(l_ctg), 
                                                        std::to_string(prev), 
                                                        std::to_string(cur), 
                                                        ref, 
                                                        std::to_string(l_ref), 
                                                        std::to_string(map_pos[cur]), 
                                                        std::to_string(map_pos[prev]), 
                                                        strand} );
                        prev = cur;
                    }
                }
                if(prev == cut[last_i] && map_pos[prev] < 0) {
                    for(; prev < tail; ++prev) {
                        if(map_pos[prev] >= 0) break;
                    }
                }
                if(strand == "+") 
                    cutsites[tid].emplace_back( std::array<std::string, 9>{
                                                ctg, 
                                                std::to_string(l_ctg), 
                                                std::to_string(prev), 
                                                std::to_string(tail + 1), 
                                                ref, 
                                                std::to_string(l_ref), 
                                                std::to_string(map_pos[prev]), 
                                                std::to_string(map_pos[tail] + 1), 
                                                strand} );
                else 
                    cutsites[tid].emplace_back( std::array<std::string, 9>{
                                                ctg, 
                                                std::to_string(l_ctg), 
                                                std::to_string(prev), 
                                                std::to_string(tail + 1), 
                                                ref, 
                                                std::to_string(l_ref), 
                                                std::to_string(map_pos[tail]), 
                                                std::to_string(map_pos[prev] + 1), 
                                                strand} );
            }
        }
    };
    MultiThreads(opt_.threads, func);
    std::vector<std::array<std::string, 9>> result;
    for(auto &cs: cutsites) {
        result.insert(result.end(), cs.begin(), cs.end());
        cs.clear();
    }
    std::sort(result.begin(), result.end(), [](std::array<std::string, 9> &a, std::array<std::string, 9> &b) {
        if(a[0] != b[0]) return a[0] < b[0];
        else return atoi(a[1].c_str()) < atoi(b[1].c_str());
    });

    std::string unfiltered_fname = opt_.prefix + ".unfiltered.alt2pri.txt";
    std::ofstream uout(unfiltered_fname);
    if(!uout.is_open()) {
        std::cerr << "[" << GetCurTime() << "] Could not open " << unfiltered_fname << " for writing\n";
        exit(EXIT_FAILURE);
    }
    for(auto items: result) {
        for(auto i: items) {
            uout << i << "\t";
        }
        uout << "\n";
    }
    uout.close();
    std::cerr << "[" << GetCurTime() << "] Dump unfiltered mapping infomation into " << unfiltered_fname << "\n";

    // filter
    std::unordered_map<std::string, std::vector<std::array<std::string, 9>>> pri_unfilt;
    for(auto items: result) {
        pri_unfilt[items[4]].emplace_back(items);
    }
    result.clear();
    for(auto &pri: pri_unfilt) {
        std::sort(pri.second.begin(), pri.second.end(), 
        [](std::array<std::string, 9> &a, std::array<std::string, 9> &b) {
            return atoi(a[6].c_str()) < atoi(a[6].c_str());
        });
        std::vector<int8_t> filt(pri.second.size(), 0);
        for(std::size_t i = 0; i < pri.second.size(); ++i) {
            for(std::size_t j = 0; j < i; ++j) {
                int is = atoi(pri.second[i][6].c_str());
                int ie = atoi(pri.second[i][7].c_str());
                int js = atoi(pri.second[j][6].c_str());
                int je = atoi(pri.second[j][7].c_str());
                if(Nested(is, ie, js, je)) {
                    filt[i] = 1;
                }
                if(Nested(js, je, is, ie)) {
                    filt[j] = 1;
                }
                if(Duplicated(is, ie, js, je)) {
                    filt[i] = 1;
                }
            }
        }
        for(std::size_t i = 0; i < pri.second.size(); ++i) {
            if(filt[i] == 0) result.emplace_back(pri.second[i]);
        }
    }
    std::sort(result.begin(), result.end(), [](std::array<std::string, 9> &a, std::array<std::string, 9> &b) {
        if(a[0] != b[0]) return a[0] < b[0];
        else return atoi(a[1].c_str()) < atoi(b[1].c_str());
    });
    std::string filtered_fname = opt_.prefix + ".filtered.alt2pri.txt";
    std::ofstream fout(filtered_fname);
    if(!fout.is_open()) {
        std::cerr << "[" << GetCurTime() << "] Could not open " << filtered_fname << " for writing\n";
        exit(EXIT_FAILURE);
    }
    for(auto items: result) {
        for(auto i: items) {
            fout << i << "\t";
        }
        fout << "\n";
    }
    fout.close();
    std::cerr << "[" << GetCurTime() << "] Dump filtered mapping infomation into " << filtered_fname << "\n";
}

void Remapping::Run() {
    if(!opt_.paffname.empty() && !opt_.altfname.empty()) {
        Load();
        LoadPaf();
        Remap();
        DumpPaf();
    }
    if(!opt_.breakfname.empty()) {
        if(Empty()) LoadPaf();
        Cut();
    }
}

int main(int argc, char **argv) {
    Options opt;
    if(ParseArgument(argc, argv, opt)) exit(EXIT_FAILURE);
    Remapping rm(opt);
    rm.Run();
}