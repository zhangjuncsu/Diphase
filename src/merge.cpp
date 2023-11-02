#include "merge.hpp"
#include "utility.hpp"

#include <regex>
#include <mutex>
#include <atomic>
#include <unordered_set>

#include "getopt.h"

int Options::Check() {
    if(pfname.empty()) {
        std::cerr << "[" << GetCurTime() << "] Please specify the file name of primary contigs [-P | --pri]\n";
        return 1;
    }
    if(threads < 1) {
        std::cerr << "[" << GetCurTime() << "] threads must be >= 1 [-t | --threads]\n";
        return 1;
    }
    return 0;
}

void Usage(Options &opt) {
    std::cerr << "Usage: merge [options] <alt2pri.sorted.bam> -P <primary.fasta>\n";
    std::cerr << "\tOptions:\n";
    std::cerr << "\t\t-P | --pri     <FILE> file name of primary contigs\n";
    std::cerr << "\t\t-p | --prefix  <STR>  prefix of output files [" << opt.prefix << "]\n";
    std::cerr << "\t\t-t | --threads <INT>  number of threads [" << opt.threads << "]\n";
    std::cerr << "\t\t-h | --help           display this message\n";
}

int ParseArgument(int argc, char **argv, Options &opt) {
    struct option long_opt[] = {
        {"pri",     required_argument,  NULL,   'P'}, 
        {"prefix",  optional_argument,  NULL,   'p'}, 
        {"threads", optional_argument,  NULL,   't'}, 
        {"help",    no_argument,        NULL,   'h'}, 
        {0, 0,  0,  }
    };
    const char *short_opt = "P:p:t:h";
    int c;
    while((c = getopt_long(argc, argv, short_opt, long_opt, NULL)) != -1) {
        if(c == 'P') opt.pfname = optarg;
        else if(c == 'p') opt.prefix = optarg;
        else if(c == 't') opt.threads = atoi(optarg);
        else if(c == 'h') {
            Usage(opt);
            exit(EXIT_SUCCESS);
        } else if(c == ':') {
            std::cerr << "[" << GetCurTime() << "] ERROR missing argument in option " << optopt << "\n";
            return 1;
        } else if(c == '?') {
            std::cerr << "[" << GetCurTime() << "] ERROR unknown option " << optopt << "\n";
            return 1;
        }
    }
    if(optind < argc) opt.bfname = argv[optind];
    return opt.Check();
}
bool debug = false;
std::vector<std::array<int64_t, 2>> Merge::GenerateSubctg(std::vector<int64_t> &poss, std::vector<char> &marks) {
    int num_i = 0, num_c = 0;
    std::vector<std::array<int64_t, 2>> ret;
    std::size_t index = 0;
    if(debug) {
        std::cerr << poss.size() << " " << marks.size() << std::endl;
        for(auto p: poss) {
            std::cerr << p << " ";
        }
        std::cerr << "\n";
        for(auto m: marks) {
            std::cerr << m << " ";
        }
        std::cerr << "\n";
    }
    while(index < marks.size()) {
        std::size_t step = 1;
        if(marks[index] == 'i') {
            if(num_c > opt_.threshold) {
                int num_ii = num_i;
                std::size_t index_i = index;
                while(index_i < marks.size()) {
                    if(marks[index_i] == 'i') {
                        num_ii++;
                    } else if(marks[index_i] == 'c') {
                        break;
                    }
                    index_i++;
                }
                step = index_i - index;
                if(num_ii > opt_.threshold) {
                    assert(marks[index - 1] == 'c' || marks[index - 1] == 's' || marks[index - 1] == 'e');
                    auto last = std::string(marks.begin(), marks.begin() + index).find_last_of('s');
                    ret.emplace_back(std::array<int64_t, 2> { poss[0], poss[last] });
                    if(marks[index - 1] == 's' || marks[index - 1] == 'e') {
                        ret.emplace_back(std::array<int64_t, 2> { poss[index - 1], poss[poss.size() - 1]});
                    } else if(marks[index - 1] == 'c') {
                        ret.emplace_back(std::array<int64_t, 2> { (poss[index - 1] + poss[index]) / 2, poss[poss.size() - 1]});
                    }
                    break;
                }
            }
            num_i += 1;
            num_c = 0;
        } else if(marks[index] == 'c') {
            if(num_i > opt_.threshold) {
                int num_cc = num_c;
                std::size_t index_i = index;
                while(index_i < marks.size()) {
                    if(marks[index_i] == 'c') {
                        num_cc++;
                    } else if(marks[index_i] == 'i') {
                        break;
                    }
                    index_i++;
                }
                step = index_i - index;
                if(num_cc > opt_.threshold) {
                    assert(marks[index - 1] == 'i' || marks[index - 1] == 's' || marks[index - 1] == 'e');
                    if(marks[index - 1] == 's' || marks[index - 1] == 'e') {
                        ret.emplace_back(std::array<int64_t, 2> { poss[0], poss[index - 1]});
                    } else if(marks[index - 1] == 'i') {
                        ret.emplace_back(std::array<int64_t, 2> { poss[0], (poss[index - 1] + poss[index]) / 2});
                    }
                    break;
                }
                index_i = index;
                // while(index_i < marks.size()) {
                //     if(marks[index_i] == 'c') num_cc--;
                //     if(num_cc == 0) break;
                //     index_i++;
                // }
                // if(index_i + 1 < poss.size()) ret.emplace_back(std::array<int64_t, 2> { poss[index_i + 1], poss[poss.size() - 1] });
            }
            num_c += 1;
            num_i = 0;
        }
        if(step != 1) {
            num_i = 0;
            num_c = 0;
        }
        index += step;
    }
    if(ret.empty()) {
        ret.emplace_back(std::array<int64_t, 2>{ poss[0], poss[poss.size() - 1] });
    }
    return ret;
}

std::unordered_map<std::string, std::string> Merge::MergeCtg(std::vector<bam1_t*> &records) {
    std::unordered_set<std::string> used;
    std::unordered_map<std::string, std::vector<std::array<int, 3>>> ctg_info;
    std::vector<std::string> ctg_names;
    std::vector<std::array<int64_t, 2>> intervals;
    // std::cerr << "1" << std::endl;
    
    for(auto &record: records) {
        std::string ref = reader_.Header()->target_name[record->core.tid];
        // if(ref == "ptg000004l") debug = true;
        // else debug = false;
        std::string query = bam_get_qname(record);
        uint32_t *cigar = bam_get_cigar(record);
        if(record->core.qual < 60) continue;
        if(record->core.n_cigar > 0) {
            if(bam_cigar_op(cigar[0]) == BAM_CHARD_CLIP || bam_cigar_op(cigar[record->core.n_cigar - 1]) == BAM_CHARD_CLIP) continue;
        }
        if(used.find(query) != used.end()) continue;
        used.insert(query);
        std::size_t len = rs_.GetSeqLength(ref);
        if(ctg_info[ref].size() != len + 1) {
            std::cerr << ref << std::endl;
            ctg_info[ref].resize(len + 1);
        }
        int64_t ref_pos = record->core.pos;
        int64_t start = ref_pos;
        ctg_info[ref][ref_pos][0] = START;
        for(uint32_t i = 0; i < record->core.n_cigar; ++i) {
            int size = bam_cigar_oplen(cigar[i]);
            int op = bam_cigar_op(cigar[i]);
            if(op == BAM_CMATCH || op == BAM_CEQUAL) {
                for(int j = 0; j < size; ++j) {
                    ctg_info[ref][ref_pos + j][1] += 1;
                }
                ref_pos += size;
            } else if(op == BAM_CDIFF) {
                for(int j = 0; j < size; ++j) {
                    ctg_info[ref][ref_pos + j][2] += 1;
                }
                ref_pos += size;
            } else if(op == BAM_CDEL) {
                ref_pos += size;
            } else if(op == BAM_CINS) {

            } else if(op == BAM_CSOFT_CLIP) {

            }
        }
        ctg_info[ref][ref_pos][0] = END;
        ctg_names.emplace_back(query);
        intervals.emplace_back(std::array<int64_t, 2> { start, ref_pos });
    }
    // std::cerr << "2" << std::endl;
    
    std::unordered_set<std::string> filtered;
    for(std::size_t i = 0; i < intervals.size(); ++i) {
        for(std::size_t j = 0; j < i; ++j) {
            int64_t ast = intervals[i][0];
            int64_t ae = intervals[i][1];
            int64_t bst = intervals[j][0];
            int64_t be = intervals[j][1];
            if(ast == bst && ae == be) {
                filtered.insert(ctg_names[i]);
            }
            if(ast < bst && be < ae) {
                filtered.insert(ctg_names[j]);
            }
            if(bst < ast && ae < be) {
                filtered.insert(ctg_names[i]);
            }
        }
    }
    intervals.clear();
    ctg_names.clear();
    used.clear();
    // std::cerr << "3" << std::endl;
    
    std::unordered_map<std::string, std::string> subctgs;
    std::string prev_ref;
    int64_t prev_end = 0;
    int64_t last_end = 0;
    for(auto &record: records) {
        if(debug) std::cerr << "a" << std::endl;
        std::string ref = reader_.Header()->target_name[record->core.tid];
        std::string query = bam_get_qname(record);
        if(query == "atg007113l") debug = true;
        else debug = false;
        if(filtered.find(query) != filtered.end()) continue;
        if(record->core.qual < 60) continue;
        uint32_t *cigar = bam_get_cigar(record);
        if(record->core.n_cigar > 0) {
            if(bam_cigar_op(cigar[0]) == BAM_CHARD_CLIP || bam_cigar_op(cigar[record->core.n_cigar - 1]) == BAM_CHARD_CLIP) continue;
        }
        if(used.find(query) != used.end()) continue;
        if(debug) std::cerr << "b" << std::endl;
        used.insert(query);
        int64_t start = record->core.pos;
        int64_t end = start;
        if(prev_ref.empty() && ref != prev_ref) {
            prev_ref = ref;
            prev_end = start;
        }
        if(start > prev_end) {
            prev_end = start;
        }
        if(debug) std::cerr << query << std::endl;
        if(debug) std::cerr << "c" << std::endl;
        int beg = 0;
        int end_tmp = beg;
        std::vector<int64_t> pos;
        std::vector<char> marks;
        uint8_t *seq = bam_get_seq(record);
        std::string query_seq;
        query_seq.reserve(record->core.l_qseq);
        for(int32_t i = 0; i < record->core.l_qseq; ++i) {
            query_seq += seq_nt16_str[bam_seqi(seq, i)];
        }
        std::vector<int64_t> mappos(record->core.l_qseq + 1);
        if(debug) std::cerr << "d" << std::endl;
        if(debug) std::cerr << "prev end " << prev_end << " last end " << last_end << std::endl;
        for(uint32_t i = 0; i < record->core.n_cigar; ++i) {
            int size = bam_cigar_oplen(cigar[i]);
            int op = bam_cigar_op(cigar[i]);
            if(op == BAM_CMATCH || op == BAM_CEQUAL) {
                for(int j = 0; j < size; ++j) {
                    mappos[end_tmp + j] = end + j;
                    if(end + j < prev_end) continue;
                    else if(end + j == prev_end) {
                        pos.emplace_back(end_tmp + j);
                        marks.emplace_back('S');
                        continue;
                    }
                    if(ctg_info[ref][end + j][0] == START) {
                        pos.emplace_back(end_tmp + j);
                        marks.emplace_back('s');
                    } else if(ctg_info[ref][end + j][0] == END) {
                        pos.emplace_back(end_tmp + j);
                        marks.emplace_back('e');
                    } else if(ctg_info[ref][end + j][1] > 0 && ctg_info[ref][end + j][2] > 0) {
                        if(ctg_info[ref][end + j][1] < ctg_info[ref][end + i][2]) {
                            pos.emplace_back(end_tmp + j);
                            marks.emplace_back('i');
                        } else {
                            pos.emplace_back(end_tmp + j);
                            marks.emplace_back('c');
                        }
                    }
                }
                end += size;
                end_tmp += size;
            } else if(op == BAM_CDIFF) {
                for(int j = 0; j < size; ++j) {
                    mappos[end_tmp + j] = end + j;
                    if(end + j < prev_end) continue;
                    else if(end + j == prev_end) {
                        pos.emplace_back(end_tmp + j);
                        marks.emplace_back('S');
                        continue;
                    }
                    if(ctg_info[ref][end + j][0] == START) {
                        pos.emplace_back(end_tmp + j);
                        marks.emplace_back('s');
                    } else if(ctg_info[ref][end + j][0] == END) {
                        pos.emplace_back(end_tmp + j);
                        marks.emplace_back('e');
                    } else if(ctg_info[ref][end + j][1] > 0 && ctg_info[ref][end + j][2] > 0) {
                        if(ctg_info[ref][end + j][2] < ctg_info[ref][end + j][1]) {
                            pos.emplace_back(end_tmp + j);
                            marks.emplace_back('c');
                        } else {
                            pos.emplace_back(end_tmp + j);
                            marks.emplace_back('i');
                        }
                    }
                }
                end += size;
                end_tmp += size;
            } else if(op == BAM_CDEL) {
                end += size;
            } else if(op == BAM_CINS) {
                end_tmp += size;
            } else if(op == BAM_CSOFT_CLIP) {
                if(end_tmp == beg) {
                    end_tmp += size;
                } else {
                    end_tmp += size;
                    // pos.emplace_back(end_tmp + size);
                    // marks.emplace_back('E');
                }
            }
        }
        pos.emplace_back(end_tmp);
        marks.emplace_back('E');
        mappos[end_tmp] = end;
        if(debug) std::cerr << "e " << " " << start << " " << query << std::endl;
        auto region = GenerateSubctg(pos, marks);
        if(debug) std::cerr << "f" << std::endl;
        if(region.size() > 2) {
            std::cerr << ref << " " << query << " " << region.size() << std::endl;
            for(auto r: region) std::cerr << r[0] << "-" << r[1] << std::endl;
        }
        assert(region.size() <= 2);
        if(region.size() == 0) continue;
        if(debug) std::cerr << "region size " << region.size() << std::endl;
        // if(region.size() == 2) std::cerr << ref << " " << query << std::endl;
        for(auto r: region) {
            if(debug) std::cerr << r[0] << "-" << r[1] << std::endl;
            std::string head = query + ":" + std::to_string(r[0]) + "-" + std::to_string(r[1]);
            subctgs[head] = query_seq.substr(r[0], r[1] - r[0]);
        }
        if(debug) std::cerr << "g" << std::endl;
        prev_end = mappos[region.front()[1]];
        if(last_end == 0) {
            last_end = mappos[region.back()[1]];
        } else if(last_end > mappos[region.back()[1]]) {
            prev_end = last_end;
        } else {
            last_end = mappos[region.back()[1]];
        }
        if(debug) std::cerr << "h" << std::endl;
    }
    // std::cerr << "4 " << records.size() << std::endl;
    for(auto &record: records) {
        bam_destroy1(record);
    }
    // std::cerr << "5 " << subctgs.size() << std::endl;
    return subctgs;
}

void Merge::MergeCtgs() {
    std::string ofname = opt_.prefix + ".zipped.fasta";
    std::ofstream out(ofname);
    if(!out.is_open()) {
        std::cerr << "[" << GetCurTime() << "] Could not open " << ofname << " for writing\n";
        exit(EXIT_FAILURE);
    }
    
    auto refs = reader_.LoadRefName();
    std::atomic<std::size_t> index { 0 };
    std::mutex mtx;
    auto func = [&](std::size_t tid) {
        for(auto cur = index.fetch_add(1); cur < refs.size(); cur = index.fetch_add(1)) {
            std::string ref = refs[cur];
            auto records = reader_.Load(ref);
            auto seqs = MergeCtg(records);
            {
                std::lock_guard<std::mutex> lock(mtx);
                for(auto &seq: seqs) {
                    out << ">" << seq.first << "\n";
                    out << seq.second << "\n";
                }
            }
        }
    };
    MultiThreads(opt_.threads, func);
    out.close();
}

void Merge::Run() {
    std::ios::sync_with_stdio(false);
    reader_.Initialize(opt_.bfname);
    rs_.Load(opt_.pfname);
    MergeCtgs();
}

int main(int argc, char **argv) {
    Options opt;
    if(ParseArgument(argc, argv, opt) == 1) {
        Usage(opt);
        exit(EXIT_FAILURE);
    } else {
        Merge m(opt);
        m.Run();
    }
    return 0;
}