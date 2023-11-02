#include "coverage.hpp"
#include "utility.hpp"
#include "read_store.hpp"

#include "getopt.h"

#include <sys/stat.h>

#include <regex>
#include <numeric>
#include <fstream>
#include <iostream>
#include <algorithm>

void Usage(Options &opt) {
    std::cerr << "Usage: cov [options] in.bam -r ref -o output\n";
    std::cerr << "\toptions:\n";
    std::cerr << "\t\t-o | --output  <FILE> output file name\n";
    std::cerr << "\t\t-r | --ref     <FILE> reference name\n";
    std::cerr << "\t\t-c | --ctg     <STR>  output the coverage of STR. if not set, output all\n";
    std::cerr << "\t\t-t | --threads <INT>  number of threads [" << opt.threads << "]\n";
    std::cerr << "\t\t-h | --help           display this message\n";
}

int Options::Check() {
    // if(rfname.empty()) {
    //     std::cerr << "[" << GetCurTime() << "] Please specify a reference file name [-r | --ref]\n";
    //     return 1;
    // }
    if(threads < 1) {
        std::cerr << "[" << GetCurTime() << "] threads should be >= 1 [-t | --threads]\n";
        return 1;
    }
    return 0;
}

int ParseArgument(int argc, char **argv, Options &opt) {
    struct option long_opt[] = {
        {"output",  required_argument,  NULL,   'o'}, 
        {"ref",     required_argument,  NULL,   'r'}, 
        {"ctg",     optional_argument,  NULL,   'c'}, 
        {"var",     required_argument,  NULL,   'v'}, 
        {"threads", optional_argument,  NULL,   't'}, 
        {"help",    no_argument,        NULL,   'h'}, 
        {0, 0,  0,  0}
    };
    const char *short_opt = "o:r:c:v:t:h";
    int c;
    while((c = getopt_long(argc, argv, short_opt, long_opt, NULL)) != -1) {
        if(c == 'o') opt.ofname = optarg;
        else if(c == 'r') opt.rfname = optarg;
        else if(c == 'c') opt.cname = optarg;
        else if(c == 'v') opt.vfname = optarg;
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

void Coverage::Init() {
    ReadStore rs;
    rs.Load(opt_.rfname);
    for(std::size_t i = 0; i < rs.Size(); ++i) {
        cov_[rs.GetNameById(i)].resize(rs.GetSeqLength(i));
    }
}

void Coverage::LoadVariants() {
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

int Coverage::FilterHiCPairWithSNP(std::array<bam1_t*, 2> &pair) {
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

void Coverage::GetCoverage() {
    for(int32_t i = 0; i < reader_.Header()->n_targets; ++i) {
        cov_[reader_.Header()->target_name[i]].resize(reader_.Header()->target_len[i] + 1, 0);
    }
    std::size_t SIZE = 100000;
    std::vector<std::array<bam1_t*, 2>> pairs(SIZE);
    for(std::size_t i = 0; i < SIZE; ++i) {
        pairs[i][0] = bam_init1();
        pairs[i][1] = bam_init1();
    }
    std::size_t count = 0;
    while(1) {
        std::size_t size = 0;
        for(; size < SIZE; ++size) {
            if(reader_.LoadPair(pairs[size]) < 0) break;
        }
        if(size == 0) break;
        count += size;
        for(std::size_t i = 0; i < size; ++i) {
            if(FilterHiCPairWithSNP(pairs[i]) == 1) continue;
            std::string rname1 = reader_.Header()->target_name[pairs[i][0]->core.tid];
            std::string rname2 = reader_.Header()->target_name[pairs[i][1]->core.tid];
            if(rname1 != rname2) continue;
            if(pairs[i][0]->core.tid != pairs[i][1]->core.tid) continue;
            if(!opt_.cname.empty() && rname1.substr(0, opt_.cname.size()) != opt_.cname) continue;
            std::size_t st1 = pairs[i][0]->core.pos;
            std::size_t en1 = st1 + bam_cigar2rlen(pairs[i][0]->core.n_cigar, bam_get_cigar(pairs[i][0]));
            std::size_t st2 = pairs[i][1]->core.pos;
            std::size_t en2 = st2 + bam_cigar2rlen(pairs[i][1]->core.n_cigar, bam_get_cigar(pairs[i][1]));
            
            std::size_t pos1 = pairs[i][0]->core.pos;
            std::size_t pos2 = pairs[i][1]->core.pos;
            if(pos1 < pos2) {
                pos2 += bam_cigar2rlen(pairs[i][1]->core.n_cigar, bam_get_cigar(pairs[i][1]));
                // std::cerr << pos1 << " " << pos2 << " " << cov_[rname1].size() << std::endl;
                cov_[rname1][pos1] += 1;
                cov_[rname2][pos2] -= 1;
            } else {
                pos1 += bam_cigar2rlen(pairs[i][0]->core.n_cigar, bam_get_cigar(pairs[i][0]));
                cov_[rname1][pos2] += 1;
                cov_[rname1][pos1] -= 1;
            }
        }
        if(count % 1000000 == 0) {
            std::cerr << "[" << GetCurTime() << "] Parse " << count << " pairs\n";
        }
    }
    
    std::cerr << "[" << GetCurTime() << "] Parsed " << count << " pairs\n";
    for(auto &ctg: cov_) {
        for(std::size_t i = 1; i < ctg.second.size(); ++i) {
            ctg.second[i] += ctg.second[i - 1];
        }
    //     auto cov = ctg.second;
    //     if(*std::max_element(cov.begin(), cov.end()) < 10) continue;
    //     double average = std::accumulate(cov.begin(), cov.end(), 0.0) / cov.size();
    //     std::vector<long> positions;
    //     for(int div = 5; div <= 15; ++div) {
    //         double cutoff = average / div;
    //         std::vector<int> delta;
    //         for(std::size_t i = 0; i < cov.size(); ++i) {
    //             if(cov[i] < cutoff) {
    //                 delta.push_back(10);
    //             } else {
    //                 delta.push_back(-10);
    //             }
    //         }
    //         std::size_t sz = delta.size();
    //         int max_so_far = -100000, max_ending_here = 0;
    //         long start = sz / 1000, end = sz / 1000, s = sz / 1000;
    //         for(long i = sz / 1000; i < 99 * sz / 100; ++i) {
    //             max_ending_here += delta.at(i);
    //             if(max_so_far < max_ending_here) {
    //                 max_so_far = max_ending_here;
    //                 start = s;
    //                 end = i;
    //             }
    //             if(max_ending_here < 0) {
    //                 max_ending_here = 0;
    //                 s = i + 1;
    //             }
    //         }
    //         if(start >= sz / 1000 + 50000 && end <= 99 * sz / 100 - 50000) {
    //             // long p = (start + end) / 2;
    //             // long st, en;
    //             // for(st = p; st > 0; --st) {
    //             //     if(cov[st] > cov[p]) break;
    //             // }
    //             // for(en = p; en < cov.size(); ++en) {
    //             //     if(cov[en] > cov[p]) break;
    //             // }
    //             // positions.push_back((st + en) / 2);
    //             positions.push_back((start + end) / 2);
    //         }
    //     }
    //     if(positions.size() > 0) {
    //         int max = *std::max_element(cov.begin(), cov.end());
    //         long consensus_pos = std::accumulate(positions.begin(), positions.end(), 0.0) / positions.size();
    //         std::size_t st, en;
    //         for(st = consensus_pos; st > 0; --st) {
    //             if(cov[st] > cov[consensus_pos] + 10) break;
    //         }
    //         for(en = consensus_pos; en < cov.size(); ++en) {
    //             if(cov[en] > cov[consensus_pos] + 10) break;
    //         }
    //         std::cout << ctg.first << "\t" << consensus_pos << "\t" << cov[consensus_pos] << "/" << max << " ";
    //         std::cout << st << "-" << en << " " << cov[st] << " " << cov[en] << "\n";
    //         // for(auto p: positions) {
    //         //     std::cerr << p << " ";
    //         // }
    //         std::cerr << "\n";
    //     }
    }

    if(mkdir(opt_.ofname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1) {
        std::cerr << "[" << GetCurTime() << "] Could not create folder " << opt_.ofname << "\n";
        exit(EXIT_FAILURE);
    }
    for(auto ctg: cov_) {
        std::string fname = opt_.ofname + "/" + ctg.first + ".cov";
        std::ofstream out(fname);
        if(!out.is_open()) {
            std::cerr << "[" << GetCurTime() << "] Could not open " << fname << " for writing\n";
            exit(EXIT_FAILURE);
        }
        for(auto c: ctg.second) {
            out << c <<"\n";
        }
        out.close();
    }

    for(std::size_t i = 0; i < SIZE; ++i) {
        bam_destroy1(pairs[i][0]);
        bam_destroy1(pairs[i][1]);
    }
}

void Coverage::Run() {
    std::ios::sync_with_stdio(false);
    reader_.Initialize(opt_.bfname);
    LoadVariants();
    GetCoverage();
}

int main(int argc, char **argv) {
    Options opt;
    if(ParseArgument(argc, argv, opt) == 1) {
        Usage(opt);
        exit(EXIT_FAILURE);
    } else {
        Coverage c(opt);
        c.Run();
    }
    return 0;
}