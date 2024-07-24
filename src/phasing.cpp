#include "phasing.hpp"
#include "utility.hpp"
#include "logger.hpp"

#include <getopt.h>

#include <unordered_set>
#include <iostream>
#include <fstream>
#include <atomic>
#include <mutex>

int Options::Check() {
    if(cfname.empty()) {
        LOG(WARNING)("Please specify a contig file name [-c | --ctg]\n");
        // std::cerr << "[" << GetCurTime() << "] Please specify a contig file name [-c | --ctg]\n";
        return 1;
    }
    if(mfname.empty()) {
        LOG(WARNING)("Please specify a binary matrix file name [-b | --binmat]\n");
        // std::cerr << "[" << GetCurTime() << "] Please specify a binary matrix file name [-b | --binmat]\n";
        return 1;
    }
    if(pfname.empty()) {
        LOG(WARNING)("Please specify a pair file name [-r | --pair]\n");
        // std::cerr << "[" << GetCurTime() << "] Please specify a pair file name [-r | --pair]\n";
        return 1;
    }
    if(threads < 1) {
        LOG(WARNING)("threads must be >= 1 [-t | --threads]\n");
        // std::cerr << "[" << GetCurTime() << "] threads must be >= 1 [-t | --threads]\n";
        return 1;
    }
    if(iteration < 0) {
        LOG(WARNING)("iteration must be >= 0 [-i | --iteration]\n");
        // std::cerr << "[" << GetCurTime() << "] iteration must be >= 0 [-i | --iteration]\n";
        return 1;
    }
    return 0;
}

void USAGE(Options &opt) {
    std::cerr << "Usage:\t";
    std::cerr << "phase [options] -c contig.fasta -b binmat -r pair.txt\n\n";
    std::cerr << "\tOptions:\n";
    std::cerr << "\t\t-c | --ctg       <FILE> congit file name\n";
    std::cerr << "\t\t-b | --binmat    <FILE> binary matrix file name\n";
    std::cerr << "\t\t-r | --pair      <FILE> pair file name\n";
    std::cerr << "\t\t-p | --prefix    <STR>  prefix of output file [" << opt.prefix << "]\n";
    std::cerr << "\t\t-s | --seed      <INT>  seed for random function [" << opt.seed << "]\n";
    std::cerr << "\t\t-t | --threads   <INT>  number of threads [" << opt.threads << "]\n";
    std::cerr << "\t\t-i | --iteration <INT>  iterations for phasing [" << opt.iteration << "]\n";
    std::cerr << "\t\t-h | --help             display this message\n";
}

int ParseArgument(int argc, char **argv, Options &opt) {
    struct option long_opt[] = {
        {"ctg",         required_argument,  NULL,   'c'}, 
        {"binmat",      required_argument,  NULL,   'b'}, 
        {"pair",        required_argument,  NULL,   'r'}, 
        {"prefix",      optional_argument,  NULL,   'p'}, 
        {"seed",        optional_argument,  NULL,   's'}, 
        {"threads",     optional_argument,  NULL,   't'}, 
        {"iteration",   optional_argument,  NULL,   'i'}, 
        {"help",        no_argument,        NULL,   'h'}, 
        {0, 0,  0,  0}
    };
    const char *short_opt = "c:b:r:p:s:t:i:h";
    int c;
    while((c = getopt_long(argc, argv, short_opt, long_opt, NULL)) != -1) {
        if(c == 'c') opt.cfname = optarg;
        else if(c == 'b') opt.mfname = optarg;
        else if(c == 'r') opt.pfname = optarg;
        else if(c == 'p') opt.prefix = optarg;
        else if(c == 's') opt.seed = atoi(optarg);
        else if(c == 'i') opt.iteration = atoi(optarg);
        else if(c == 't') opt.threads = atoi(optarg);
        else if(c == 'h') {
            USAGE(opt);
            exit(EXIT_SUCCESS);
        } else if(c == ':') {
            LOG(WARNING)("missing option argument in %c\n", optopt);
            // std::cerr << "[" << GetCurTime() << "] ERROR missing option argument in " << optopt << "\n";
            return 1;
        } else if(c == '?') {
            LOG(WARNING)("unknown option in %c\n", optopt);
            // std::cerr << "[" << GetCurTime() << "] ERROR unknown option in " << optopt << "\n";
            return 1;
        }
    }
    return opt.Check();
}

std::vector<std::array<std::string, 2>> &Phasing::LoadPair() {
    std::ifstream pin(opt_.pfname);
    if(!pin.is_open()) {
        LOG(ERROR)("Could not open %s for reading\n", opt_.pfname.c_str());
        // std::cerr << "[" << GetCurTime() << "] Could not open " << opt_.pfname << " for reading\n";
        // exit(EXIT_FAILURE);
    }
    std::vector<std::array<std::string, 2>> pairs;
    std::string line;
    while(std::getline(pin, line)) {
        auto items = SplitString(line, '\t');
        pairs.emplace_back(std::array<std::string, 2>{items[0], items[1]});
    }
    return pairs;
}

void Phasing::GroupCluster(std::vector<std::array<std::string, 2>> &pairs) {
    std::string prev_pri;
    Group g;
    for(auto p: pairs) {
        auto items = SplitString(p[0], ':');
        if(prev_pri.empty() || items[0] != prev_pri) {
            if(!prev_pri.empty()) group_.push_back(g);
            prev_pri = items[0];
            g.num_ctg = 0;
            g.group_name = items[0];
            g.ovlps.clear();
        }
        g.num_ctg++;
        Overlap ovlp;
        ovlp.first = rs.GetIdByName(items[0]);
        ovlp.second = rs.GetIdByName(items[1]);
        int init = rand() & 1;
        ovlp.label = init;
        ovlp.status = INIT;
        g.ovlps.emplace_back(ovlp);
    }
}

void Phasing::RandomPhasing() {
    srand(opt_.seed);
    int thres = 5;
    std::ofstream out(opt_.prefix + ".phase.out");
    std::mutex mtx;
    std::atomic<std::size_t> index { 0 };
    auto func = [&] (std::size_t tid) {
        for(auto cur = index.fetch_add(1); cur < group_.size(); cur = index.fetch_add(1)) {
            Group &g = group_[cur];
            auto &ovlp = g.ovlps;
            int num_strong = 0;
            for(int iter = 0; iter < opt_.iteration; ++iter) {
                num_strong = 0;
                for(std::size_t i = 0; i < ovlp.size(); ++i) {
                    double d1 = 0.000001, d2 = 0.000001;
                    for(std::size_t j = 0; j < ovlp.size(); ++j) {
                        if(i == j) continue;
                        if(ovlp[i].label == ovlp[j].label) {
                            d1 += matrix_.Get(ovlp[i].first, ovlp[j].first);
                            d1 += matrix_.Get(ovlp[i].second, ovlp[j].second);
                            d2 += matrix_.Get(ovlp[i].first, ovlp[j].second);
                            d2 += matrix_.Get(ovlp[i].second, ovlp[j].first);
                        } else {
                            d1 += matrix_.Get(ovlp[i].first, ovlp[j].second);
                            d1 += matrix_.Get(ovlp[i].second, ovlp[j].first);
                            d2 += matrix_.Get(ovlp[i].first, ovlp[j].first);
                            d2 += matrix_.Get(ovlp[i].second, ovlp[j].second);
                        }
                    }
                    if(d1 >= d2 * thres) {
                        ovlp[i].status = STRONG;
                        num_strong++;
                    } else if(d2 >= d1 * thres) {
                        ovlp[i].label ^= 1;
                        ovlp[i].status = STRONG;
                        num_strong++;
                    } else {
                        ovlp[i].status = WEAK;
                    }
                }
                if(num_strong == ovlp.size()) break;
            }
            {
                std::lock_guard<std::mutex> lock(mtx);
                std::cerr << "[" << GetCurTime() << "] " << num_strong << "/" << ovlp.size() << " contigs in " << g.group_name << " are clustered strongly\n";
                for(std::size_t i = 0; i < ovlp.size(); ++i) {
                    if(ovlp[i].label == 0) {
                        out << rs.GetNameById(ovlp[i].first) << "\t" << rs.GetNameById(ovlp[i].second) << "\t" << ovlp[i].first << "\t" << ovlp[i].second << "\t" << ovlp[i].status << "\n";
                    } else {
                        out << rs.GetNameById(ovlp[i].second) << "\t" << rs.GetNameById(ovlp[i].first) << "\t" << ovlp[i].second << "\t" << ovlp[i].first << "\t" << ovlp[i].status << "\n";
                    }
                }
                LOG(INFO)("Group %s has been phased\n", g.group_name.c_str());
                // std::cerr << "[" << GetCurTime() << "] Group " << g.group_name << " has been phased\n";
            }
        }
    };
    MultiThreads(opt_.threads, func);
}

void Phasing::Run() {
    rs.Load(opt_.cfname);
    matrix_.Load(opt_.mfname);
    auto pair = LoadPair();
    GroupCluster(pair);
    RandomPhasing();
}

int main(int argc, char **argv) {
    Options opt;
    if(ParseArgument(argc, argv, opt) == 1) exit(EXIT_FAILURE);
    Phasing p(opt);
    p.Run();
}