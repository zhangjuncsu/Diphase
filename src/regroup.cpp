#include "regroup.hpp"
#include "utility.hpp"
#include "matrix.hpp"

#include "getopt.h"

void USAGE(Option &opt) {
    std::cerr << "Usage: regroup -b <bam_file> -o <group> [options]\n";
    std::cerr << "\tOptions:\n";
    std::cerr << "\t\t -b | --bam     <FILE> input bam file\n";
    std::cerr << "\t\t -o | --output  <FILE> output file\n";
    std::cerr << "\t\t -i | --iter    <INT>  number of iterations [" << opt.iter << "]\n";
    std::cerr << "\t\t -s | --seed    <INT>  random seed [" << opt.seed << "]\n";
    std::cerr << "\t\t -t | --threads <INT>  number of threads [" << opt.threads << "]\n";
    std::cerr << "\t\t      --b1      <FILE> bed file of hap1\n";
    std::cerr << "\t\t      --b2      <FILE> bed file of hap2\n";
    std::cerr << "\t\t -h | --help           print this help message\n";
}

int ParseArgument(int argc, char **argv, Option &opt) {
    struct option long_opt[] = {
        {"bam",     required_argument,  0,  'b'},
        {"output",  required_argument,  0,  'o'},
        {"iter",    required_argument,  0,  'i'},
        {"seed",    required_argument,  0,  's'},
        {"threads", required_argument,  0,  't'},
        {"help",    no_argument,        0,  'h'},
        {"b1",      required_argument,  0,  301},
        {"b2",      required_argument,  0,  302},
        {"more",    no_argument,        0,  303},
        {0, 0, 0, 0}
    };
    const char *opt_string = "b:o:t:i:s:p:h";
    int c;
    while((c = getopt_long(argc, argv, opt_string, long_opt, nullptr)) != -1) {
        if(c == 'b') opt.bfname = optarg;
        else if(c == 'o') opt.ofname = optarg;
        else if(c == 'p') opt.pfname = optarg;
        else if(c == 'i') opt.iter = std::stoi(optarg);
        else if(c == 's') opt.seed = std::stoi(optarg);
        else if(c == 't') opt.threads = std::stoi(optarg);
        else if(c == 301) opt.b1fname = optarg;
        else if(c == 302) opt.b2fname = optarg;
        else if(c == 303) opt.more = true;
        else if(c == 'h') {
            USAGE(opt);
            return 0;
        } else if(c == ':') {
            std::cerr << "Missing argument for option: " << optopt << "\n";
            return 1;
        } else if(c == '?') {
            std::cerr << "Unknown option: " << optopt << "\n";
            return 1;
        }
    }
    if(argc == 1) {
        USAGE(opt);
        return 1;
    }
    return 0;
}

void Regroup::Count(const std::string &b1fname, const std::string &b2fname) {
    std::ifstream ifs1(b1fname);
    if(!ifs1.is_open()) {
        LOG(ERROR)("Cannot open file: %s", b1fname.c_str());
    }
    std::string line;
    while(std::getline(ifs1, line)) {
        auto items = SplitString(line, '\t');
        if(items.size() < 4) {
            LOG(ERROR)("Invalid bed file: %s", b1fname.c_str());
        }
        std::string head = items[0] + ":" + items[1] + "-" + items[2];
        if(count_.find(head) == count_.end()) {
            count_[head] = 1;
        } else {
            count_[head] += 1;
        }
    }
    ifs1.close();
    std::ifstream ifs2(b2fname);
    if(!ifs2.is_open()) {
        LOG(ERROR)("Cannot open file: %s", b2fname.c_str());
    }
    while(std::getline(ifs2, line)) {
        auto items = SplitString(line, '\t');
        if(items.size() < 4) {
            LOG(ERROR)("Invalid bed file: %s", b2fname.c_str());
        }
        std::string head = items[0] + ":" + items[1] + "-" + items[2];
        if(count_.find(head) == count_.end()) {
            count_[head] = 1;
        } else {
            count_[head] += 1;
        }
    }
    ifs2.close();
}

void Regroup::LoadPair() {
    std::ifstream ifs(opt_.pfname);
    if(!ifs.is_open()) {
        LOG(ERROR)("Cannot open file: %s", opt_.pfname.c_str());
    }
    std::string line;
    while(std::getline(ifs, line)) {
        auto items = SplitString(line, '\t');
        pairs_.push_back(items[0]);
        pairs_.push_back(items[1]);
    }
    ifs.close();

}

std::vector<std::string> Regroup::LoadBed(const std::string &fname) {
    std::ifstream ifs(fname);
    if(!ifs.is_open()) {
        LOG(ERROR)("Cannot open file: %s", fname.c_str());
    }
    std::vector<std::string> tmp;
    std::string line;
    while(std::getline(ifs, line)) {
        auto items = SplitString(line, '\t');
        if(items.size() < 4) {
            LOG(ERROR)("Invalid bed file: %s", fname.c_str());
        }
        std::string head = items[0] + ":" + items[1] + "-" + items[2];
        if(std::find(tmp.begin(), tmp.end(), items[3]) == tmp.end()) {
            tmp.push_back(items[3]);
        }
        // if(count_[head] == 2) continue;
        if(std::find(pairs_.begin(), pairs_.end(), head) == pairs_.end()) continue;
        mapped_[head] = items[3];
        if(index_.find(items[3]) == index_.end()) {
            index_[items[3]] = names_.size();
            names_.push_back(items[3]);
        }
    }
    ifs.close();
    return tmp;
}

void Regroup::GenerateMatrix() {
    // Matrix matrix(index_.size(), index_.size());
    matrix_.Resize(index_.size(), index_.size());
    BamReader reader;
    reader.Initialize(opt_.bfname);

    std::array<bam1_t*, 2> reads;
    reads[0] = bam_init1();
    reads[1] = bam_init1();
    long count = 0;
    // int mapq = 0;
    while(1) {
        int r = reader.LoadPair(reads);
        if(r < 0) break;
        // if(reads[0]->core.qual < mapq || reads[1]->core.qual < mapq) continue;
        std::string name1 = bam_get_qname(reads[0]);
        std::string name2 = bam_get_qname(reads[1]);
        if(name1 != name2) {
            LOG(ERROR)("Read names are not the same: %s, %s", name1.c_str(), name2.c_str());
        }
        std::string ref1 = reader.Header()->target_name[reads[0]->core.tid];
        std::string ref2 = reader.Header()->target_name[reads[1]->core.tid];
        std::size_t index1 = -1, index2 = -1;
        if(mapped_.find(ref1) == mapped_.end() || mapped_.find(ref2) == mapped_.end()) continue;
        index1 = index_[mapped_[ref1]];
        index2 = index_[mapped_[ref2]];
        if(matrix_.AddLink(index1, index2) < 0) {
            LOG(ERROR)("Failed to add link: %s, %s", ref1.c_str(), ref2.c_str());
        }
        if(matrix_.AddLink(index2, index1) < 0) {
            LOG(ERROR)("Failed to add link: %s, %s", ref2.c_str(), ref1.c_str());
        }
        count += 1;
        if(count % 1000000 == 0) {
            LOG(INFO)("Processed %ld\n", count);
        }
    }
    LOG(INFO)("process %zd reads\n", count);
    bam_destroy1(reads[0]);
    bam_destroy1(reads[1]);
}

void Regroup::GenerateMatrixMore() {
    matrix_.Resize(index_.size(), index_.size());
    BamReader reader;
    reader.Initialize(opt_.bfname);

    int mapq = 1;
    long count = 0;
    std::size_t SIZE = 1000;
    while(1) {
        std::vector<std::vector<bam1_t*>> reads(SIZE);
        std::size_t size = reader.LoadBatchPair(reads);
        if(size == 0) break;
        for(std::size_t k = 0; k < size; ++k) {
            for(std::size_t i = 0; i < reads[k].size(); ++i) {
                if(reads[k][i]->core.qual < mapq) continue;
                std::string namei = bam_get_qname(reads[k][i]);
                std::string refi = reader.Header()->target_name[reads[k][i]->core.tid];
                for(std::size_t j = 0; j < reads[k].size(); ++j) {
                    if(reads[k][j]->core.qual < mapq) continue;
                    if(i == j) continue;
                    std::string namej = bam_get_qname(reads[k][j]);
                    std::string refj = reader.Header()->target_name[reads[k][j]->core.tid];
                    if(namei != namej) {
                        LOG(ERROR)("Read names are not the same: %s, %s", namei.c_str(), namej.c_str());
                    }
                    std::size_t indexi = -1, indexj = -1;
                    indexi = index_[mapped_[refi]];
                    indexj = index_[mapped_[refj]];
                    if(matrix_.AddLink(indexi, indexj) < 0) {
                        LOG(ERROR)("Failed to add link: %s, %s", refi.c_str(), refj.c_str());
                    }
                    if(matrix_.AddLink(indexj, indexi) < 0) {
                        LOG(ERROR)("Failed to add link: %s, %s", refj.c_str(), refi.c_str());
                    }
                }
            }
            count += 1;
            if(count % 1000000 == 0) {
                LOG(INFO)("Processed %ld reads\n", count);
            }
        }
    }
}

void Regroup::LoadGroup() {
    std::vector<std::string> tmp;
    for(auto &it: index_) {
        tmp.push_back(it.first);
    }
    std::sort(tmp.begin(), tmp.end());
    if(tmp.size() % 2 != 0) {
        LOG(ERROR)("Invalid number of groups: %ld\n", tmp.size());
    }
    for(std::size_t i = 0; i < tmp.size(); i += 2) {
        Group group;
        group.first = index_[tmp[i]];
        group.second = index_[tmp[i + 1]];
        if(SplitString(tmp[i], '_')[0] != SplitString(tmp[i + 1], '_')[0]) {
            LOG(ERROR)("Invalid group: %s, %s", tmp[i].c_str(), tmp[i + 1].c_str());
        }
        // int r = rand() & 1;
        // group.label_first = r;
        // group.label_second = r ^ 1;
        group.label_first = 0;
        group.label_second = 1;
        groups_.push_back(group);
    }
    for(std::size_t i = 0; i < groups_.size(); ++i) {
        for(std::size_t j = 0; j < groups_.size(); ++j) {
            if(j == i) continue;
            int c1 = (int)matrix_.Get(groups_[i].first, groups_[j].first);
            int c2 = (int)matrix_.Get(groups_[i].first, groups_[j].second);
            int c3 = (int)matrix_.Get(groups_[i].second, groups_[j].first);
            int c4 = (int)matrix_.Get(groups_[i].second, groups_[j].second);
            if((c1 < 5 && c2 < 5 && c3 < 5 && c4 < 5)) {
                matrix_.Set(groups_[i].first, groups_[j].first, 0);
                matrix_.Set(groups_[i].first, groups_[j].second, 0);
                matrix_.Set(groups_[i].second, groups_[j].first, 0);
                matrix_.Set(groups_[i].second, groups_[j].second, 0);
            }
        }
    }
}

int Regroup::GetMax() const {
    int max = 0;
    for(std::size_t i = 0; i < groups_.size(); ++i) {
        for(std::size_t j = 0; j < groups_.size(); ++j) {
            if(j == i) continue;
            int c1 = (int)matrix_.Get(groups_[i].first, groups_[j].first);
            int c2 = (int)matrix_.Get(groups_[i].first, groups_[j].second);
            int c3 = (int)matrix_.Get(groups_[i].second, groups_[j].first);
            int c4 = (int)matrix_.Get(groups_[i].second, groups_[j].second);
            if(groups_[i].label_first == groups_[j].label_first) {
                if(c1 + c4 > c2 + c3) max += 1;
            } else {
                if(c2 + c3 > c1 + c4) max += 1;
            }
        }
    }
    return max;
}

void Regroup::LocalMaximum() {
    for(std::size_t i = 0; i < groups_.size(); ++i) {
        int count = 0;
        for(std::size_t j = 0; j < groups_.size(); ++j) {
            if(j == i) continue;
            int c1 = (int)matrix_.Get(groups_[i].first, groups_[j].first);
            int c2 = (int)matrix_.Get(groups_[i].first, groups_[j].second);
            int c3 = (int)matrix_.Get(groups_[i].second, groups_[j].first);
            int c4 = (int)matrix_.Get(groups_[i].second, groups_[j].second);
            if(groups_[i].label_first == groups_[j].label_first) {
                if(c1 + c4 > c2 + c3) count -= 1;
                else if(c2 + c3 > c1 + c4) count += 1;
            } else {
                if(c1 + c4 > c2 + c3) count += 1;
                else if(c2 + c3 > c1 + c4) count -= 1;
            }
        }
        if(count > 0) groups_[i].label_first ^= 1;
    }
}

void Regroup::Phase() {
    std::array<int, 2> all { 0, 0 };
    all[0] = GetMax();
    int best = 0;
    for(int iter = 0; iter < opt_.iter; ++iter) {
        while(1) {
            LocalMaximum();
            all[1] = GetMax();
            // if(iter == 0) {
            //     std::cerr << "all[0]: " << all[0] << " all[1]: " << all[1] << std::endl;
            // }
            if(all[1] > all[0]) all[0] = all[1];
            else break;
        }
        // LOG(INFO)("Iteration: %d, max: %d\n", iter, all[0]);
        if(all[0] > best) {
            for(std::size_t i = 0; i < groups_.size(); ++i) {
                groups_[i].label_second = groups_[i].label_first;
            }
            best = all[0];
        }
        // for(std::size_t t = 0; t < 100; ++t) {
            std::size_t random = rand() % groups_.size();
            groups_[random].label_first ^= 1;
            for(std::size_t i = 0; i < groups_.size(); ++i) {
                if(i == random) continue;
                int c1 = (int)matrix_.Get(groups_[random].first, groups_[i].first);
                int c2 = (int)matrix_.Get(groups_[random].first, groups_[i].second);
                int c3 = (int)matrix_.Get(groups_[random].second, groups_[i].first);
                int c4 = (int)matrix_.Get(groups_[random].second, groups_[i].second);
                if(c1 + c2 + c3 + c4 > 0) groups_[i].label_first ^= 1;
            }
        // }
        all[0] = GetMax();
        all[1] = 0;
    }
    for(std::size_t i = 0; i < groups_.size(); ++i) {
        groups_[i].label_first = groups_[i].label_second;
    }
    std::ofstream ofs(opt_.ofname);
    if(!ofs.is_open()) {
        LOG(ERROR)("Cannot open file: %s\n", opt_.ofname.c_str());
    }
    // LOG(INFO)("Best: %d\n", best);
    std::vector<std::string> tmp_name;
    for(std::size_t i = 0; i < groups_.size(); ++i) {
        tmp_name.emplace_back(names_[groups_[i].first]);
        tmp_name.emplace_back(names_[groups_[i].second]);
        if(groups_[i].label_first == 0) {
            ofs << names_[groups_[i].first] << "\t" << names_[groups_[i].second] << "\n";
        } else {
            ofs << names_[groups_[i].second] << "\t" << names_[groups_[i].first] << "\n";
        }
    }
    for(std::size_t i = 0; i < hap1_.size(); ++i) {
        if(std::find(tmp_name.begin(), tmp_name.end(), hap1_[i]) == tmp_name.end()) {
            ofs << hap1_[i] << "\t" << hap2_[i] << "\n";
        }
    }
    ofs.close();
    // for(std::size_t i = 0; i < groups_.size(); ++i) {
    //     std::cout << names_[groups_[i].first] << " " << names_[groups_[i].second] << " ";
    // }
    // std::cout << std::endl;
    // for(std::size_t i = 0; i < groups_.size(); ++i) {
    //     for(std::size_t j = 0; j < groups_.size(); ++j) {
    //         int c1 = (int)matrix_.Get(groups_[i].first, groups_[j].first);
    //         int c2 = (int)matrix_.Get(groups_[i].first, groups_[j].second);
    //         std::cout << c1 << " " << c2 << " ";
    //     }
    //     std::cout << std::endl;
    //     for(std::size_t j = 0; j < groups_.size(); ++j) {
    //         int c3 = (int)matrix_.Get(groups_[i].second, groups_[j].first);
    //         int c4 = (int)matrix_.Get(groups_[i].second, groups_[j].second);
    //         std::cout << c3 << " " << c4 << " ";
    //     }
    //     std::cout << std::endl;
    // }
}

void Regroup::Run() {
    srand(opt_.seed);
    Count(opt_.b1fname, opt_.b2fname);
    LoadPair();
    hap1_ = LoadBed(opt_.b1fname);
    hap2_ = LoadBed(opt_.b2fname);
    if(opt_.more) GenerateMatrixMore();
    else GenerateMatrix();
    LoadGroup();
    Phase();
}

int main(int argc, char **argv) {
    Option opt;
    if(ParseArgument(argc, argv, opt) != 0) return 1;
    Regroup regroup(opt);
    regroup.Run();
    return 0;
}