#include "bam2bed.hpp"
#include "bam_reader.hpp"
#include "file_io.hpp"
#include "utility.hpp"
#include "logger.hpp"

#include <regex>
#include <atomic>
#include <cassert>

#include "getopt.h"

void Bam2Bed::LoadBed() {
    for(auto fname: { opt_.b1fname, opt_.b2fname }) {
        GzFileReader reader(fname);
        if(reader.Valid()) {
            std::string prev = "";
            long end = 0;
            std::string line = reader.GetNoEmptyLine();
            while(!line.empty()) {
                auto items = SplitString(line, '\t');
                assert(items.size() == 4);
                if(prev == "" || items[3] != prev) {
                    prev = items[3];
                    end = 0;
                }
                std::string ctg = items[0] + ":" + items[1] + "-" + items[2];
                map_[ctg][0] = items[3];
                map_[ctg][1] = std::to_string(end);
                // map_[ctg][items[3]] = end;
                end += atol(items[2].c_str()) - atol(items[1].c_str());
                line = reader.GetNoEmptyLine();
            }
        } else {
            LOG(ERROR)("Failed to open %s for reading\n", fname.c_str());
            // std::cerr << "[" << GetCurTime() << "] Failed to open " << fname << " for reading\n";
            // exit(EXIT_FAILURE);
        }
    }
}

void Bam2Bed::ToBed() {
    BamReader reader;
    reader.Initialize(opt_.hfname);
    std::size_t SIZE = 100000;
    std::vector<std::array<bam1_t*, 2>> reads(SIZE);
    for(std::size_t i = 0; i < SIZE; ++i) {
        reads[i][0] = bam_init1();
        reads[i][1] = bam_init1();
    }
    std::vector<std::array<std::string, 6>> bed;
    std::regex pattern("(\\w+):(\\d+)-(\\d+)");
    std::smatch result;
    // std::vector<std::array<std::string, 6>> tmp(SIZE);
    while(1) {
        std::size_t size = 0;
        for(std::size_t i = 0; i < SIZE; ++i) {
            int r = reader.LoadPair(reads[i]);
            if(r < 0) break;
            size += 1;
        }
        if(size == 0) break;
        std::vector<std::array<std::string, 6>> tmp(size);
        std::atomic<std::size_t> index { 0 };
        auto func = [&](std::size_t tid) {
            for(auto cur = index.fetch_add(1); cur < size; cur = index.fetch_add(1)) { 
                std::string name1 = reader.Header()->target_name[reads[cur][0]->core.tid];
                long st1 = atol(map_[name1][1].c_str()) + reads[cur][0]->core.pos;
                long en1 = st1 + bam_cigar2rlen(reads[cur][0]->core.n_cigar, bam_get_cigar(reads[cur][0]));
                std::string name2 = reader.Header()->target_name[reads[cur][1]->core.tid];
                long st2 = atol(map_[name2][1].c_str()) + reads[cur][1]->core.pos;
                long en2 = st2 + bam_cigar2rlen(reads[cur][1]->core.n_cigar, bam_get_cigar(reads[cur][1]));
                tmp[cur][0] = map_[name1][0];
                tmp[cur][1] = std::to_string(st1);
                tmp[cur][2] = std::to_string(en1);
                tmp[cur][3] = map_[name2][0];
                tmp[cur][4] = std::to_string(st2);
                tmp[cur][5] = std::to_string(en2);
            }
        };
        MultiThreads(opt_.threads, func);
        bed.insert(bed.end(), tmp.begin(), tmp.end());
    }
    for(std::size_t i = 0; i < SIZE; ++i) {
        bam_destroy1(reads[i][0]);
        bam_destroy1(reads[i][1]);
    }
    std::string fname = opt_.prefix + ".bed";
    std::ofstream out(fname);
    if(!out.is_open()) {
        LOG(ERROR)("Failed to open %s for writing\n", fname.c_str());
        // std::cerr << "[" << GetCurTime() << "] could not open " << fname << " for writing\n";
        // exit(EXIT_FAILURE);
    }
    for(auto a: bed) {
        for(auto i: a) out << i << "\t";
        out << "\n";
    }
    out.close();
}

void Bam2Bed::Run() {
    std::ios::sync_with_stdio(false);
    LoadBed();
    ToBed();
}

int Option::Check() {
    if(b1fname.empty() || b2fname.empty()) {
        LOG(WARNING)("Missing one or two bed file(s)\n");
        // std::cerr << "[" << GetCurTime() << "] Missing one or two bed file(s)\n";
        return 1;
    }
    if(hfname.empty()) {
        LOG(WARNING)("Missing Hi-C mapping file\n");
        // std::cerr << "[" << GetCurTime() << "] Missing Hi-C mapping file\n";
        return 1;
    }
    if(threads < 1) {
        LOG(WARNING)("threads should be >= 1\n");
        // std::cerr << "[" << GetCurTime() << "] threads should be >= 1\n";
        return 1;
    }
    return 0;
}

void Usage(Option &opt) {
    std::cerr << "Usage:\tbam2bed -1 <phased0.bed> -2 <phased1.bed> -b <hic.bam>\n";
    std::cerr << "\t\t-1 | --bed1       <FILE> phased0 bed file name\n";
    std::cerr << "\t\t-2 | --bed2       <FILE> phased1 bed file name\n";
    std::cerr << "\t\t-b | --bam        <FILE> Hi-C mapping file\n";
    std::cerr << "\t\t-p | --prefix     <STR>  prefix of output file [" << opt.prefix << "]\n";
    std::cerr << "\t\t-t | --threads    <INT>  number of threads [" << opt.threads << "]\n";
    std::cerr << "\t\t-h | --help              display this message\n";
}

int ParseArgument(int argc, char **argv, Option &opt) {
    struct option long_opt[] = {
        {"bed1",    required_argument,  NULL,   '1'}, 
        {"bed2",    required_argument,  NULL,   '2'}, 
        {"bam",     required_argument,  NULL,   'b'}, 
        {"prefix",  required_argument,  NULL,   'p'}, 
        {"threads", required_argument,  NULL,   't'}, 
        {"help",    no_argument,        NULL,   'h'}, 
        {0, 0,  0,  0}
    };
    const char *short_opt = "1:2:b:p:t:h";
    int c;
    while((c = getopt_long(argc, argv, short_opt, long_opt, NULL)) != -1) {
        if(c == '1') opt.b1fname = optarg;
        else if(c == '2') opt.b2fname = optarg;
        else if(c == 'b') opt.hfname = optarg;
        else if(c == 'p') opt.prefix = optarg;
        else if(c == 't') opt.threads = atoi(optarg);
        else if(c == 'h') { Usage(opt); exit(EXIT_SUCCESS); }
        else if(c == ':') {
            LOG(WARNING)("Missing argument in option %c\n", optopt);
            // std::cerr << "[" << GetCurTime() << "] ERROR missing argument in option " << optopt << "\n";
            return 1;
        } else if(c == '?') {
            LOG(WARNING)("Unknown option %c\n", optopt);
            // std::cerr << "[" << GetCurTime() << "] ERROR unknown option " << optopt << "\n";
            return 1;
        }
    }
    return opt.Check();
}

int main(int argc, char **argv) {
    Option opt;
    if(ParseArgument(argc, argv, opt) == 1) {
        Usage(opt);
        exit(EXIT_FAILURE);
    }
    Bam2Bed b(opt);
    b.Run();
    return 0;
}