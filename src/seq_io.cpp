#include "seq_io.hpp"
#include "utility.hpp"
#include "logger.hpp"

#include <cassert>

bool FastaReader::Next(Unit &u) {
    return GetHead(u.head, u.sub_head) && GetSeq(u.seq);
}

std::size_t FastaReader::NextBatch(std::vector<Unit>& units, std::size_t size, int min_length) {
    units.clear();
    Unit u;
    std::size_t sz = 0;
    while(sz < size && Next(u)) {
        if(u.seq.size() < min_length) continue;
        units.emplace_back(u);
        sz += u.seq.size();
    }
    return sz;
}

bool FastaReader::GetHead(std::string &head, std::string &sub_head) {
    std::string line = in_.QueryNoEmptyLine();
    in_.ConsumeStrippedLine();

    if(!line.empty() && line[0] == '>') {
        std::string::size_type s = 1;
        while(s < line.size() && ::isspace(line[s])) s++;

        std::string::size_type e = std::min(s, line.size());
        while(e < line.size() && !::isspace(line[e])) e++;

        if(e > s) {
            head = line.substr(s, e - s);
            while(e < line.size() && ::isspace(line[e])) e++;
            sub_head = line.substr(e);
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}

bool FastaReader::GetSeq(std::string &seq) {
    seq = "";
    std::string line = in_.QueryNoEmptyLine();
    while(!line.empty() && line[0] != '>') {
        seq += line;
        in_.ConsumeStrippedLine();
        line = in_.QueryNoEmptyLine();
    }
    return true;
}

bool FastqReader::Next(Unit& u) {
    assert(Valid());
    bool ret = GetHead(u.head, u.sub_head) && GetSeq(u.seq) && GetComment() && GetQuality(u.quality);
    if(u.seq.size() != u.quality.size()) {
        LOG(ERROR)("The length of sequence and quatity of read %s %s is not equal\n", u.head.c_str(), u.sub_head.c_str());
        // std::cerr << "[" << GetCurTime() << "] The length of sequence and quatity of read " << u.head << " " << u.sub_head << " is not equal\n";
        // exit(EXIT_FAILURE);
    }
    return ret;
}

std::size_t FastqReader::NextBatch(std::vector<Unit> &units, std::size_t size, int min_length) {
    units.clear();
    Unit u;
    std::size_t sz = 0;
    while(sz < size && Next(u)) {
        if(u.seq.size() < min_length) continue;
        units.emplace_back(u);
        sz += u.seq.size();
    }
    return sz;
}

bool FastqReader::GetHead(std::string &head, std::string &sub_head) {
    std::string line = in_.GetStrippedLine();
    if(!line.empty() && line[0] == '@') {
        std::string::size_type s = 1;
        while(s < line.size() && ::isspace(line[s])) s++;

        std::string::size_type e = std::min(s, line.size());
        while(e < line.size() && !::isspace(line[e])) e++;

        if(e > s) {
            head = line.substr(s, e - s);
            while(e < line.size() && ::isspace(line[e])) e++;
            sub_head = line.substr(e);
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}