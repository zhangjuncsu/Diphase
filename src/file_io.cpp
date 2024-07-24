#include "file_io.hpp"
#include "utility.hpp"
#include "logger.hpp"

std::string StrippedString(const std::string &str) {
    std::string::size_type s = 0;
    while(s < str.size() && ::isspace(str[s])) s++;

    std::string::size_type e = str.size();
    while(e > s && ::isspace(str[e - 1])) e--;
    
    return str.substr(s, e - s);
}

std::size_t StdioReader::GetLines(std::vector<std::string> &lines) {
    std::size_t sz = 0;
    for(; sz < lines.size(); ++sz) {
        bool ret = (bool)std::getline(std::cin, lines[sz]);
        if(!ret) break;
    }
    return sz;
}

std::string GzFileReader::GetLine() {
    std::string str;
    const int BUF_SIZE = 1024 * 2;
    char buf[BUF_SIZE];

    char *line = gzgets(in_, buf, BUF_SIZE);
    while(line != nullptr) {
        str += line;
        if(str.back() == '\n') break;
        line = gzgets(in_, buf, BUF_SIZE);
    }
    return str;
}

bool GzFileReader::GetLine(std::string &line) {
    line.clear();
    line = GetLine();
    return !line.empty();
}

std::size_t GzFileReader::GetLines(std::vector<std::string> &lines) {
    std::size_t sz = 0;
    while(sz < lines.size()) {
        lines[sz] = GetLine();
        if(lines[sz].empty()) break;
        sz++;
    }
    return sz;
}

std::string GzFileReader::GetStrippedLine() {
    return StrippedString(GetLine());
}

std::string GzFileReader::GetNoEmptyLine() {
    if(!HasNoEmptyLine()) {
        std::string line = GetStrippedLine();
        while(line.empty() && !IsEnd()) {
            line = GetStrippedLine();
        }
        return line;
    } else {
        next_line_pos_ = -1;
        return next_line_;
    }
}

std::vector<std::string> &GetLineFromFile(const std::string &fname) {
    GzFileReader in(fname);
    std::vector<std::string> lines;
    if(in.Valid()) {
        std::string str = in.GetNoEmptyLine();
        while(!str.empty()) {
            lines.emplace_back(str);
            str = in.GetNoEmptyLine();
        }
    } else {
        LOG(ERROR)("Failed to open file %s for reading\n", fname.c_str());
        // std::cerr << "[" << GetCurTime() << "] Failed to open file " << fname << " for reading\n";
        // exit(EXIT_FAILURE);
    }
    return lines;
}