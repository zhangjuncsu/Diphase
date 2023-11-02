#pragma once

#include <time.h>
#include <string>
#include <vector>
#include <thread>
#include <algorithm>

template<typename F>
void MultiThreads(int threads, F func) {
    std::vector<std::thread> worker;
    for(int i = 0; i < threads; ++i) worker.emplace_back(func, i);
    for(std::size_t i = 0; i < worker.size(); ++i) worker[i].join();
}

inline std::vector<std::string> SplitString(const std::string &str, char delim) {
    auto is_split_point = [delim](char c) { return c == delim; };
    auto is_not_split_point = [delim](char c) { return c != delim; };
    
    std::vector<std::string> substrs;
    auto begin = std::find_if(str.begin(), str.end(), is_not_split_point);
    while(begin != str.end()) {
        auto end = std::find_if(begin, str.end(), is_split_point);
        substrs.emplace_back(begin, end);
        begin = std::find_if(end, str.end(), is_not_split_point);
    }
    return substrs;
}

inline std::string GetCurTime() {
    time_t now;
    time(&now);
    char tmp[64];
    strftime(tmp, sizeof(tmp), "%Y-%m-%d %H:%M:%S", localtime(&now));
    return std::string(tmp);
}