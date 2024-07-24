#pragma once

#include <time.h>
#include <string>
#include <vector>
#include <thread>
#include <algorithm>

#include <sys/time.h>
#include <sys/resource.h>

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

static double real_time_start;

double GetCPUTime() {
    struct rusage rusage;
    getrusage(RUSAGE_SELF, &rusage);
    return rusage.ru_utime.tv_sec + rusage.ru_stime.tv_sec + 1e-6 * (rusage.ru_utime.tv_usec + rusage.ru_stime.tv_usec);
}

static inline double RealTime() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + 1e-6 * tv.tv_usec;
}

void ResetRealTime() {
    real_time_start = RealTime();
}

double GetRealTime() {
    return RealTime() - real_time_start;
}

long GetPeakMemory() {
    struct rusage rusage;
    getrusage(RUSAGE_SELF, &rusage);
    return rusage.ru_maxrss;
}

double PeakMemoryGB() {
    return GetPeakMemory() / 1048576.0;
}

double GetCpuUsage() {
    return (GetCPUTime() + 1e-9) / (GetRealTime() + 1e-9);
}