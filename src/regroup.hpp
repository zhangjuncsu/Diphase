#pragma once

#include "bam_reader.hpp"
#include "logger.hpp"
#include "matrix.hpp"

#include <unordered_map>

struct Option {
    std::string bfname;
    std::string ofname;
    std::string b1fname;
    std::string b2fname;
    std::string pfname;

    int threads { 1 };
    int iter { 10000 };
    int seed { 1000 };

    bool more { false };
};

class Regroup {
public:
    Regroup(Option& opt): opt_(opt) {}

    struct Group {
        std::size_t first, second;
        int label_first, label_second;
    };

    void Count(const std::string &b1fname, const std::string &b2fname);
    void LoadPair();
    std::vector<std::string> LoadBed(const std::string &fname);
    void LoadGroup();
    void GenerateMatrix();
    void GenerateMatrixMore();
    int GetMax() const;
    void LocalMaximum();
    void Phase();
    void Run();

private:
    Option& opt_;

    Matrix matrix_;
    std::vector<Group> groups_;
    std::vector<std::string> hap1_;
    std::vector<std::string> hap2_;
    std::vector<std::string> names_;
    std::vector<std::string> pairs_;
    std::unordered_map<std::string, int> count_;
    std::unordered_map<std::string, std::size_t> index_;
    std::unordered_map<std::string, std::string> mapped_;
};