#pragma once

#include "seq_io.hpp"
#include "sequence.hpp"

#include <unordered_map>

class ReadStore {
public:
    ReadStore() = default;
    struct Unit {
        Unit() = default;
        Unit(const std::string &seq, SeqReader::ID id): seq(seq), id(id) {}

        DNASeq seq;
        SeqReader::ID id;
    };

    void Load(const std::string &fname, const std::string &type = "", int min_length = 0);
    void LoadFasta(const std::string &fname, int min_length = 0);
    void LoadFastq(const std::string &fname, int min_length = 0);

    SeqReader::ID GetIdByName(const std::string &name) const;
    std::string GetNameById(SeqReader::ID id) const { return names_[id]; }
    DNASeq GetSeq(SeqReader::ID id) const { return items_[id].seq; }
    DNASeq GetSeq(const std::string &name) const { return GetSeq(GetIdByName(name)); }
    DNASeq GetSeq(const std::string &name) { return GetSeq(GetIdByName(name)); }
    std::size_t GetSeqLength(SeqReader::ID id) const { return GetSeq(id).Size(); }
    std::size_t GetSeqLength(const std::string &name) const { return GetSeq(name).Size(); }

    Unit& GetUnit(SeqReader::ID id) { assert(id < items_.size()); return items_[id]; }
    std::vector<Unit> &GetUnits() { return items_; }
    // std::vector<Unit> &GetUnits() const { return items_; }

    void SaveId2Name(const std::string &fname) const;
    std::size_t Size() const { return names_.size(); }

    static SeqReader *OpenFile(const std::string &fname);

private:
    std::vector<Unit> items_;
    std::vector<std::string> names_;
    std::unordered_map<std::string, SeqReader::ID> name2id_;
};