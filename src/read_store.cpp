#include "read_store.hpp"
#include "utility.hpp"

SeqReader *ReadStore::OpenFile(const std::string &fname) {
    SeqReader *reader;
    std::size_t last_dot_pos = fname.rfind('.');
    std::string fname_nogz = fname.substr(last_dot_pos) == ".gz"? fname.substr(0, last_dot_pos): fname;
    std::string type = fname_nogz.substr(fname_nogz.rfind('.') + 1);
    if(type == "fasta" || type == "fa") {
        reader = new FastaReader(fname);
    } else if(type == "fastq" || type == "fq") {
        reader = new FastqReader(fname);
    } else {
        std::cerr << "[" << GetCurTime() << "] Unrecognized file type " << type << "\n";
        exit(EXIT_FAILURE);
    }
    if(reader->Valid()) {
        std::cerr << "[" << GetCurTime() << "] Failed to open " << fname << " for reading\n";
        exit(EXIT_FAILURE);
    }
    return reader;
}

void ReadStore::Load(const std::string &fname, const std::string &type, int min_length) {
    std::string t = type;
    if(t == "") {
        std::string fname_nogz = fname;
        std::size_t last_dot_pos = fname_nogz.rfind('.');
        if(fname_nogz.substr(last_dot_pos) == ".gz") {
            fname_nogz = fname_nogz.substr(0, last_dot_pos);
        }
        t = fname_nogz.substr(fname_nogz.rfind('.') + 1);
    }
    if(t == "fasta" || t == "fa") {
        LoadFasta(fname, min_length);
    } else if(t == "fastq" || t == "fq") {
        LoadFastq(fname, min_length);
    } else {
        std::cerr << "[" << GetCurTime() << "] Unrecognized file type " << fname << "\n";
        exit(EXIT_FAILURE);
    }
}

void ReadStore::LoadFasta(const std::string &fname, int min_length) {
    FastaReader reader(fname);
    SeqReader::Unit u;
    if(reader.Valid()) {
        while(reader.Next(u)) {
            assert(!u.head.empty());
            if(u.seq.size() < min_length) continue;
            SeqReader::ID id = names_.size();
            if(names_.capacity() == names_.size()) {
                names_.reserve(names_.capacity() * 1.5);
                name2id_.reserve(names_.capacity());
            }
            names_.emplace_back(u.head);
            name2id_[u.head] = id;
            items_.emplace_back(u.seq, id);
        }
        if(!reader.IsEnd()) {
            std::cerr << "[" << GetCurTime() << "] Not all reads in file " << fname << " are loaded\n";
        }
    } else {
        std::cerr << "[" << GetCurTime() << "] Failed to open file " << fname << "\n";
        exit(EXIT_FAILURE);
    }
    std::cerr << "[" << GetCurTime() << "] Load " << names_.size() << " reads from file " << fname << "\n";
}

void ReadStore::LoadFastq(const std::string &fname, int min_length) {
    FastqReader reader(fname);
    SeqReader::Unit u;
    if(reader.Valid()) {
        while(reader.Next(u)) {
            assert(!u.head.empty());
            if(u.seq.size() < min_length) continue;
            SeqReader::ID id = names_.size();
            if(names_.capacity() == names_.size()) {
                names_.reserve(names_.capacity() * 1.5);
                name2id_.reserve(names_.capacity());
            }
            names_.emplace_back(u.head);
            name2id_[u.head] = id;
            items_.emplace_back(u.head, id);
        }
        if(!reader.IsEnd()) {
            std::cerr << "[" << GetCurTime() << "] Not all reads in file " << fname << " are loaded\n";
        }
    } else {
        std::cerr << "[" << GetCurTime() << "] Failed to open file " << fname << "\n";
        exit(EXIT_FAILURE);
    }
    std::cerr << "[" << GetCurTime() << "] Load " << names_.size() << " reads in file " << fname << "\n";
}

SeqReader::ID ReadStore::GetIdByName(const std::string &name) const {
    auto iter = name2id_.find(name);
    if(iter == name2id_.end()) {
        std::cerr << "[" << GetCurTime() << "] Name " << name << " not found\n";
        return -1;
    } else {
        return iter->second;
    }
}

void ReadStore::SaveId2Name(const std::string &fname) const {
    GzFileWriter writer(fname);
    for(std::size_t i = 0; i < names_.size(); ++i) {
        writer.Write(std::to_string(i) + ' ');
        writer.Write(names_[i] + '\n');
    }
}