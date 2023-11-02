#pragma once

#include "file_io.hpp"

class SeqReader {
public:
    using ID = uint32_t;
    struct Unit {
        std::string head;
        std::string sub_head;
        std::string seq;
        std::string quality;
        ID id;

        Unit() = default;
        ~Unit() = default;
    };

    SeqReader() = default;
    virtual ~SeqReader() = default;

    virtual bool Valid() const = 0;
    virtual bool IsEnd() const = 0;

    virtual bool Next(Unit &u) = 0;
    virtual std::size_t NextBatch(std::vector<Unit> &u, std::size_t size, int min_length = 0) = 0;
};

class FastaReader: public SeqReader {
public:
    FastaReader() = default;
    FastaReader(const std::string &fname): in_(fname) {}

    virtual ~FastaReader() = default;

    bool Valid() const { return in_.Valid(); }
    bool IsEnd() const { return in_.IsEnd(); }

    virtual bool Next(Unit &u);
    virtual std::size_t NextBatch(std::vector<Unit> &u, std::size_t size, int min_length = 0);

private:
    bool GetHead(std::string &head, std::string &sub_head);
    bool GetSeq(std::string &seq);

    GzFileReader in_;
};

class FastqReader: public SeqReader {
public:
    FastqReader() = default;
    FastqReader(const std::string& fname): in_(fname) {}
    virtual ~FastqReader() = default;

    bool Valid() const { return in_.Valid(); }
    bool IsEnd() const { return in_.IsEnd(); }

    virtual bool Next(Unit &u);
    virtual std::size_t NextBatch(std::vector<Unit> &u, std::size_t size, int min_length = 0);

private:
    bool GetHead(std::string &head, std::string &sub_head);
    bool GetSeq(std::string &seq) { seq = in_.GetStrippedLine(); return true; }
    bool GetComment() { std::string line = in_.GetStrippedLine(); return line[0] == '+'; }
    bool GetQuality(std::string &quality) { quality = in_.GetStrippedLine(); return true; }

    GzFileReader in_;
};

class SeqWriter {
public:
    SeqWriter() = default;
    virtual ~SeqWriter() = default;
    
    virtual void Write(const SeqReader::Unit &u) = 0;
    virtual void WriteFormat(const SeqReader::Unit &u, int num = 80) = 0;
};

class FastaWriter: public SeqWriter {
public:
    FastaWriter() = default;
    FastaWriter(const std::string &fanme): out_(fanme) {}
    virtual ~FastaWriter() = default;

    void Write(const SeqReader::Unit &u) {
        out_.Write(">");
        out_.Write(u.head);
        if(!u.sub_head.empty()) {
            out_.Write(" ");
            out_.Write(u.sub_head);
        }
        out_.Write("\n");
        out_.Write(u.seq + "\n");
    }

    void WriteFormat(const SeqReader::Unit &u, int num = 80) {
        out_.Write(">");
        out_.Write(u.head);
        if(!u.sub_head.empty()) {
            out_.Write(" ");
            out_.Write(u.sub_head);
        }
        out_.Write("\n");
        for(std::size_t i = 0; i < u.seq.size(); i += num) {
            if(i + num > u.seq.size()) num = u.seq.size() - i;
            out_.Write(u.seq.substr(i, num) + "\n");
        }
    }

private:
    GzFileWriter out_;
};