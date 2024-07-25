#pragma once

#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include <zlib.h>

class Reader {
public:
    Reader() = default;
    virtual ~Reader() = default;
    virtual bool Valid() const = 0;
    virtual bool IsEnd() const = 0;
    virtual std::size_t GetLines(std::vector<std::string> &lines) = 0;
};

class StdioReader: public Reader {
public:
    StdioReader() { std::ios::sync_with_stdio(false); }
    virtual ~StdioReader() = default;

    bool Valid() const { return true; }
    bool IsEnd() const { return false; }

    std::size_t GetLines(std::vector<std::string> &lines);
};

class GzFileReader: public Reader {
public:
    GzFileReader() = default;
    GzFileReader(const std::string &fname) { in_ = gzopen(fname.c_str(), "rb"); }
    ~GzFileReader() { if(in_ != nullptr) gzclose(in_); }

    bool Valid() const { return in_ != nullptr; }
    bool IsEnd() const { return gzeof(in_); }

    void ConsumeStrippedLine() { next_line_pos_ = -1; }
    bool HasNoEmptyLine() const { return next_line_pos_ >= 0; }
    int64_t Tell() { return !HasNoEmptyLine()? gztell(in_): next_line_pos_; }

    std::string GetStrippedLine();
    std::string GetNoEmptyLine();

    std::string GetLine();
    bool GetLine(std::string &line);
    std::size_t GetLines(std::vector<std::string> &lines);

    const std::string &QueryNoEmptyLine() {
        if(!HasNoEmptyLine()) {
            int64_t pos = Tell();
            next_line_ = GetNoEmptyLine();
            next_line_pos_ = pos;
        }
        return next_line_;
    }

    void Seek(int64_t pos) {
        next_line_.clear();
        next_line_pos_ = pos;
        gzseek(in_, pos, SEEK_SET);
    }

private:
    gzFile in_ { nullptr };
    std::string next_line_ { "" };
    int64_t next_line_pos_ { -1 };
};

class Writer {
public:
    Writer() = default;
    virtual ~Writer() = default;
    virtual void Write(const std::string &str) = 0;
    virtual void Flush() = 0;
};

inline Writer &operator<<(Writer &writer, const std::string &str) {
    writer.Write(str);
    return writer;
}

template <typename T>
inline Writer &operator<<(Writer &writer, T &t) {
    std::ostringstream oss;
    oss << t;
    writer << oss.str();
    return writer;
}

class StdioWriter: public Writer {
public:
    StdioWriter() = default;
    virtual ~StdioWriter() = default;

    virtual void Write(const std::string &str) { std::cout << str; };
    virtual void Flush() {}
};

class GzFileWriter: public Writer {
public:
    GzFileWriter() = default;
    GzFileWriter(const std::string &fname) {
        compress_ = fname.substr(fname.find_last_of(".")) == ".gz";
        if(compress_) {
            out_compress_ = gzopen(fname.c_str(), "wb");
        } else {
            out_decompress_ = fopen(fname.c_str(), "w");
        }
    }

    ~GzFileWriter() {
        if(out_compress_ != nullptr) gzclose(out_compress_);
        if(out_decompress_ != nullptr) fclose(out_decompress_);
    }

    void Write(const std::string &str) {
        if(compress_) {
            gzputs(out_compress_, str.c_str());
        } else {
            fputs(str.c_str(), out_decompress_);
        }
    }

    void Flush() {}

private:
    bool compress_;
    gzFile out_compress_ { nullptr };
    FILE *out_decompress_ { nullptr };
};

std::vector<std::string> GetLineFromFile(const std::string &fname);