#include "logger.hpp"

#include <iomanip>

#include <time.h>

Logger LOGGER;

Logger::Logger() {
    to_console_ = std::ofstream("/dev/stderr");
    if(!to_console_.is_open()) {
        throw std::runtime_error("Failed to open /dev/stderr");
    }
    if(!to_console_.good()) {
        throw std::runtime_error("Failed to write to /dev/stderr");
    }

    if(!log_file_.empty()) {
        to_file_.open(log_file_);
        if(!to_file_.is_open()) {
            throw std::runtime_error("Failed to open " + log_file_);
        }
        if(!to_file_.good()) {
            throw std::runtime_error("Failed to write to " + log_file_);
        }
    }
}

Logger::~Logger() {
    to_console_.close();
    if(to_file_.is_open()) to_file_.close();
}

void Logger::SetLogFile(const std::string &fname) {
    log_file_ = fname;
    if(to_file_.is_open()) to_file_.close();
    to_file_.open(log_file_);
    if(!to_file_.is_open()) {
        throw std::runtime_error("Failed to open " + log_file_);
    }
    if(!to_file_.good()) {
        throw std::runtime_error("Failed to write to " + log_file_);
    }
}

std::string Logger::GetTimeStamps() {
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
    return oss.str();
}

std::string Logger::GetLevelStr(LEVEL level) {
    switch(level) {
        case DEBUG: return "DEBUG";
        case INFO: return "INFO";
        case WARNING: return "WARNING";
        case ERROR: return "ERROR";
        case FATAL: return "FATAL";
        default: return "INFO";
    }
}