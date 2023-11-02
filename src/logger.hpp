#pragma once

// logger class for logging information to console and file simultaneously 
// if file is specified by user

#include <mutex>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include <stdio.h>
#include <stdarg.h>

class Logger {
public:
    enum LEVEL {
        DEBUG = 0,
        INFO,
        WARNING,
        ERROR,
        FATAL
    };

    Logger();
    ~Logger();

    void SetLogLevel(LEVEL level) { level_ = level; }
    void SetLogFile(const std::string &fname);

    // template<typename... Args>
    void Log(LEVEL level, const char* const format, va_list args) {
        if(level < level_) return;
        std::lock_guard<std::mutex> lock(mutex_);
        std::stringstream ss;
        ss << "[" << GetTimeStamps() << "] [" << GetLevelStr(level) << "] ==> ";
        char buffer[1024];
        vsnprintf(buffer, 1024, format, args);
        ss << buffer;
        to_console_ << ss.str();
        to_console_.flush();
        if(to_file_.is_open()) {
            to_file_ << ss.str();
            to_file_.flush();
        }
        if(level >= ERROR) exit(EXIT_FAILURE);
    }

private:
    std::string GetTimeStamps();
    std::string GetLevelStr(LEVEL level);

    std::string log_file_;
    LEVEL level_;
    std::mutex mutex_;
    std::ofstream to_file_ { nullptr };
    std::ofstream to_console_;
};

class LoggerWrapper {
public:
    LoggerWrapper(Logger &logger, Logger::LEVEL level): logger_(logger), level_(level) {}
    ~LoggerWrapper() {}

    // template<typename... Args>
    void operator()(const char* const format, ...) {
        va_list args;
        va_start(args, format);
        logger_.Log(level_, format, args);
        va_end(args);
    }

private:
    Logger &logger_;
    Logger::LEVEL level_;
};

#define LOG(s) LoggerWrapper(LOGGER, Logger::s)

extern Logger LOGGER;