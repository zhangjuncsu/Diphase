#include "hic.hpp"
#include "logger.hpp"
#include "utility.hpp"
// #include "thread_pool.hpp"

#include <regex>
#include <mutex>
#include <atomic>
#include <bitset>

#include "getopt.h"

bool debug = false;

std::size_t GetMemoryUsage() {
    int vmrss = 0;
    int pid = getpid();
    char fname[64];
    sprintf(fname, "/proc/%d/status", pid);
    FILE *fp = fopen(fname, "r");
    if(fp == NULL) return 0;
    const int BUFSIZE = 1024;
    char buf[BUFSIZE];
    while(fgets(buf, BUFSIZE - 1, fp)) {
        if(strncmp(buf, "VmRSS:", 6) == 0) {
            sscanf(buf + 6, "%d", &vmrss);
            break;
        }
    }
    fclose(fp);
    return vmrss;
}

int Option::Check() const {
    if(vfname.empty()) {
        LOG(WARNING)("No variant file is specified");
        return 1;
    }
    if(rfname.empty()) {
        LOG(WARNING)("No reference file is specified");
        return 1;
    }
    // if(h1fname.empty()) {
    //     LOG(WARNING)("No Hi-C mate-pair 1 file is specified");
    //     return 1;
    // }
    // if(h2fname.empty()) {
    //     LOG(WARNING)("No Hi-C mate-pair 2 file is specified");
    //     return 1;
    // }
    if(k < 1) {
        LOG(WARNING)("Invalid k-mer size: %d", k);
        return 1;
    }
    if(threads < 1) {
        LOG(WARNING)("Invalid number of threads: %d", threads);
        return 1;
    }
    if(cnt_thres < 1) {
        LOG(WARNING)("Invalid threshold of k-mer count: %d", cnt_thres);
        return 1;
    }
    return 0;
}

inline uint64_t Hash64_64(uint64_t key) {
    key = ~key + (key << 21);
    key = key ^ key >> 24;
    key = (key + (key << 3)) + (key << 8);
    key = key ^ key >> 14;
    key = (key + (key << 2)) + (key << 4);
    key = key ^ key >> 28;
    key = key + (key << 31);
    return key;
}

inline uint64_t GetKmerDirection(std::array<uint64_t, 4> &kmer) {
    if(kmer[1] != kmer[3]) return kmer[1] < kmer[3]? 0: 1;
    else if (kmer[0] != kmer[2]) return kmer[0] < kmer[2]? 0: 1;
    else return (uint64_t)-1;
}

inline uint64_t HashLongKmer(std::array<uint64_t, 4> &kmer, uint64_t &dir, uint64_t k) {
    dir = GetKmerDirection(kmer);
    if(dir == (uint64_t)-1) return dir;
    if(k <= 32) return ((kmer[dir << 1 | 0] << 32) | (kmer[dir << 1 | 1]));
    return Hash64_64(kmer[dir << 1 | 0]) + Hash64_64(kmer[dir << 1 | 1]);
}

void HiC::LoadVariants() {
    for(auto vfname: opt_.vfname){
        GzFileReader reader(vfname);
        if(reader.Valid()) {
            std::string line = reader.GetNoEmptyLine();
            while(!line.empty()) {
                if(line[0] != '#') {
                    auto items = SplitString(line, '\t');
                    if(items[6][1] == 'L') continue; // filter low quality variants
                    // if(items[3].size() != 1 && items[4].size() != 1) continue;
                    long pos = atol(items[1].c_str()) - 1;
                    for(auto i = 0; i < items[3].size(); ++i) {
                        variants_[items[0]].emplace_back(pos + i);
                    }
                }
                line = reader.GetNoEmptyLine();
            }
        } else {
            LOG(ERROR)("Failed to open %s for reading", vfname.c_str());
        }
    }
}

template<typename RandomAccessIterator, typename KeySort>
void HiC::RadixSort(RandomAccessIterator first, RandomAccessIterator last, uint8_t max_bit, KeySort ks) {
    if(first >= last) return;
    std::vector<typename std::iterator_traits<RandomAccessIterator>::value_type> tmp(first, last);
    auto begin = tmp.begin();
    auto end = tmp.end();

    uint64_t buckets[0x100]{};
    uint8_t shift = 0;
    for(; shift < max_bit; shift += 8) {
        uint64_t count[0x100]{};
        for(auto iter = first; iter != last; ++iter) {
            ++count[ks(*iter) >> shift & 0xff];
        }
        for(uint64_t i = 0, j = 0; i < 0x100; j += count[i++]) {
            buckets[i] = j;
        }
        for(auto iter = first; iter != last; ++iter) {
            *(begin + buckets[ks(*iter) >> shift & 0xff]++) = *iter;
        }
        std::swap(begin, first);
        std::swap(end, last);
    }
}

uint32_t HiC::Index::Find(uint64_t key, const uint64_t **dst) const {
    auto iter = index_.find(key << 1);
    if(iter == index_.end()) return 0;
    if(iter->first & 1) {
        *dst = &iter->second;
        return 1;
    }
    *dst = &offset_[iter->second >> 32];
    return static_cast<uint32_t>(iter->second);
}

uint64_t HiC::Hash(uint64_t key, uint64_t mask) {
    key = ((~key) + (key << 21)) & mask;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask;
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask;
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}

std::vector<HiC::Kmer> HiC::Sketch(const ReadStore::Unit &unit) {
    if(unit.seq.Size() < opt_.k) return std::vector<Kmer>();
    uint64_t mask = (1ULL << opt_.k) - 1;

    std::regex pattern("(\\w+):(\\d+)-(\\d+)");
    std::smatch result;

    uint64_t shift = opt_.k - 1;
    
    uint64_t id = static_cast<uint64_t>(unit.id) << 32;
    std::vector<Kmer> kmers;
    std::string name = rs_ref_.GetNameById(unit.id);
    int offset = 0;
    if(std::regex_search(name, result, pattern)) {
        name = result[1].str();
        offset = atoi(result[2].str().c_str());
    }
    if(variants_.find(name) == variants_.end()) return std::vector<Kmer>();
    auto vars = variants_[name];
    std::size_t prev = 0;
    for(auto pos: vars) {
        uint32_t l = 0;
        std::array<uint64_t, 4> minimer = {0, 0, 0, 0};
        uint64_t dir;
        if(pos < offset) continue;
        if(pos - opt_.k + 1 > offset + unit.seq.Size()) break;
        for(auto i = pos - offset - opt_.k + 1; i < pos - offset + opt_.k; ++i) {
            if(i < prev) continue;
            if(i >= unit.seq.Size()) break;
            uint8_t c = unit.seq[i];
            if(c < 4) {
                // minimer[0] & minimer[1] are forward strand
                // minimer[2] & minimer[3] are reverse strand
                minimer[0] = (minimer[0] << 1 | (c & 1)) & mask;
                minimer[1] = (minimer[1] << 1 | (c >> 1)) & mask;
                minimer[2] = minimer[2] >> 1 | (uint64_t)(1 - (c & 1)) << shift;
                minimer[3] = minimer[3] >> 1 | (uint64_t)(1 - (c >> 1)) << shift;
                if(++l >= opt_.k) {
                    uint64_t hash = HashLongKmer(minimer, dir, opt_.k);
                    if(dir == (uint64_t)-1) continue;
                    kmers.emplace_back(hash, id | (i - (opt_.k - 1U)) << 1 | dir);
                }
            } else {
                l = 0;
                minimer[0] = minimer[1] = minimer[2] = minimer[3] = 0;
            }
        }
        prev = pos - offset + opt_.k;
    }
    return kmers;
}

std::vector<HiC::Kmer> HiC::GetKmer(const ReadStore::Unit &unit) {
    if(unit.seq.Size() < opt_.k) return std::vector<Kmer>();
    uint64_t mask = (1ULL << opt_.k) - 1;

    uint64_t shift = opt_.k - 1;
    std::array<uint64_t, 4> minimer = {0, 0, 0, 0};
    uint32_t l = 0;
    uint64_t id = static_cast<uint64_t>(unit.id) << 32;
    std::vector<Kmer> kmers;
    for(auto i = 0; i < unit.seq.Size(); ++i) {
        uint8_t c = unit.seq[i];
        uint64_t dir;
        if(c < 4) {
            // minimer[0] & minimer[1] are forward strand
            // minimer[2] & minimer[3] are reverse strand
            minimer[0] = (minimer[0] << 1 | (c & 1)) & mask;
            minimer[1] = (minimer[1] << 1 | (c >> 1)) & mask;
            minimer[2] = minimer[2] >> 1 | (uint64_t)(1 - (c & 1)) << shift;
            minimer[3] = minimer[3] >> 1 | (uint64_t)(1 - (c >> 1)) << shift;
            if(++l >= opt_.k) {
                uint64_t hash = HashLongKmer(minimer, dir, opt_.k);
                if(dir == (uint64_t)-1) continue;
                kmers.emplace_back(hash, id | (i - (opt_.k - 1U)) << 1 | dir);
             }
        } else {
            l = 0;
            minimer[0] = minimer[1] = minimer[2] = minimer[3] = 0;
        }
    }
    return kmers;
}

// std::vector<HiC::Kmer> HiC::Sketch(const ReadStore::Unit &unit) {
//     if(unit.seq.Size() < opt_.k) return std::vector<Kmer>();
//     uint64_t mask = (1ULL << (opt_.k * 2)) - 1;

//     std::regex pattern("(\\w+):(\\d+)-(\\d+)");
//     std::smatch result;

//     uint64_t shift = (opt_.k - 1) * 2;
    
//     uint64_t id = static_cast<uint64_t>(unit.id) << 32;
//     std::vector<Kmer> kmers;
//     std::string name = rs_ref_.GetNameById(unit.id);
//     int offset = 0;
//     if(std::regex_search(name, result, pattern)) {
//         name = result[1].str();
//         offset = atoi(result[2].str().c_str());
//     }
//     if(variants_.find(name) == variants_.end()) return std::vector<Kmer>();
//     auto vars = variants_[name];
//     std::size_t prev = 0;
//     for(auto pos: vars) {
//         uint32_t l = 0;
//         uint64_t minimer[2] = {0, 0};
//         if(pos < offset) continue;
//         if(pos - opt_.k + 1 > offset + unit.seq.Size()) break;
//         for(auto i = pos - offset - opt_.k + 1; i < pos - offset + opt_.k; ++i) {
//             if(i < prev) continue;
//             if(i >= unit.seq.Size()) break;
//             uint8_t c = unit.seq[i];
//             if(c < 4) {
//                 minimer[0] = (minimer[0] << 2 | c) & mask;
//                 minimer[1] = (minimer[1] >> 2) | (3ULL ^ c) << shift;
//                 if(++l >= opt_.k) {
//                     int z = minimer[0] < minimer[1]? 0: 1;
//                     kmers.emplace_back(Hash(minimer[z], mask), id | (i - (opt_.k - 1U)) << 1 | z);
//                     // if(Hash(minimer[z], mask) == 3058399470777975378 || Hash(minimer[z], mask) == 731064931610552558 || Hash(minimer[z], mask) == 4447043543039867056) {
//                     //     std::cerr << "kmer: " << (*unit.seq.ToString()).substr(i - opt_.k + 1, opt_.k) << " hash: " << Hash(minimer[z], mask) << " id: " << (id >> 32) << " name: " << rs_ref_.GetNameById((id >> 32)) << " pos: " << i - opt_.k + 1 << " strand: " << z << " var: " << pos << std::endl;
//                     // }
//                 }
//             } else {
//                 l = 0;
//                 minimer[0] = minimer[1] = 0;
//             }
//         }
//         prev = pos - offset + opt_.k;
//     }
//     return kmers;
// }

// std::vector<HiC::Kmer> HiC::GetKmer(const ReadStore::Unit &unit) {
//     if(unit.seq.Size() < opt_.k) return std::vector<Kmer>();
//     uint64_t mask = (1ULL << (opt_.k * 2)) - 1;

//     uint64_t shift = (opt_.k - 1) * 2;
//     uint64_t minimer[2] = {0, 0};
//     uint32_t l = 0;
//     uint64_t id = static_cast<uint64_t>(unit.id) << 32;
//     std::vector<Kmer> kmers;
//     for(auto i = 0; i < unit.seq.Size(); ++i) {
//         uint8_t c = unit.seq[i];
//         if(c < 4) {
//             minimer[0] = (minimer[0] << 2 | c) & mask;
//             minimer[1] = (minimer[1] >> 2) | (3ULL ^ c) << shift;
//             if(++l >= opt_.k) {
//                 int z = minimer[0] < minimer[1]? 0: 1;
//                 kmers.emplace_back(Hash(minimer[z], mask), id | (i - (opt_.k - 1U)) << 1 | z);
//              }
//         } else {
//             l = 0;
//             minimer[0] = minimer[1] = 0;
//         }
//     }
//     return kmers;
// }

void HiC::BuildIndexAll() {
    for(auto &iter: index_all_) {
        iter.index_.clear();
        iter.offset_.clear();
    }

    std::size_t i = 0;
    std::vector<std::vector<Kmer>> minimizers(index_all_.size());
    while(i < rs_ref_.Size()) {
        // LOG(INFO)("vmRSS: %zd\n", GetMemoryUsage());
        std::size_t j = i;
        std::size_t byte = 0;
        while(j < rs_ref_.Size() && byte < 50000000 /*(1ULL << 32)*/) {
            byte += rs_ref_.GetUnit(j).seq.Size();
            ++j;
        }
        uint64_t mask = index_all_.size() - 1;
        std::atomic<std::size_t> index { i };
        std::vector<std::vector<Kmer>> tmp(j - i);
        auto func_sketch = [&](std::size_t tid) {
            for(auto k = index.fetch_add(1); k < j; k = index.fetch_add(1)) {
                tmp[k - i] = GetKmer(rs_ref_.GetUnit(k));
            }
        };
        MultiThreads(opt_.threads, func_sketch);
        for(auto &iter: tmp) {
            for(auto &kmer: iter) {
                auto &minimizer = minimizers[kmer.hash & mask];
                if(minimizer.capacity() == minimizer.size()) minimizer.reserve(minimizer.capacity() * 1.5);
                minimizer.emplace_back(kmer);
            }
        }
        i = j;
        std::vector<std::vector<Kmer>>().swap(tmp);
    }
    LOG(INFO)("counting finished vmRSS: %zd\n", GetMemoryUsage());
    std::vector<std::pair<std::size_t, std::size_t>> pairs(minimizers.size());
    std::atomic<std::size_t> index { 0 };
    auto func_add = [&](std::size_t tid) {
        for(auto i = index.fetch_add(1); i < minimizers.size(); i = index.fetch_add(1)) {
            if(minimizers[i].empty()) {
                pairs[i] = std::make_pair(0, 0);
                continue;
            }
            RadixSort(minimizers[i].begin(), minimizers[i].end(), opt_.k * 2, Kmer::SortByHash);
            minimizers[i].emplace_back(-1, -1);
            std::size_t num_hashes = 0;
            std::size_t num_keys = 0;
            for(uint64_t j = 1, c = 1; j < minimizers[i].size(); ++j, ++c) {
                if(minimizers[i][j - 1].hash != minimizers[i][j].hash) {
                    if(c > 1) num_hashes += c;
                    ++num_keys;
                    c = 0;
                }
            }
            pairs[i] = std::make_pair(num_hashes, num_keys);
        }
    };
    MultiThreads(opt_.threads, func_add);
    // LOG(INFO)("vmRSS: %zd\n", GetMemoryUsage());
    for(std::size_t k = 0; k < minimizers.size(); ++k) {
        if(minimizers[k].empty()) continue;
        index_all_[k].index_.reserve(pairs[k].second);
        // index_all_[k].offset_.reserve(pairs[k].first);
        for(uint64_t l = 1, c = 1; l < minimizers[k].size(); ++l, ++c) {
            if(minimizers[k][l - 1].hash != minimizers[k][l].hash) {
                if(c == 1) {
                    index_all_[k].index_.emplace(minimizers[k][l - 1].hash << 1 | 1, minimizers[k][l - 1].pos);
                } else {
                    index_all_[k].index_.emplace(minimizers[k][l - 1].hash << 1, index_all_[k].offset_.size() << 32 | c);
                    // for(uint64_t m = l - c; m < l; ++m) {
                    //     index_all_[k].offset_.emplace_back(minimizers[k][m].pos);
                    // }
                }
                c = 0;
            }
         }
        std::vector<Kmer>().swap(minimizers[k]);
    }
    std::vector<std::pair<std::size_t, std::size_t>>().swap(pairs);
    std::vector<std::vector<Kmer>>().swap(minimizers);
}

void HiC::BuildIndexSingle() {
    for(auto &iter: index_) {
        iter.index_.clear();
        iter.offset_.clear();
    }

    std::size_t i = 0;
    std::vector<std::vector<Kmer>> minimizers(index_.size());
    while(i < rs_ref_.Size()) {
        std::size_t j = i;
        std::size_t byte = 0;
        while(j < rs_ref_.Size() && byte < 50000000 /*(1ULL << 32)*/) {
            byte += rs_ref_.GetUnit(j).seq.Size();
            ++j;
        }
        uint64_t mask = index_.size() - 1;
        std::atomic<std::size_t> index { i };
        std::vector<std::vector<Kmer>> tmp(j - i);
        auto func_sketch = [&](std::size_t tid) {
            for(auto k = index.fetch_add(1); k < j; k = index.fetch_add(1)) {
                tmp[k - i] = GetKmer(rs_ref_.GetUnit(k));
            }
        };
        MultiThreads(opt_.threads, func_sketch);
        for(auto &iter: tmp) {
            for(auto &kmer: iter) {
                auto &minimizer = minimizers[kmer.hash & mask];
                if(minimizer.capacity() == minimizer.size()) minimizer.reserve(minimizer.capacity() * 1.5);
                minimizer.emplace_back(kmer);
            }
        }
        i = j;
        std::vector<std::vector<Kmer>>().swap(tmp);
    }
    
    std::vector<std::pair<std::size_t, std::size_t>> pairs(minimizers.size());
    std::atomic<std::size_t> index { 0 };
    auto func_add = [&](std::size_t tid) {
        for(auto i = index.fetch_add(1); i < minimizers.size(); i = index.fetch_add(1)) {
            if(minimizers[i].empty()) {
                pairs[i] = std::make_pair(0, 0);
                continue;
            }
            RadixSort(minimizers[i].begin(), minimizers[i].end(), opt_.k * 2, Kmer::SortByHash);
            minimizers[i].emplace_back(-1, -1);
            std::size_t num_hashes = 0;
            std::size_t num_keys = 0;
            for(uint64_t j = 1, c = 1; j < minimizers[i].size(); ++j, ++c) {
                if(minimizers[i][j - 1].hash != minimizers[i][j].hash) {
                    if(c > 1) num_hashes += c;
                    if(c == 1) ++num_keys;
                    c = 0;
                }
            }
            pairs[i] = std::make_pair(num_hashes, num_keys);
        }
    };
    MultiThreads(opt_.threads, func_add);
    // LOG(INFO)("vmRSS: %zd\n", GetMemoryUsage());
    for(std::size_t k = 0; k < minimizers.size(); ++k) {
        if(minimizers[k].empty()) continue;
        index_[k].index_.reserve(pairs[k].second);
        for(uint64_t l = 1, c = 1; l < minimizers[k].size(); ++l, ++c) {
            if(minimizers[k][l - 1].hash != minimizers[k][l].hash) {
                if(c == 1) {
                    index_[k].index_.emplace(minimizers[k][l - 1].hash << 1 | 1, minimizers[k][l - 1].pos);
                }
                c = 0;
            }
         }
        std::vector<Kmer>().swap(minimizers[k]);
    }
    std::vector<std::pair<std::size_t, std::size_t>>().swap(pairs);
    std::vector<std::vector<Kmer>>().swap(minimizers);
}

void HiC::BuildIndex() {
    for(auto &iter: index_) {
        iter.index_.clear();
        iter.offset_.clear();
    }

    std::size_t i = 0;
    std::vector<std::vector<Kmer>> minimizers(index_.size());
    while(i < rs_ref_.Size()) {
        std::size_t j = i;
        std::size_t byte = 0;
        while(j < rs_ref_.Size() && byte < (1ULL << 32)) {
            byte += rs_ref_.GetUnit(j).seq.Size();
            ++j;
        }
        uint64_t mask = index_.size() - 1;
        std::atomic<std::size_t> index { i };
        std::vector<std::vector<Kmer>> tmp(j - i);
        auto func_sketch = [&](std::size_t tid) {
            for(auto k = index.fetch_add(1); k < j; k = index.fetch_add(1)) {
                tmp[k - i] = Sketch(rs_ref_.GetUnit(k));
                // if(rs_ref_.GetNameById(k).substr(0, 11) == "h1tg000002l") {
                //     std::cerr << "id: " << k << " name: " << rs_ref_.GetNameById(k) << " size: " << tmp[k - i].size() << std::endl;
                // }
            }
        };
        MultiThreads(opt_.threads, func_sketch);
        for(auto &iter: tmp) {
            for(auto &kmer: iter) {
                auto &minimizer = minimizers[kmer.hash & mask];
                if(minimizer.capacity() == minimizer.size()) minimizer.reserve(minimizer.capacity() * 1.5);
                minimizer.emplace_back(kmer);
            }
        }
        i = j;
        std::vector<std::vector<Kmer>>().swap(tmp);
    }
    
    std::vector<std::pair<std::size_t, std::size_t>> pairs(minimizers.size());
    std::atomic<std::size_t> index { 0 };
    auto func_add = [&](std::size_t tid) {
        for(auto i = index.fetch_add(1); i < minimizers.size(); i = index.fetch_add(1)) {
            if(minimizers[i].empty()) {
                pairs[i] = std::make_pair(0, 0);
                return;
            }
            RadixSort(minimizers[i].begin(), minimizers[i].end(), opt_.k * 2, Kmer::SortByHash);
            minimizers[i].emplace_back(-1, -1);
            std::size_t num_hashes = 0;
            std::size_t num_keys = 0;
            for(uint64_t j = 1, c = 1; j < minimizers[i].size(); ++j, ++c) {
                if(minimizers[i][j - 1].hash != minimizers[i][j].hash) {
                    if(c > 1) num_hashes += c;
                    ++num_keys;
                    c = 0;
                }
            }
            pairs[i] = std::make_pair(num_hashes, num_keys);
        }
    };
    MultiThreads(opt_.threads, func_add);
    for(std::size_t k = 0; k < minimizers.size(); ++k) {
        if(minimizers[k].empty()) continue;
        index_[k].index_.reserve(pairs[k].second);
        index_[k].offset_.reserve(pairs[k].first);
        for(uint64_t l = 1, c = 1; l < minimizers[k].size(); ++l, ++c) {
            if(minimizers[k][l - 1].hash != minimizers[k][l].hash) {
                if(c == 1) {
                    index_[k].index_.emplace(minimizers[k][l - 1].hash << 1 | 1, minimizers[k][l - 1].pos);
                } else {
                    index_[k].index_.emplace(minimizers[k][l - 1].hash << 1, index_[k].offset_.size() << 32 | c);
                    for(uint64_t m = l - c; m < l; ++m) {
                        index_[k].offset_.emplace_back(minimizers[k][m].pos);
                    }
                    // if(c >= 10000) {
                    //     std::cerr << "hash: " << minimizers[k][l - 1].hash << " c: " << c << std::endl;
                    // }
                }
                c = 0;
            }
        }
        std::vector<Kmer>().swap(minimizers[k]);
    }
    std::vector<std::pair<std::size_t, std::size_t>>().swap(pairs);
    std::vector<std::vector<Kmer>>().swap(minimizers);
}

void HiC::StatIndex() {
    std::array<std::size_t, 3> count = {0, 0, 0};
    for(auto &iter: index_) {
        for(auto &item: iter.index_) {
            if(item.first & 1) count[0]++;
            count[1] += 1;
            count[2] += item.first & 1? 1: static_cast<uint32_t>(item.second);
        }
    }
    LOG(INFO)("distinct kmers: %zd singletons: %zd total: %zd\n", count[1], count[0], count[2]);
    
    count = {0, 0, 0};
    for(auto &iter: index_all_) {
        for(auto &item: iter.index_) {
            if(item.first & 1) count[0]++;
            count[1] += 1;
            count[2] += item.first & 1? 1: static_cast<uint32_t>(item.second);
        }
    }
    LOG(INFO)("distinct kmers: %zd singletons: %zd total: %zd\n", count[1], count[0], count[2]);

    count = {0, 0, 0};
    for(auto &iter: index_) {
        for(auto &item: iter.index_) {
            if(!(item.first & 1)) continue;
            const uint64_t *list = nullptr;
            uint64_t hash_key = item.first >> 1;
            int cnt = index_all_[hash_key & (index_.size() - 1)].Find(hash_key, &list);
            if(cnt == 0) ++count[0];
            else if(cnt == 1) ++count[1];
            else ++count[2];
        }
    }
    LOG(INFO)("zero %zd one %zd more %zd\n", count[0], count[1], count[2]);
    
    // std::ofstream ofs("ctg.txt");
    count = {0, 0, 0};
    std::vector<std::size_t> ctg (rs_ref_.Size(), 0);
    std::unordered_map<std::string, std::vector<std::size_t>> ctg_name;
    for(auto &iter: index_all_) {
        for(auto &item: iter.index_) {
            if(item.first & 1) {
                const uint64_t *list = nullptr;
                uint64_t hash_key = item.first >> 1;
                int cnt = index_[hash_key & (index_.size() - 1)].Find(hash_key, &list);
                if(cnt == 0) {
                    ++count[0];
                    uint64_t val = item.second;
                    ctg[val >> 32]++;
                    // if(rs_ref_.GetNameById(val >> 32) == "CM039035.1") {
                    //     ofs << ((uint32_t)val >> 1)  << "\n";
                    // }
                    ctg_name[rs_ref_.GetNameById(val >> 32)].emplace_back((uint32_t)val >> 1);
                }
                else if(cnt == 1) ++count[1];
                else ++count[2];
            }
        }
    }
    LOG(INFO)("zero %zd one %zd more %zd\n", count[0], count[1], count[2]);
    // std::cerr << "ctg: \n";
    // for(std::size_t i = 0; i < ctg.size(); ++i) {
    //     if(rs_ref_.GetNameById(i).substr(0, 2) != "CM") continue;
    //     std::cerr << rs_ref_.GetNameById(i) << ": " << ctg[i] << "\n";
    // }
    // for(auto ctg: ctg_name) {
    //     std::sort(ctg.second.begin(), ctg.second.end());
    //     std::size_t prev = ctg.second[0];
    //     for(std::size_t i = 1; i < ctg.second.size(); ++i) {
    //         if(ctg.second[i] - ctg.second[i - 1] > 30) {
    //             ofs << ctg.first << "\t" << prev << "\t" << ctg.second[i - 1] << "\t" << ctg.second[i - 1] - prev << "\n";
    //             prev = ctg.second[i];
    //         }
    //         if(i == ctg.second.size() - 1) {
    //             ofs << ctg.first << "\t" << prev << "\t" << ctg.second[i] << "\t" << ctg.second[i] - prev << "\n";
    //         }
    //     }
    // }
}

int HiC::ExactMatch(const std::string &query, std::size_t q_beg, const std::string &target, std::size_t t_beg, uint64_t rev, uint64_t dir) {
    std::size_t i = 0;
    SerialDNATable ftable;
    if(rev == 0) {
        if(dir == 0) {
            for(i = 0; i < query.size() && q_beg < query.size() && t_beg < target.size(); ++i) {
                // std::cerr << "rev == 0, dir == 0 query[" << q_beg << "]: " << query[q_beg] << " target[" << t_beg << "]: " << target[t_beg] << std::endl;
                if(ftable[query[q_beg++]] != ftable[target[t_beg++]]) return i;
            }
        } else {
            for(i = 0; i < query.size() && q_beg > 0 && t_beg > 0; ++i) {
                // std::cerr << "rev == 0, dir == 1 query[" << q_beg << "]: " << query[q_beg] << " target[" << t_beg << "]: " << target[t_beg] << std::endl;
                if(ftable[query[q_beg--]] != ftable[target[t_beg--]]) return i;
            }
        }
    } else {
        ComplementDNATable rtable;
        if(dir == 0) {
            for(i == 0; i < query.size() && q_beg < query.size() && t_beg < target.size(); ++i) {
                // std::cerr << "rev == 1, dir == 0 query[" << q_beg << "]: " << query[q_beg] << " target[" << target.size() - t_beg - 1 << "]: " << target[target.size() - t_beg - 1] << std::endl;
                // std::cerr << "rev == 1, dir == 0 ftable[query[" << q_beg << "]]: " << (int)ftable[query[q_beg]] << " rtable[target[" << target.size() - t_beg - 1 << "]]: " << (int)rtable[target[target.size() - t_beg - 1]] << std::endl;
                if(ftable[query[q_beg]] != rtable[target[target.size() - t_beg - 1]]) return i;
                ++q_beg;
                ++t_beg;
            }
        } else {
            for(i = 0; i < query.size() && q_beg > 0 && t_beg > 0; ++i) {
                // std::cerr << "rev == 1, dir == 1 query[" << q_beg << "]: " << query[q_beg] << " target[" << target.size() - t_beg - 1 << "]: " << target[target.size() - t_beg - 1] << std::endl;
                // std::cerr << "rev == 1, dir == 1 ftable[query[" << q_beg << "]]: " << (int)ftable[query[q_beg]] << " rtable[target[" << target.size() - t_beg - 1 << "]]: " << (int)rtable[target[target.size() - t_beg - 1]] << std::endl;
                if(ftable[query[q_beg]] != rtable[target[target.size() - t_beg - 1]]) return i;
                --q_beg;
                --t_beg;
            }
        }
    }
    return i;
}

void HiC::GetLongestHit(const std::string &seq, std::size_t i, int strand, std::vector<HiC::Candidate> &cand, const uint64_t *list, int cnt, uint64_t &suffix) {
    suffix = (uint64_t)-1;
    Candidate c;
    std::vector<Candidate> tmp;
    // if(debug) {
    //     std::cerr << "GetLongestHit i: " << i << " cnt: " << cnt << " strand: " << strand << std::endl;
    //     std::cerr << "GetLongestHit seq: " << seq << std::endl;
    // }
    for(int j = 0; j < cnt; ++j) {
        uint64_t rev = (list[j] & 0x1) != strand;
        uint64_t ref_pos = (uint32_t)list[j] >> 1;
        ref_pos += opt_.k - 1;
        uint64_t id = list[j] >> 32;
        uint64_t ref_len = rs_ref_.GetSeqLength(id);
        // if(debug) {
        //     std::cerr << "GetLongestHit rev: " << rev << " ref_pos: " << ref_pos << " id: " << id << " ref_len: " << ref_len << std::endl;
        //     std::cerr << "GetLongestHit ref: " << (*rs_ref_.GetSeq(id).ToString()).substr(ref_pos - 30, 250) << std::endl;
        // }
        if(rev) ref_pos = ref_len - 1 - (ref_pos - opt_.k + 1);
        c.offset = i | ((uint64_t)opt_.k << 32);
        c.ref = ref_pos >= i? (ref_pos - i): (i - ref_pos) + ((uint32_t)1 << 31);
        c.ref = (rev << 63) | (id << 32) | c.ref;
        int k_len = ExactMatch(seq, i + 1, *rs_ref_.GetSeq(id).ToString(), ref_pos + 1, rev, 0);
        if(k_len < suffix) suffix = k_len;
        c.offset += ((uint64_t)k_len << 32) + k_len;
        if(i >= opt_.k && ref_pos >= opt_.k) {
            int l_len = ExactMatch(seq, i - opt_.k, *rs_ref_.GetSeq(id).ToString(), ref_pos - opt_.k, rev, 1);
            c.offset += ((uint64_t)l_len << 32);
            // if(debug) std::cerr << "i: " << i << " k_len: " << k_len << " l_len: " << l_len << " suffix: " << suffix << " offset: " << c.offset << std::endl;
        }
        tmp.emplace_back(c);
    }
    std::sort(tmp.begin(), tmp.end(), [](const Candidate &lhs, const Candidate &rhs) {
        return lhs.offset < rhs.offset;
    });
    uint64_t max_p = 0, max_p_occ = 0;
    for(int j = 0, k = 1; k <= cnt; ++k) {
        if(k == cnt || tmp[k].offset != tmp[j].offset) {
            if((max_p >> 32) < (tmp[j].offset >> 32)) {
                max_p = tmp[j].offset;
                max_p_occ = k - j;
            } else if(((max_p >> 32) == (tmp[j].offset >> 32)) && ((k - j) > max_p_occ)) {
                max_p = tmp[j].offset;
                max_p_occ = k - j;
            }
            j = k;
        }
    }
    for(int j = 0; j < cnt; ++j) {
        if(tmp[j].offset == max_p) {
            cand.emplace_back(tmp[j]);
        }
    }
}

template <typename RandomAccessIterator>
uint64_t HiC::Collect(RandomAccessIterator first, RandomAccessIterator last) {
    if(first >= last) return 0;
    if(--last == first) return (first->offset >> 32);
    uint64_t cur_beg, cur_end, beg, end, ovlp = 0, tlen = 0;
    cur_end = (uint32_t)last->offset;
    cur_beg = cur_end + 1 - (last->offset >> 32);
    for(auto iter = last - 1; iter >= first; --iter) {
        end = (uint32_t)iter->offset;
        beg = end + 1 - (iter->offset >> 32);
        if(std::max(cur_beg, beg) <= std::min(cur_end, end)) {
            cur_beg = std::min(cur_beg, beg);
        } else {
            ovlp += (cur_end + 1 - cur_beg);
            cur_beg = beg;
            cur_end = end;
        }
    }
    ovlp += (cur_end + 1 - cur_beg);
    tlen = (uint32_t)last->offset + 1 - cur_beg;
    tlen -= ovlp;
    tlen = tlen << 16;
    return (ovlp | tlen);
}

void HiC::InterprePos(HiC::Candidate &cand, uint64_t &rev, uint64_t &id, uint64_t &ref_p, uint64_t &self_p, uint64_t &exact_len, uint64_t *total_len) {
    rev = cand.ref >> 63;
    id = cand.ref << 1 >> 33;
    self_p = (uint32_t)cand.offset;
    exact_len = (cand.offset >> 32) & ((uint64_t)65535);
    if(total_len != NULL) {
        (*total_len) = (cand.offset >> 48) + exact_len;
    }
    if(cand.ref << 32 >> 63) {
        ref_p = self_p - (cand.ref & 0x7fffffff);
    } else {
        ref_p = self_p + (uint32_t)cand.ref;
    }
}

#define ISUPDATE(ml, mr, cl, cr) (((ml) < (cl)) || ((ml) == (cl) && (mr) < (cr)))

std::vector<HiC::Candidate> HiC::CompressMappedPos(std::vector<HiC::Candidate> &candidates, uint64_t thres) {
    if(candidates.empty()) return candidates;
    uint64_t rev, id, ref_p, self_p, elen, tlen;
    uint64_t max_beg = 0, max_end = 0, max_i, max_occ, cur_beg, cur_end, ovlp;
    uint64_t second_i = (uint64_t)-1, second_occ;
    uint64_t max_elen, sec_elen;
    double max_erate, sec_erate, erate;
    std::sort(candidates.begin(), candidates.end(), [](Candidate &lhs, Candidate &rhs) {
        return lhs.offset < rhs.offset;
    });
    max_elen = 0;
    max_i = (uint64_t)-1;
    max_occ = 0;
    max_erate = -1;
    for(std::size_t i = 0, j = 1; j <= candidates.size(); ++j) {
        if(j == candidates.size() || candidates[j].offset != candidates[i].offset) {
            InterprePos(candidates[i], rev, id, ref_p, self_p, elen, &tlen);
            if(debug) {
                std::cerr << "CompressMappedPos i: " << i << " j: " << j << " id: " << id << " rev: " << rev << " ref_p: " << ref_p << " self_p: " << self_p << " elen: " << elen << " tlen: " << tlen << std::endl;
            }
            erate = (double)elen / tlen;
            if(ISUPDATE(max_elen, max_erate, elen, erate)) {
                max_elen = elen;
                max_erate = erate;
                max_end = self_p;
                max_beg = self_p + 1 - tlen;
                max_i = i;
                max_occ = j - i;
            }
            i = j;
        }
    }
    sec_elen = 0;
    second_i = (uint64_t)-1;
    second_occ = 0;
    sec_erate = -1;
    if(debug) {
        std::cerr << "CompressMappedPos max_i: " << max_i << " max_occ: " << max_occ << " max_beg: " << max_beg << " max_end: " << max_end << " max_elen: " << max_elen << " max_erate: " << max_erate << std::endl;
    }
    for(std::size_t i = 0, j = 1; j <= candidates.size(); ++j) {
        if(j == candidates.size() || candidates[i].offset != candidates[j].offset) {
            if(i != max_i) {
                InterprePos(candidates[i], rev, id, ref_p, self_p, elen, &tlen);
                if(debug) std::cerr << "CompressMappedPos i: " << i << " j: " << j << " id: " << id << " rev: " << rev << " ref_p: " << ref_p << " self_p: " << self_p << " elen: " << elen << " tlen: " << tlen << std::endl;
                erate = (double)elen / tlen;
                cur_end = self_p;
                cur_beg = self_p + 1 - tlen;
                if(std::max(cur_beg, max_beg) <= std::min(cur_end, max_end)) {
                    ovlp = std::min(cur_end, max_end) - std::max(cur_beg, max_beg) + 1;
                    if(ovlp == std::min(max_end + 1 - max_end, tlen)) {
                        i = j;
                        continue;
                    }
                    if(ovlp > ((max_end + 1 - max_end) * 0.8) && elen > (max_elen * 0.8)) {
                        return std::vector<Candidate>();
                    }
                    if(ovlp > ((max_end + 1 - max_end) * 0.15) + 1) {
                        i = j;
                        continue;
                    }
                }
                if(ISUPDATE(sec_elen, sec_erate, elen, erate)) {
                    sec_elen = elen;
                    sec_erate = erate;
                    second_i = i;
                    second_occ = j - i;
                }
            }
            i = j;
        }
    }
    if(debug) {
        std::cerr << "CompressMappedPos second_i: " << second_i << " second_occ: " << second_occ << " sec_elen: " << sec_elen << " sec_erate: " << sec_erate << std::endl;
    }
    std::vector<Candidate> res;
    if(second_i == (uint64_t)-1) {
        for(uint64_t i = max_i; i < max_i + max_occ; ++i) {
            res.emplace_back(candidates[i]);
        }
    } else {
        if(max_i < second_i) {
            for(uint64_t i = max_i; i < max_i + max_occ; ++i) {
                res.emplace_back(candidates[i]);
            }
            for(uint64_t i = second_i; i < second_i + second_occ; ++i) {
                res.emplace_back(candidates[i]);
            }
        } else {
            for(uint64_t i = second_i; i < second_i + second_occ; ++i) {
                res.emplace_back(candidates[i]);
            }
            for(uint64_t i = max_i; i < max_i + max_occ; ++i) {
                res.emplace_back(candidates[i]);
            }
        }
    }
    std::vector<Candidate>().swap(candidates);
    return res;
}

std::vector<HiC::Candidate> HiC::GetAlignment(const std::string &seq, std::size_t id) {
    std::array<uint64_t, 4> minimer = {0, 0, 0, 0};
    uint64_t mask = (1ULL << opt_.k) - 1;
    uint64_t index_mask = index_.size() - 1;
    uint64_t shift = opt_.k - 1;
    uint64_t suffix;
    std::vector<Candidate> candidates;
    static SerialDNATable dna_table;
    if(debug) {
        std::cerr << "[" << GetCurTime() << "] GetAlignment id: " << id << " seq: " << seq << "\n";
    }
    for(std::size_t i = 0, l = 0; i < seq.size(); ++i) {
        uint8_t c = dna_table[seq[i]];
        uint64_t dir;
        if(c < 4 && c >= 0) {
            // minimer[0] & minimer[1] are forward minimizers
            // minimer[2] & minimer[3] are reverse minimizers
            minimer[0] = (minimer[0] << 2 | c) & mask;
            minimer[1] = (minimer[1] >> 2) | ((3 - c) << shift);
            minimer[2] = (minimer[2] << 2 | c) & mask;
            minimer[3] = (minimer[3] >> 2) | ((3 - c) << shift);
            
            if(++l >= opt_.k) {
                const uint64_t *list = nullptr;
                const uint64_t *list1 = nullptr;
                uint64_t hash_key = HashLongKmer(minimer, dir, opt_.k);
                int cnt1 = index_all_[hash_key & index_mask].Find(hash_key, &list1);
                if(cnt1 <= 0 || cnt1 > opt_.cnt_thres) continue;
                int cnt = index_[hash_key & index_mask].Find(hash_key, &list);
                if(debug) {
                    for(int j = 0; j < cnt; ++j) {
                        std::cerr << "GetAlignment kmer: " << seq.substr(i - opt_.k + 1, opt_.k) << " hash_key: " << i - opt_.k + 1 << " " << hash_key << " cnt: " << cnt << " id: " << (list[j] >> 32) << " pos: " << ((uint32_t)list[j] >> 1) << " strand: " << (list[j] & 1) << "\n";
                    }
                }
                if(cnt <= 0 || cnt > opt_.cnt_thres) continue;
                // if(debug) std::cerr << "kmer: " << seq.substr(i - opt_.k + 1, opt_.k) << " hash_key: " << hash_key << " " << hash_key << " cnt: " << cnt << std::endl;
                GetLongestHit(seq, i, dir, candidates, list, cnt, suffix);
                // if(debug) {
                //     std::cerr << "[" << GetCurTime() << "] GetAlignment i: " << i << " cnt: " << cnt << " suffix: " << suffix << " kmer: " << seq.substr(i - opt_.k + 1, opt_.k) << " strand: " << strand << " dir: " << hash_key << "\n";
                // }
                if(suffix != (uint64_t)-1) {
                    if(suffix + 1 >= opt_.k) {
                        l = 0;
                        minimer[0] = minimer[1] = 0;
                        i += suffix - (opt_.k - 1);
                    } else {
                        l = opt_.k - suffix - 1;
                    }
                }
            }
        } else {
            l = 0;
            minimer[0] = minimer[1] = 0;
        }
    }
    if(debug) {
        std::cerr << "[" << GetCurTime() << "] GetAlignment size: " << candidates.size() << "\n";
    }
    if(candidates.empty()) return candidates;
    std::sort(candidates.begin(), candidates.end(), [](const Candidate &lhs, const Candidate &rhs) {
        return lhs.ref < rhs.ref;
    });
    if(debug) {
        for(auto c: candidates) {
            uint64_t rev, id, ref_p, self_p, elen, tlen;
            InterprePos(c, rev, id, ref_p, self_p, elen, &tlen);
            std::cerr << "GetAlignment rev: " << rev << " id: " << id << " ref_p: " << ref_p << " self_p: " << self_p << " elen: " << elen << " tlen: " << tlen << " " << (c.offset >> 32) << " " << (uint32_t)c.offset << std::endl;
        }
    }
    uint64_t cur_ref_p, thres = seq.size() * 0.01 + 1, index_beg, ovlp;
    uint64_t rev, rid, ref_pos, self_pos, cnt;
    std::vector<Candidate> res;
    for(std::size_t i = 0; i < candidates.size(); ++i) {
        cur_ref_p = candidates[i].ref;
        index_beg = i;
        while(i < candidates.size() && ((candidates[i].ref >> 31) == (cur_ref_p >> 31)) && (candidates[i].ref - cur_ref_p <= thres)) i++;
        if(i - index_beg > 1) std::sort(candidates.begin() + index_beg, candidates.begin() + i, [](Candidate &lhs, Candidate &rhs) {
            return (uint32_t)lhs.offset < (uint32_t)rhs.offset;
        });
        // ovlp = Collect(candidates.data() + index_beg, i - index_beg);
        ovlp = Collect(candidates.begin() + index_beg, candidates.begin() + i);
        if(debug) std::cerr << "ovlp: " << ovlp << std::endl;
        res.emplace_back(candidates[i - 1]);
        res.back().offset = (res.back().offset << 32) >> 32;
        res.back().offset += ((uint64_t)ovlp << 32);
        i -= 1;
    }
    if(debug) {
        std::cerr << "[" << GetCurTime() << "] GetAlignment res.size(): " << res.size() << "\n";
    }
    res = CompressMappedPos(res, (opt_.k * 0.1) > 0? (opt_.k * 0.1): 1);
    if(debug) {
        std::cerr << "[" << GetCurTime() << "] GetAlignment res.size(): " << res.size() << "\n";
        for(auto c: res) {
            uint64_t rev, id, ref_p, self_p, elen, tlen;
            InterprePos(c, rev, id, ref_p, self_p, elen, &tlen);
            std::cerr << "GetAlignment rev: " << rev << " id: " << id << " ref_p: " << ref_p << " self_p: " << self_p << " elen: " << elen << " tlen: " << tlen << " " << (c.offset >> 32) << " " << (uint32_t)c.offset << std::endl;
        }
    }
    return res;
}

// std::vector<HiC::Candidate> HiC::GetAlignment(const std::string &seq, std::size_t id) {
//     uint64_t minimer[2] = {0, 0};
//     uint64_t mask = (1ULL << (opt_.k * 2)) - 1;
//     uint64_t index_mask = index_.size() - 1;
//     uint64_t shift = (opt_.k - 1) * 2;
//     uint64_t suffix;
//     std::vector<Candidate> candidates;
//     static SerialDNATable dna_table;
//     if(debug) {
//         std::cerr << "[" << GetCurTime() << "] GetAlignment id: " << id << " seq: " << seq << "\n";
//     }
//     for(std::size_t i = 0, l = 0; i < seq.size(); ++i) {
//         uint8_t c = dna_table[seq[i]];
//         if(c < 4 && c >= 0) {
//             minimer[0] = (minimer[0] << 2 | c) & mask;
//             minimer[1] = (minimer[1] >> 2) | (3ULL ^ c) << shift;
//             int strand = minimer[0] < minimer[1]? 0: 1;
//             if(++l >= opt_.k) {
//                 const uint64_t *list = nullptr;
//                 const uint64_t *list1 = nullptr;
//                 uint64_t hash_key = Hash(minimer[strand], mask);
//                 int cnt1 = index_all_[hash_key & index_mask].Find(hash_key, &list1);
//                 if(cnt1 <= 0 || cnt1 > opt_.cnt_thres) continue;
//                 int cnt = index_[hash_key & index_mask].Find(hash_key, &list);
//                 if(debug) {
//                     for(int j = 0; j < cnt; ++j) {
//                         std::cerr << "GetAlignment kmer: " << seq.substr(i - opt_.k + 1, opt_.k) << " hash_key: " << i - opt_.k + 1 << " " << minimer[strand] << " cnt: " << cnt << " id: " << (list[j] >> 32) << " pos: " << ((uint32_t)list[j] >> 1) << " strand: " << (list[j] & 1) << "\n";
//                     }
//                 }
//                 if(cnt <= 0 || cnt > opt_.cnt_thres) continue;
//                 // if(debug) std::cerr << "kmer: " << seq.substr(i - opt_.k + 1, opt_.k) << " hash_key: " << hash_key << " " << minimer[strand] << " cnt: " << cnt << std::endl;
//                 GetLongestHit(seq, i, strand, candidates, list, cnt, suffix);
//                 // if(debug) {
//                 //     std::cerr << "[" << GetCurTime() << "] GetAlignment i: " << i << " cnt: " << cnt << " suffix: " << suffix << " kmer: " << seq.substr(i - opt_.k + 1, opt_.k) << " strand: " << strand << " hash: " << hash_key << "\n";
//                 // }
//                 if(suffix != (uint64_t)-1) {
//                     if(suffix + 1 >= opt_.k) {
//                         l = 0;
//                         minimer[0] = minimer[1] = 0;
//                         i += suffix - (opt_.k - 1);
//                     } else {
//                         l = opt_.k - suffix - 1;
//                     }
//                 }
//             }
//         } else {
//             l = 0;
//             minimer[0] = minimer[1] = 0;
//         }
//     }
//     if(debug) {
//         std::cerr << "[" << GetCurTime() << "] GetAlignment size: " << candidates.size() << "\n";
//     }
//     if(candidates.empty()) return candidates;
//     std::sort(candidates.begin(), candidates.end(), [](const Candidate &lhs, const Candidate &rhs) {
//         return lhs.ref < rhs.ref;
//     });
//     if(debug) {
//         for(auto c: candidates) {
//             uint64_t rev, id, ref_p, self_p, elen, tlen;
//             InterprePos(c, rev, id, ref_p, self_p, elen, &tlen);
//             std::cerr << "GetAlignment rev: " << rev << " id: " << id << " ref_p: " << ref_p << " self_p: " << self_p << " elen: " << elen << " tlen: " << tlen << " " << (c.offset >> 32) << " " << (uint32_t)c.offset << std::endl;
//         }
//     }
//     uint64_t cur_ref_p, thres = seq.size() * 0.01 + 1, index_beg, ovlp;
//     uint64_t rev, rid, ref_pos, self_pos, cnt;
//     std::vector<Candidate> res;
//     for(std::size_t i = 0; i < candidates.size(); ++i) {
//         cur_ref_p = candidates[i].ref;
//         index_beg = i;
//         while(i < candidates.size() && ((candidates[i].ref >> 31) == (cur_ref_p >> 31)) && (candidates[i].ref - cur_ref_p <= thres)) i++;
//         if(i - index_beg > 1) std::sort(candidates.begin() + index_beg, candidates.begin() + i, [](Candidate &lhs, Candidate &rhs) {
//             return (uint32_t)lhs.offset < (uint32_t)rhs.offset;
//         });
//         // ovlp = Collect(candidates.data() + index_beg, i - index_beg);
//         ovlp = Collect(candidates.begin() + index_beg, candidates.begin() + i);
//         if(debug) std::cerr << "ovlp: " << ovlp << std::endl;
//         res.emplace_back(candidates[i - 1]);
//         res.back().offset = (res.back().offset << 32) >> 32;
//         res.back().offset += ((uint64_t)ovlp << 32);
//         i -= 1;
//     }
//     if(debug) {
//         std::cerr << "[" << GetCurTime() << "] GetAlignment res.size(): " << res.size() << "\n";
//     }
//     res = CompressMappedPos(res, (opt_.k * 0.1) > 0? (opt_.k * 0.1): 1);
//     if(debug) {
//         std::cerr << "[" << GetCurTime() << "] GetAlignment res.size(): " << res.size() << "\n";
//         for(auto c: res) {
//             uint64_t rev, id, ref_p, self_p, elen, tlen;
//             InterprePos(c, rev, id, ref_p, self_p, elen, &tlen);
//             std::cerr << "GetAlignment rev: " << rev << " id: " << id << " ref_p: " << ref_p << " self_p: " << self_p << " elen: " << elen << " tlen: " << tlen << " " << (c.offset >> 32) << " " << (uint32_t)c.offset << std::endl;
//         }
//     }
//     return res;
// }

void HiC::Get53List(std::vector<HiC::Candidate> &cand, std::vector<HiC::Candidate> &l5, std::vector<HiC::Candidate> &l3) {
    l3.clear();
    l5.clear();
    uint64_t rev, id, ref_p, self_p, elen, tlen, cur_beg, num;
    // uint64_t beg_5 = (uint64_t)-1;
    // // if(debug) {
    // //     std::cerr << "Get53List cand.size(): " << cand.size() << std::endl;
    // // }
    // for(std::size_t i = 0, j = 1, num = 0; j <= cand.size(); ++j) {
    //     if(j == cand.size() || cand[j].offset != cand[i].offset) {
    //         InterprePos(cand[i], rev, id, ref_p, self_p, elen, &tlen);
    //         // if(debug) {
    //         //     std::cerr << "Get53List rev: " << rev << " id: " << id << " ref_p: " << ref_p << " self_p: " << self_p << " elen: " << elen << " tlen: " << tlen << std::endl;
    //         // }
    //         cur_beg = self_p + 1 - tlen;
    //         num++;
    //         if(cur_beg <= beg_5) {
    //             beg_5 = cur_beg;
    //             for(std::size_t k = i; k < j; ++k) {
    //                 l5.emplace_back(cand[k]);
    //                 l3.emplace_back(cand[k]);
    //             }
    //         } else {
    //             for(std::size_t k = i; k < j; ++k) {
    //                 l3.emplace_back(cand[k]);
    //             }
    //         }
    //         i = j;
    //     }
    // }
    uint64_t max_len = (uint64_t)1;
    // if(debug) {
    //     std::cerr << "Get53List cand.size(): " << cand.size() << std::endl;
    // }
    std::sort(cand.begin(), cand.end(), [](const Candidate &lhs, const Candidate &rhs) {
        return (lhs.offset >> 32) > (rhs.offset >> 32);
    });
    for(std::size_t i = 0, j = 1, num = 0; j <= cand.size(); ++j) {
        if(j == cand.size() || cand[j].offset != cand[i].offset) {
            InterprePos(cand[i], rev, id, ref_p, self_p, elen, &tlen);
            // if(debug) {
            //     std::cerr << "Get53List rev: " << rev << " id: " << id << " ref_p: " << ref_p << " self_p: " << self_p << " elen: " << elen << " tlen: " << tlen << std::endl;
            // }
            num++;
            if(max_len < tlen) {
                max_len = tlen;
                for(std::size_t k = i; k < j; ++k) {
                    l5.emplace_back(cand[k]);
                    l3.emplace_back(cand[k]);
                }
            } else {
                for(std::size_t k = i; k < j; ++k) {
                    l3.emplace_back(cand[k]);
                }
            }
            i = j;
        }
    }
}

void HiC::SetPos(std::vector<HiC::Candidate> &cand1, std::vector<HiC::Candidate> &cand2, uint64_t rid, HiC::PeHit &ph) {
    if(cand1.empty() || cand2.empty()) return;
    uint64_t rev, id, ref_p, self_p, elen, tlen;
    std::vector<Candidate> l1_5;
    std::vector<Candidate> l1_3;
    std::vector<Candidate> l2_5;
    std::vector<Candidate> l2_3;
    Get53List(cand1, l1_5, l1_3);
    Get53List(cand2, l2_5, l2_3); 
    if(l1_5.empty() || l2_5.empty()) return;
    ph.id = rid;
    ph.len = 0;
    for(std::size_t i = 0; i < l1_5.size(); ++i) {
        InterprePos(l1_5[i], rev, id, ref_p, self_p, elen, &tlen);
        if(ref_p + 1 < tlen) continue;
        ref_p = ref_p + 1 - tlen;
        if(rev) ref_p = rs_ref_.GetSeqLength(id) - 1 - ref_p;
        ph.s = (rev << 63) | (id << 32) | ref_p;
        ph.len = tlen;
        ph.len = ph.len << 32;
        // if(debug) {
        //     std::cerr << "SetPos l1_5 id: " << id << " ref_p: " << ref_p << " rev: " << rev << " elen: " << elen << " tlen: " << tlen << std::endl;
        // }
    }
    for(std::size_t i = 0; i < l2_5.size(); ++i) {
        InterprePos(l2_5[i], rev, id, ref_p, self_p, elen, &tlen);
        if(ref_p + 1 < tlen) continue;
        ref_p = ref_p + 1 - tlen;
        if(rev) ref_p = rs_ref_.GetSeqLength(id) - 1 - ref_p;
        ph.e = (rev << 63) | (id << 32) | ref_p;
        ph.len |= tlen;
        // if(debug) {
        //     std::cerr << "SetPos l2_5 id: " << id << " ref_p: " << ref_p << " rev: " << rev << " elen: " << elen << " tlen: " << tlen << std::endl;
        // }
    }
    if(ph.s == (uint64_t)-1 || ph.e == (uint64_t)-1) {
        ph.id = ph.s = ph.e = ph.len = (uint64_t)-1;
    }
}

std::vector<HiC::PeHit> HiC::DedupHits(std::vector<HiC::PeHit> &hits) {
    if(hits.empty()) return hits;
    std::sort(hits.begin(), hits.end(), [](const PeHit &lhs, const PeHit &rhs) {
        return lhs.s < rhs.s;
    });
    PeHit cur;
    std::vector<PeHit> res;
    for(std::size_t i = 0, j = 1; j <= hits.size(); ++j) {
        if(j == hits.size() || hits[j].s != hits[i].s) {
            if(j - i > 1) std::sort(hits.begin() + i, hits.begin() + j, [](const PeHit &lhs, const PeHit &rhs) {
                return lhs.e < rhs.e;
            });
            cur.s = cur.e = (uint64_t)-1;
            while(i < j) {
                if(hits[i].e != cur.e || hits[i].s != cur.s) {
                    cur = hits[i];
                    if((hits[i].s << 1 >> 33) != (hits[i].e) << 1 >> 33) {
                        res.emplace_back(hits[i]);
                    }
                }
                i++;
            }
            i = j;
        }
    }
    hits.clear();
    return res;
}

void HiC::Run() {
    assert(opt_.h1fname.size() == opt_.h2fname.size());

    LOG(INFO)("Load reference...\n");
    for(auto &rfname: opt_.rfname) {
        rs_ref_.Load(rfname);
    }
    LOG(INFO)("Load reference finished. CPU time: %.3f sec, real time: %.3f sec, peak RSS: %.3f GB\n", GetCPUTime(), GetRealTime(), PeakMemoryGB());

    LOG(INFO)("Load variants...\n");
    LoadVariants();
    LOG(INFO)("Load variants finished. CPU time: %.3f sec, real time: %.3f sec, peak RSS: %.3f GB\n", GetCPUTime(), GetRealTime(), PeakMemoryGB());

    LOG(INFO)("Index begining...\n");
    BuildIndex();
    // BuildIndexSingle();
    LOG(INFO)("Index finished. CPU time: %.3f sec, real time: %.3f sec, peak RSS: %.3f GB\n", GetCPUTime(), GetRealTime(), PeakMemoryGB());
    BuildIndexAll();

    StatIndex();
    
    // std::cerr << "h2tg000337l:9322-312589: " << rs_ref_.GetIdByName("h2tg000337l:9322-312589") << std::endl;
    // std::cerr << "h1tg000086l:987216-1066088: " << rs_ref_.GetIdByName("h1tg000086l:987216-1066088") << std::endl;

    std::size_t rid = 0;
    std::vector<PeHit> hits;
    hic_name_.reserve(1000000);
    for(std::size_t i = 0; i < opt_.h1fname.size(); ++i) {
        SeqReader *reader1 = ReadStore::OpenFile(opt_.h1fname[i]);
        SeqReader *reader2 = ReadStore::OpenFile(opt_.h2fname[i]);
        while(true) {
            SeqReader::Unit u1, u2;
            std::vector<SeqReader::Unit> unit1, unit2;
            unit1.reserve(1ULL << 20);
            unit2.reserve(1ULL << 20);
            uint64_t total_size = 0;
            while(reader1->Next(u1) && reader2->Next(u2) && total_size < (1ULL << 31)) {
                if(u1.seq.size() < opt_.k || u2.seq.size() < opt_.k) continue;
                assert(u1.head == u2.head);
                u1.id = u2.id = rid++;
                // if(u1.head == "A00675:51:HJVCKDRXX:1:2117:12933:4883") {
                //     std::cerr << "id: " << u1.id << " name: " << u1.head << std::endl;
                //     std::cerr << "seq1: " << u1.seq << std::endl;
                //     std::cerr << "qual: " << u1.quality << std::endl;
                //     std::cerr << "seq2: " << u2.seq << std::endl;
                //     std::cerr << "qual: " << u2.quality << std::endl;
                //     return;
                // }
                // continue;
                unit1.emplace_back(u1);
                unit2.emplace_back(u2);
                hic_name_.emplace_back(u1.head);
                total_size += u1.seq.size() + u2.seq.size();
            }
            if(total_size == 0) break;
            std::vector<PeHit> tmp(unit1.size());
            std::atomic<std::size_t> index { 0 };
            auto func = [&](std::size_t tid) {
                for(auto cur = index.fetch_add(1); cur < unit1.size(); cur = index.fetch_add(1)) {
                    auto seq1 = unit1[cur];
                    auto seq2 = unit2[cur];
                    // if(seq1.id < 252325148 || seq1.id > 252325148) continue;
                    // if(seq1.id == 252325148) debug = true;
                    // else debug = false;
                    assert(seq1.id == seq2.id);
                    auto cand1 = GetAlignment(seq1.seq, seq1.id);
                    // std::cerr << "[" << GetCurTime() << "] Run 1 id: " << seq1.id << " name: " << seq1.head << " size: " << cand1.size() << " vmRSS: " << GetMemoryUsage() << std::endl;
                    if(cand1.empty()) continue;
                    // debug = true;
                    auto cand2 = GetAlignment(seq2.seq, seq2.id);
                    // std::cerr << "[" << GetCurTime() << "] Run 2 id: " << seq2.id << " name: " << seq2.head << " size: " << cand2.size() << " vmRSS: " << GetMemoryUsage() << std::endl;
                    if(cand2.empty()) continue;
                    tmp[cur].s = tmp[cur].e = tmp[cur].id = tmp[cur].len = (uint64_t)-1;
                    SetPos(cand1, cand2, seq1.id, tmp[cur]);
                    // std::cerr << "[" << GetCurTime() << "] Run id: " << seq1.id << " name: " << seq1.head << " id1: " << (tmp[cur].s << 1 >> 33) << " strand: " << (tmp[cur].s >> 63) << " pos1: " << (uint32_t)tmp[cur].s << " len1: " << (tmp[cur].len >> 32) << " id2: " << (tmp[cur].e << 1 >> 33) << " strand: " << (tmp[cur].e >> 63) << " pos2: " << (uint32_t)tmp[cur].e << " len2: " << (uint32_t)tmp[cur].len << std::endl;
                }
            };
            MultiThreads(opt_.threads, func);
            if(hits.capacity() < (hits.size() + tmp.size())) {
                hits.reserve(hits.size() + tmp.size());
            }
            for(auto &ph: tmp) {
                if(ph.id == (uint64_t)-1) continue;
                hits.emplace_back(ph);
            }
            unit1.clear();
            unit2.clear();
        }
    }
    hits = DedupHits(hits);
    LOG(INFO)("%zd / %zd hits, vmRSS: %.3f GB\n", hits.size(), rid, PeakMemoryGB());

    std::string ofname = opt_.prefix_ + ".hic.txt";
    std::ofstream out(ofname);
    if(!out.is_open()) {
        LOG(ERROR)("Failed to open file %s\n", ofname.c_str());
    }
    // hic name | name1 | strand1 | pos1 | len1 | name2 | strand2 | pos2 | len2
    for(auto &h: hits) {
        // std::cerr << h.id << " " << hic_name_[h.id] << " " << (h.s << 1 >> 33) << " " << (h.s >> 63) << " " << (uint32_t)h.s << " " << (h.len >> 32) << " " << (h.e << 1 >> 33) << " " << (h.e >> 63) << " " << (uint32_t)h.e << " " << (uint32_t)h.len << std::endl;
        out << hic_name_[h.id] << "\t" << rs_ref_.GetNameById((h.s << 1 >> 33)) << "\t" << (h.s >> 63) << "\t" << (uint32_t)h.s << "\t" << (h.len >> 32) << "\t" << rs_ref_.GetNameById((h.e << 1 >> 33)) << "\t" << (h.e >> 63) << "\t" << (uint32_t)h.e << "\t" << (uint32_t)h.len << "\n";
    }
    out.close();
    LOG(INFO)("Finished. CPU time: %.3f sec, real time: %.3f sec, peak RSS: %.3f GB\n", GetCPUTime(), GetRealTime(), PeakMemoryGB());
}

std::vector<std::string> GetFiles(const std::string &str) {
    std::vector<std::string> files;
    auto begin = str.find_first_not_of(' ');
    while(begin != std::string::npos) {
        auto end = str.find_first_of(',', begin);
        if(end == std::string::npos) end = str.size();
        assert(begin < end);
        files.emplace_back(str.substr(begin, end - begin));
        if(end == str.size()) break;
        begin = str.find_first_not_of(' ', end + 1);
    }
    return files;
}

void Usage(Option &opt) {
    std::cerr << "Usage: HiC [options]" << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << "  -v, --variants  FILE     variants file" << std::endl;
    std::cerr << "  -r, --reference FILE     reference file" << std::endl;
    std::cerr << "  -t, --threads   INT      number of threads [" << opt.threads << "]" << std::endl;
    std::cerr << "  -k, --kmer      INT      kmer size [" << opt.k << "]" << std::endl;
    std::cerr << "  -s, --threshold INT      threshold of kmer count [" << opt.cnt_thres << "]" << std::endl;
    std::cerr << "  -h, --help               print help information" << std::endl;
    std::cerr << "      --hic1      FILE     Hi-C mate-pair 1 file" << std::endl;
    std::cerr << "      --hic2      FILE     Hi-C mate-pair 2 file" << std::endl;
}

int ParseArgument(int argc, char **argv, Option &opt) {
    struct option option_long[] = {
        {"variants",    required_argument,  0,  'v'},
        {"reference",   required_argument,  0,  'r'}, 
        {"threads",     optional_argument,  0,  't'},
        {"kmer",        optional_argument,  0,  'k'}, 
        {"threshold",   optional_argument,  0,  's'}, 
        {"hic1",        required_argument,  0,  301},
        {"hic2",        required_argument,  0,  302}, 
        {"help",        no_argument,        0,  'h'},
        {0,             0,                  0,  0}
    };
    const char *option_short = "v:r:t:k:s:h";
    int c;
    while((c = getopt_long(argc, argv, option_short, option_long, NULL)) != -1) {
        if(c == 'v') opt.vfname = GetFiles(std::string(optarg));
        else if(c == 'r') opt.rfname = GetFiles(std::string(optarg));
        else if(c == 't') opt.threads = atoi(optarg);
        else if(c == 'k') opt.k = atoi(optarg);
        else if(c == 's') opt.cnt_thres = atoi(optarg);
        else if(c == 301) opt.h1fname = GetFiles(std::string(optarg));
        else if(c == 302) opt.h2fname = GetFiles(std::string(optarg));
        else if(c == 'h') {
            Usage(opt);
            return 0;
        } else if(c == '?') {
            LOG(WARNING)("Unknown option: %c\n", (char)optopt);
            return 1;
        } else if(c == ':') {
            LOG(WARNING)("Missing option argument: %c\n", (char)optopt);
            return 1;
        }
    }
    return opt.Check();
}

int main(int argc, char **argv) {
    LOG(INFO)("Program begin.\n");
    ResetRealTime();
    // LOG(INFO)("Program begin. CPU time: %.3f sec, real time: %.3f sec, peak RSS: %.3f GB\n", GetCPUTime(), GetRealTime(), PeakMemoryGB());
    Option opt;
    if(ParseArgument(argc, argv, opt) == 1) {
        Usage(opt);
        return 0;
    }
    HiC hic(opt);
    hic.Run();
    LOG(INFO)("Program finished. CPU time: %.3f sec, real time: %.3f sec, peak RSS: %.3f GB\n", GetCPUTime(), GetRealTime(), PeakMemoryGB());
    return 0;
}