#include "sequence.hpp"

SerialDNATable DNASeq::serial_dna_table_;
ComplementDNATable DNASeq::complement_dna_table_;

std::string Seq::ReverseComplement(const std::string &seq) {
    std::string s(seq.rbegin(), seq.rend());
    std::for_each(s.begin(), s.end(), [](char &c) { c = DNASeq::Complement(c); });
    return s;
}

void Seq::ReverseComplementInPlace(std::string &seq) {
    std::reverse(seq.begin(), seq.end());
    std::for_each(seq.begin(), seq.end(), [](char &c) { c = DNASeq::Complement(c); });
}