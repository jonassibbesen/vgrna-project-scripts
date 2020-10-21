
#ifndef VGRNA_PROJECT_SCRIPTS_UTILS_HPP
#define VGRNA_PROJECT_SCRIPTS_UTILS_HPP

#include <string>
#include <vector>
#include <sstream>
#include <assert.h>

#include "SeqLib/BamRecord.h"
#include "SeqLib/GenomicRegion.h"
#include "SeqLib/GenomicRegionCollection.h"

using namespace std;

// Precision used when comparing double variables.
static const double double_precision = numeric_limits<double>::epsilon() * 100;

// Compare double variables using above precision.
bool doubleCompare(const double a, const double b) {

    assert(isfinite(a));
    assert(isfinite(b));

    return ((a == b) or (abs(a - b) < abs(min(a, b)) * double_precision));
}

template<class T>
ostream & operator<<(ostream & os, const vector<T> & values) {

    auto values_it = values.cbegin();

    if (values_it == values.cend()) {

        return os;
    }

    os << *values_it;
    ++values_it;

    while (values_it != values.cend()) {

        os << " " << *values_it;
        ++values_it;
    }

    return os;
}

vector<string> splitString(const string & str, const char delim) {

    stringstream ss(str);
    vector<string> elems;

    for (string item; getline(ss, item, delim);) {

        elems.push_back(item);
    }

    return elems;
}

void printScriptHeader(int argc, char * argv[]) {

    vector<string> argv_vec;
    argv_vec.assign(argv, argv + argc);

    cerr << GIT_COMMIT << endl;
    cerr << argv_vec << "\n" << endl;
}

pair<SeqLib::Cigar, uint32_t> trimCigar(const SeqLib::Cigar & cigar, const uint32_t trim_start, const uint32_t trim_length) {

    assert(trim_start + trim_length <= cigar.NumQueryConsumed());

    if (trim_start > 0 || trim_length < cigar.NumQueryConsumed()) { 

        SeqLib::Cigar trimmed_cigar;
        uint32_t num_shifted_genomic_bases = 0;

        uint32_t cur_query_length = 0;

        for (auto & field: cigar) {

            if (field.ConsumesQuery()) {

                uint32_t new_field_length = 0;

                if (cur_query_length <= trim_start && trim_start <= cur_query_length + field.Length() - 1) {

                    assert(trimmed_cigar.TotalLength() == 0);

                    new_field_length = min(static_cast<uint32_t>(trim_length), cur_query_length + field.Length() - trim_start);
                    assert(new_field_length > 0);

                    if (field.ConsumesReference()) {

                        num_shifted_genomic_bases += (trim_start - cur_query_length);
                    }

                } else if (trimmed_cigar.size() > 0) {

                    new_field_length = min(static_cast<uint32_t>(trim_length - trimmed_cigar.NumQueryConsumed()), field.Length());
                    assert(new_field_length > 0);
                }

                assert(new_field_length <= field.Length());
                cur_query_length += field.Length();

                if (new_field_length > 0) {

                    trimmed_cigar.add(SeqLib::CigarField(field.Type(), new_field_length));
                
                    if (trimmed_cigar.NumQueryConsumed() == trim_length) {

                        break;
                    }   
                } 

            } else if (trimmed_cigar.size() > 0) {

                trimmed_cigar.add(field); 
            }
            
            if (field.ConsumesReference() && trimmed_cigar.size() == 0) {

                num_shifted_genomic_bases += field.Length();
            }
        }
    
        assert(trimmed_cigar.NumQueryConsumed() == trim_length); 
        return make_pair(trimmed_cigar, num_shifted_genomic_bases);

    } else {

        return make_pair(cigar, 0);
    }
}

void cigarToGenomicRegions(SeqLib::GRC * cigar_genomic_regions, const SeqLib::Cigar & cigar, const uint32_t chrom_idx, const uint32_t start_pos) {
    
    uint32_t cur_length = 0;

    for (auto & field: cigar) {

        if (field.Type() == 'M' || field.Type() == '=' || field.Type() == 'X') {

            cigar_genomic_regions->add(SeqLib::GenomicRegion(chrom_idx, start_pos + cur_length, start_pos + cur_length + field.Length() - 1));
            cur_length += field.Length();

        } else if (field.Type() == 'D' || field.Type() == 'N') {

            cur_length += field.Length();
        }
    }
}

SeqLib::GRC cigarToGenomicRegions(const SeqLib::Cigar & cigar, const uint32_t chrom_idx, const uint32_t start_pos) {

    SeqLib::GRC cigar_genomic_regions;
    cigarToGenomicRegions(&cigar_genomic_regions, cigar, chrom_idx, start_pos);

    return cigar_genomic_regions;
}

uint32_t cigarTypeLength(const SeqLib::Cigar & cigar, const char type) {

    uint32_t type_length = 0;

    for (auto & field: cigar) {

        if (field.Type() == type) {

            type_length += field.Length();
        }
    }

    return type_length;
}

string genomicRegionsToString(const SeqLib::GRC & genomic_regions) {

    stringstream genomic_regions_ss; 

    bool is_first = true;
    for (auto & region: genomic_regions.AsGenomicRegionVector()) {

        if (!is_first) {

            genomic_regions_ss << ",";
        }

        genomic_regions_ss << region.pos1 + 1 << "-" << region.pos2 + 1;
        is_first = false;
    }

    return genomic_regions_ss.str();
}

vector<string> parseGenotype(const string & sample) {

    auto genotype_str = splitString(sample, ':');

    if (genotype_str.front().find('/') != std::string::npos) {

        assert(genotype_str.front().find('|') == std::string::npos);
        return splitString(genotype_str.front(), '/');

    } else {

        return splitString(genotype_str.front(), '|');

    }
}

string alleleIdxToSequence(const string allele_idx_str, const vector<string> & variant) {

    if (allele_idx_str == ".") {

        return "";        
    
    } else {

        auto allele_idx = stoi(allele_idx_str); 
        
        if (allele_idx == 0) {

            return variant.at(3);

        } else {

            return splitString(variant.at(4), ',').at(allele_idx - 1);
        }
    }
}

void rightTrim(string * allele, const uint32_t trim_length) {

    if (trim_length >= allele->size()) {

        *allele = "";
    
    } else {

        *allele = allele->substr(0, allele->size() - trim_length);
    }
}

void leftTrim(string * allele, const uint32_t trim_length) {

    if (trim_length >= allele->size()) {

        *allele = "";
    
    } else {

        *allele = allele->substr(trim_length);
    }
}

uint32_t trimAlleles(string * ref_allele, string * alt_allele) {

    assert(!ref_allele->empty());
    assert(!alt_allele->empty());

    if (*ref_allele == *alt_allele) {

        *ref_allele = "";
        *alt_allele = "";

        return 1;
    }

    uint32_t right_trim_len = 0;

    for (size_t i = 0; i < min(ref_allele->size(), alt_allele->size()); ++i) {

        if (ref_allele->at(ref_allele->size() - i - 1) == alt_allele->at(alt_allele->size() - i - 1)) {

            right_trim_len++;
        
        } else {

            break;
        }
    }

    if (right_trim_len > 0) {

        rightTrim(ref_allele, right_trim_len);
        rightTrim(alt_allele, right_trim_len);
    }

    uint32_t left_trim_len = 0;

    for (size_t i = 0; i < min(ref_allele->size(), alt_allele->size()); ++i) {

        if (ref_allele->at(i) == alt_allele->at(i)) {

            left_trim_len++;
        
        } else {

            break;
        }
    }

    if (left_trim_len > 0) {

        leftTrim(ref_allele, left_trim_len);
        leftTrim(alt_allele, left_trim_len);
    }

    return left_trim_len;
}

string getAlleleType(string ref_allele, string alt_allele) {

    if (ref_allele == alt_allele) {

        return "REF";

    } else if (ref_allele.size() == alt_allele.size() and ref_allele.size() == 1) {

        return "SNV";

    } else if (ref_allele.empty()) {

        assert(!alt_allele.empty());
        return "INS";

    } else if (alt_allele.empty()) {

        assert(!ref_allele.empty());
        return "DEL";        

    } else {

        assert(!ref_allele.empty());
        assert(!alt_allele.empty());
        return "COM";
    }
}


#endif
