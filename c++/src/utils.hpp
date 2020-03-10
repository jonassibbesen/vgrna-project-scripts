
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
inline bool doubleCompare(const double a, const double b) {

    assert(isfinite(a));
    assert(isfinite(b));

    return ((a == b) or (abs(a - b) < abs(min(a, b)) * double_precision));
}

template<class T>
inline ostream & operator<<(ostream & os, const vector<T> & values) {

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

inline vector<string> splitString(const string & str, const char delim) {

    stringstream ss(str);
    vector<string> elems;

    for (string item; getline(ss, item, delim);) {

        elems.push_back(item);
    }

    return elems;
}

inline void printScriptHeader(int argc, char * argv[]) {

    vector<string> argv_vec;
    argv_vec.assign(argv, argv + argc);

    cerr << GIT_HEADER << endl;
    cerr << argv_vec << "\n" << endl;
}

inline void cigarToGenomicRegions(SeqLib::GRC * cigar_genomic_regions, const SeqLib::Cigar & cigar, const uint32_t start_pos, const uint32_t min_deletion) {
    
    uint32_t cur_length = 0;

    for (auto & field: cigar) {

        if (field.Type() == 'M' || field.Type() == '=' || field.Type() == 'X') {

            cigar_genomic_regions->add(SeqLib::GenomicRegion(0, start_pos + cur_length, start_pos + cur_length + field.Length() - 1));
            cur_length += field.Length();

        } else if (field.Type() == 'D' || field.Type() == 'N') {

            if (min_deletion > field.Length()) {

                cigar_genomic_regions->add(SeqLib::GenomicRegion(0, start_pos + cur_length, start_pos + cur_length + field.Length() - 1));
            }

            cur_length += field.Length();
        }
    }
}

inline SeqLib::GRC cigarToGenomicRegions(const SeqLib::Cigar & cigar, const uint32_t start_pos, const uint32_t min_deletion) {

    SeqLib::GRC cigar_genomic_regions;
    cigarToGenomicRegions(&cigar_genomic_regions, cigar, start_pos, min_deletion);

    return cigar_genomic_regions;
}

uint32_t numSoftClippedBases(const SeqLib::Cigar & cigar) {

    uint32_t num_soft_clipped_bases = 0;

    for (auto & field: cigar) {

        if (field.Type() == 'S') {

            num_soft_clipped_bases += field.Length();
        }
    }

    return num_soft_clipped_bases;
}

#endif
