
#ifndef VGRNA_PROJECT_SCRIPTS_UTILS_HPP
#define VGRNA_PROJECT_SCRIPTS_UTILS_HPP

#include <string>
#include <vector>
#include <sstream>
#include <assert.h>

#include "vg/io/basic_stream.hpp"
#include "google/protobuf/util/json_util.h"

using namespace std;

//------------------------------------------------------------------------------

/*
The following code was copied and modified from https://github.com/vgteam/vg
*/

inline int mapping_to_length(const vg::Mapping & m) {
    int l = 0;
    for (int i = 0; i < m.edit_size(); ++i) {
        const vg::Edit& e = m.edit(i);
        l += e.to_length();
    }
    return l;
}

inline int mapping_from_length(const vg::Mapping & m) {
    int l = 0;
    for (int i = 0; i < m.edit_size(); ++i) {
        const vg::Edit& e = m.edit(i);
        l += e.from_length();
    }
    return l;

}

inline string pb2json(const google::protobuf::Message &msg) {
    // Set options to preserve field names and not camel case them
    google::protobuf::util::JsonPrintOptions opts;
    opts.preserve_proto_field_names = true;

	string buffer;
    auto status = google::protobuf::util::MessageToJsonString(msg, &buffer, opts);
    
    if (!status.ok()) {
        throw runtime_error("Could not serialize " + msg.GetTypeName() + ": " + status.ToString());
    }
    
    return buffer;
}

inline void json2pb(google::protobuf::Message &msg, const string& buf) {
    auto status = google::protobuf::util::JsonStringToMessage(buf, &msg);
    
    if (!status.ok()) {
        // This generally will happen if someone feeds in the wrong type of JSON.
        // TODO: It would be nice to be able to find the neme of the offending non-existent field.
        throw runtime_error("Could not deserialize " + msg.GetTypeName() + ": " + status.ToString());
    }
}

//------------------------------------------------------------------------------


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

inline vector<string> parseGenotype(const string & sample) {

    auto genotype_str = splitString(sample, ':');

    if (genotype_str.front().find('/') != std::string::npos) {

        assert(genotype_str.front().find('|') == std::string::npos);
        return splitString(genotype_str.front(), '/');

    } else {

        return splitString(genotype_str.front(), '|');

    }
}

inline string alleleIdxToSequence(const int32_t allele_idx, const vector<string> & variant) {

    if (allele_idx == 0) {

        return variant.at(3);

    } else {

        return splitString(variant.at(4), ',').at(allele_idx - 1);
    }
}

inline void rightTrim(string * allele, const int32_t trim_length) {

    if (trim_length >= allele->size()) {

        *allele = "";
    
    } else {

        *allele = allele->substr(0, allele->size() - trim_length);
    }
}

inline void leftTrim(string * allele, const int32_t trim_length) {

    if (trim_length >= allele->size()) {

        *allele = "";
    
    } else {

        *allele = allele->substr(trim_length);
    }
}

inline void trimAlleles(string * ref_allele, string * alt_allele) {

    int32_t right_trim_len = 0;

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

    int32_t left_trim_len = 0;

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
}

inline string getAlleleType(string ref_allele, string alt_allele) {

    trimAlleles(&ref_allele, &alt_allele);

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
