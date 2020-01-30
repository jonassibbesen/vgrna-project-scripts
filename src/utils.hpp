
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

inline void printScriptHeader(int argc, char * argv[]) {

    vector<string> argv_vec;
    argv_vec.assign(argv, argv + argc);

    cerr << GIT_HEADER << endl;
    cerr << argv_vec << "\n" << endl;
}


#endif
