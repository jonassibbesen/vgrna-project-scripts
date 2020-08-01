
/*
convert_snaptron_to_bed
Convert introns and filter introns based on number 
of samples and read coverage.
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <assert.h>

#include "utils.hpp"

using namespace SeqLib;


int main(int argc, char* argv[]) {

    if (argc != 4) {

        cout << "Usage: convert_snaptron_to_bed <introns_tsv_name> <min_num_samples> <min_num_reads_per_sample> > introns.bed" << endl;
        return 1;
    }

    printScriptHeader(argc, argv);

    ifstream introns_istream(argv[1]);
    assert(introns_istream.is_open());

    uint32_t min_num_samples = stoi(argv[2]);
    assert(min_num_samples > 0);

    uint32_t min_num_reads_per_sample = stoi(argv[3]);
    assert(min_num_reads_per_sample > 0);

    string line;

    uint32_t num_introns = 0;
    uint32_t num_introns_filt = 0;

    while (introns_istream.good()) {

        getline(introns_istream, line, '\n');

        if (line.empty() || line.front() == '#') {

            continue;
        }

        num_introns++;

        auto line_split = splitString(line, '\t');
        auto samples_split = splitString(line_split.at(11), ',');

        assert(samples_split.front() == "");
        assert(samples_split.size() > 1);

        if ((samples_split.size() - 1) < min_num_samples) {

            ++num_introns_filt;
            continue;
        }

        uint32_t num_read_samples = 0;

        auto samples_split_it = samples_split.begin();
        assert(samples_split_it != samples_split.end());

        ++samples_split_it;
        assert(samples_split_it != samples_split.end());

        while (samples_split_it != samples_split.end()) {

            auto sample_split = splitString(*samples_split_it, ':');
            
            assert(sample_split.size() == 2);
            assert(stoi(sample_split.back()) >= 1);

            if (stoi(sample_split.back()) >= min_num_reads_per_sample) {

                ++num_read_samples;
            }

            ++samples_split_it;
        }

        if (num_read_samples < min_num_samples) {

            ++num_introns_filt;
            continue;
        }

        cout << line_split.at(1) << "\t" << stoi(line_split.at(2)) - 1 << "\t" << line_split.at(3) << "\t" << line_split.at(5) << endl;
    }

    introns_istream.close();

    cerr << "\nNumber of introns: " << num_introns << endl;
    cerr << "Number of introns filtered: " << num_introns_filt << endl;

	return 0;
}
