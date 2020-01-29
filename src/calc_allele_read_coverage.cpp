
/*
calc_allele_rsem_read_coverage
Calculates mapped read coverage across alleles in a diploid vcf file. 
The read names should end with the haplotype origin (e.g. "_h1")
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <assert.h>
#include <omp.h>

#include "SeqLib/BamReader.h"

#include "utils.hpp"

using namespace SeqLib;

int main(int argc, char* argv[]) {

    if (argc != 4) {

        cout << "Usage: calc_allele_read_coverage <vcf_input_name> <bam_input_name> <mapq_threshold> > output.txt" << endl;
        return 1;
    }

    ifstream vcf_istream(argv[1]);
    assert(vcf_istream.is_open());

    BamReader bam_reader;
    bam_reader.Open(argv[2]);
    assert(bam_reader.IsOpen());

    int32_t mapq_threshold = stoi(argv[3]);

    string line;
    BamRecord bam_record;

    int32_t num_variants = 0;

    while (vcf_istream.good()) {

        getline(vcf_istream, line);

        if (line.empty() || line.front() == '#') {

            continue;
        }

        num_variants++;

        if (num_variants % 10000 == 0) {

            cerr << num_variants << endl;
        }

        auto line_split = splitString(line, '\t');
        assert(line_split.size() == 10);

        auto genotype = parseGenotype(line_split.at(9));
        assert(genotype.size() == 2);

        GenomicRegion genomic_region(line_split.at(0), to_string(stoi(line_split.at(1)) - 1), to_string(stoi(line_split.at(1)) + line_split.at(3).size() - 1), bam_reader.Header());
        assert(bam_reader.SetRegion(genomic_region));

        int32_t num_reads_h1 = 0;
        int32_t num_reads_h2 = 0;

        while (bam_reader.GetNextRecord(bam_record)) { 

            if (bam_record.MapQuality() >= mapq_threshold && !bam_record.SecondaryFlag()) {

                auto read_origin = bam_record.Qname().substr(bam_record.Qname().size() - 2);

                if (read_origin == "h1") {

                    num_reads_h1 += 1;
                
                } else {

                    assert(read_origin == "h2");
                    num_reads_h2 += 1;
                }
            }
        }

        if (num_reads_h1 + num_reads_h2 > 0) {

            auto first_allele_type = getAlleleType(line_split.at(3), alleleIdxToSequence(stoi(genotype.front()), line_split));
            auto second_allele_type = getAlleleType(line_split.at(3), alleleIdxToSequence(stoi(genotype.back()), line_split));

            cout << num_reads_h1 << "\t" << first_allele_type << "\t" << num_reads_h2 << "\t" << second_allele_type << endl;
        }
    }

    vcf_istream.close();
    bam_reader.Close();

	return 0;
}

