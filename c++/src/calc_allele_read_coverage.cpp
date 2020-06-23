
/*
calc_allele_rsem_read_coverage
Calculates mapped read coverage across alleles in a diploid vcf file. 
The read names should end with the haplotype origin (e.g. "_h1")
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <assert.h>

#include "SeqLib/BamReader.h"
#include "SeqLib/BamRecord.h"

#include "utils.hpp"

using namespace SeqLib;

int main(int argc, char* argv[]) {

    if (argc != 4) {

        cout << "Usage: calc_allele_read_coverage <read_bam> <variant_vcf> <allele_idx (1-based)> > coverage.txt" << endl;
        return 1;
    }

    printScriptHeader(argc, argv);

    BamReader bam_reader;
    bam_reader.Open(argv[1]);
    assert(bam_reader.IsOpen());

    ifstream vcf_istream(argv[2]);
    assert(vcf_istream.is_open());

    uint32_t allele_idx = stoi(argv[3]) - 1;
    assert(allele_idx < 2);

    cout << "Count" << "\t" << "MapQ" << "\t" << "AllelePosition" << "\t" << "AlleleType" << "\t" << "RelativeAlleleLength" << endl;

    string line;
    BamRecord bam_record;

    uint32_t num_variants = 0;

    while (vcf_istream.good()) {

        getline(vcf_istream, line);

        if (line.empty() || line.front() == '#') {

            continue;
        }

        num_variants++;

        auto line_split = splitString(line, '\t');
        assert(line_split.size() == 10);

        auto genotype = parseGenotype(line_split.at(9));
        assert(genotype.size() == 2);

        GenomicRegion genomic_region(line_split.at(0), to_string(stoi(line_split.at(1)) - 1), to_string(stoi(line_split.at(1)) + line_split.at(3).size() - 1), bam_reader.Header());
        assert(bam_reader.SetRegion(genomic_region));

        unordered_map<uint32_t, uint32_t> mapq_read_counts;

        while (bam_reader.GetNextRecord(bam_record)) { 

            auto mapq_read_counts_it = mapq_read_counts.emplace(bam_record.MapQuality(), 0);
            mapq_read_counts_it.first->second++;
        }

        for (auto & mapq_count: mapq_read_counts) {

            cout << mapq_count.second;
            cout << "\t" << mapq_count.first;
            cout << "\t" << line_split.at(0) << ":" << stoi(line_split.at(1));
            cout << "\t" << getAlleleType(line_split.at(3), alleleIdxToSequence(genotype.at(allele_idx), line_split));
            cout << "\t" << getAlleleLength(line_split.at(3), alleleIdxToSequence(genotype.at(allele_idx), line_split));
            cout << endl;
        }

        if (num_variants % 10000 == 0) {

            cerr << "Number of analysed variants: " << num_variants << endl;
        }        
    }

    vcf_istream.close();
    bam_reader.Close();

    cerr << "\nTotal number of analysed variants: " << num_variants << endl;

	return 0;
}
