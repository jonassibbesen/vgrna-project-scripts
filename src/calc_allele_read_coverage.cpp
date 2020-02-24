
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

vector<string> parseGenotype(const string & sample) {

    auto genotype_str = splitString(sample, ':');

    if (genotype_str.front().find('/') != std::string::npos) {

        assert(genotype_str.front().find('|') == std::string::npos);
        return splitString(genotype_str.front(), '/');

    } else {

        return splitString(genotype_str.front(), '|');

    }
}

string alleleIdxToSequence(const uint32_t allele_idx, const vector<string> & variant) {

    if (allele_idx == 0) {

        return variant.at(3);

    } else {

        return splitString(variant.at(4), ',').at(allele_idx - 1);
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

void trimAlleles(string * ref_allele, string * alt_allele) {

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
}

string getAlleleType(string ref_allele, string alt_allele) {

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

int main(int argc, char* argv[]) {

    if (argc != 4) {

        cout << "Usage: calc_allele_read_coverage <vcf_name> <bam_name> <mapq_threshold> > coverage.txt" << endl;
        return 1;
    }

    printScriptHeader(argc, argv);

    ifstream vcf_istream(argv[1]);
    assert(vcf_istream.is_open());

    BamReader bam_reader;
    bam_reader.Open(argv[2]);
    assert(bam_reader.IsOpen());

    uint32_t mapq_threshold = stoi(argv[3]);

    unordered_map<string, uint32_t> allele_coverage_stats;

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

        uint32_t num_reads_h1 = 0;
        uint32_t num_reads_h2 = 0;

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

            stringstream allele_coverage_ss;
            allele_coverage_ss << num_reads_h1;
            allele_coverage_ss << "\t" << first_allele_type;
            allele_coverage_ss << "\t" << num_reads_h2;
            allele_coverage_ss << "\t" << second_allele_type;
        
            auto allele_coverage_stats_it = allele_coverage_stats.emplace(allele_coverage_ss.str(), 0);
            allele_coverage_stats_it.first->second++;
        }

        if (num_variants % 10000 == 0) {

            cerr << "Number of analysed variants: " << num_variants << endl;
        }        
    }

    vcf_istream.close();
    bam_reader.Close();

    cout << "Count" << "\t" << "NumReadsH1" << "\t" << "AlleleTypeH1" << "\t" << "NumReadsH2" << "\t" << "AlleleTypeH2" << endl;

    for (auto & stats: allele_coverage_stats) {

        cout << stats.second << "\t" << stats.first << endl;
    }

	return 0;
}
