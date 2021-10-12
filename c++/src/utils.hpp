
#ifndef VGRNA_PROJECT_SCRIPTS_UTILS_HPP
#define VGRNA_PROJECT_SCRIPTS_UTILS_HPP

#include <string>
#include <vector>
#include <sstream>
#include <assert.h>
#include <sys/stat.h>

#include "SeqLib/BamRecord.h"
#include "SeqLib/GenomicRegion.h"
#include "SeqLib/GenomicRegionCollection.h"
#include "htslib/vcf.h"
#include "htslib/tbx.h"


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

    if (trim_start > 0 || trim_length < cigar.NumQueryConsumed()) { 

        SeqLib::Cigar trimmed_cigar;
        uint32_t num_shifted_genomic_bases = 0;

        uint32_t cur_query_length = 0;

        for (auto & field: cigar) {

            if (field.ConsumesQuery()) {

                uint32_t new_field_length = 0;

                if (cur_query_length <= trim_start && trim_start <= cur_query_length + field.Length() - 1) {

                    assert(trimmed_cigar.TotalLength() == 0);

                    new_field_length = min(trim_length, cur_query_length + field.Length() - trim_start);
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
    
        assert(trimmed_cigar.NumQueryConsumed() <= trim_length); 
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

vector<tuple<htsFile*, bcf_hdr_t*, tbx_t*, int>> initializeVCFs(const vector<string>& vcf_filenames,
                                                                const string& sample_name,
                                                                unordered_map<string, int>& contig_to_vcf_out) {
    
    vector<tuple<htsFile*, bcf_hdr_t*, tbx_t*, int>> vcfs;
    for (const string& vcf_filename : vcf_filenames) {
        
        // make sure the VCF and tabis exist
        string tabix_filename = vcf_filename + ".tbi";
        struct stat stat_tbi, stat_vcf;
        if (stat(vcf_filename.c_str(), &stat_vcf) != 0) {
            cerr << "VCF file " << vcf_filename << " cannot be opened" << endl;
            exit(1);
        }
        if (stat(tabix_filename.c_str(), &stat_tbi) != 0) {
            cerr << "Tabix file " << tabix_filename << " cannot be opened. Must tabix index VCF file " << vcf_filename << " before running benchmark." << endl;
            exit(1);
        }
        
        // load them up
        htsFile* vcf = bcf_open(vcf_filename.c_str(), "r");
        assert(vcf);
        
        bcf_hdr_t* header = bcf_hdr_read(vcf);
        assert(header);
        
        tbx_t* tabix_index = tbx_index_load2(tabix_filename.c_str(),
                                             vcf_filename.c_str());
        assert(tabix_index);
        
        
        // find the index of the sample we want
        // TODO: there should be a way to do this using the dictionary in the
        // header, but i can't find it...
        int idx = -1;
        for (int i = 0, n = bcf_hdr_nsamples(header); i < n; ++i) {
            if (header->samples[i] == sample_name) {
                idx = i;
                break;
            }
        }
        vcfs.emplace_back(vcf, header, tabix_index, idx);
        
        // record which contigs occur in this VCF
        int num_seq_names = 0;
        const char** contig_names = tbx_seqnames(tabix_index, &num_seq_names);
        for (int j = 0; j < num_seq_names; ++j) {
            string contig = contig_names[j];
            contig_to_vcf_out[contig] = vcfs.size() - 1;
        }
        free(contig_names);
    }
    
    return vcfs;
}

inline pair<uint32_t, uint32_t> countIndelsAndSubs(const SeqLib::BamHeader & bam_header, const SeqLib::GRC & regions,
                                                   htsFile * vcf, bcf_hdr_t * bcf_header, tbx_t * tabix_index,
                                                   int sample_idx) {
    
    uint32_t indel_bps = 0;
    uint32_t subs_bps = 0;
    
    for (const SeqLib::GenomicRegion& region : regions) {
        
        string contig = region.ChrName(bam_header);
        
        // init iteration variables
        int tid = tbx_name2id(tabix_index, contig.c_str());
        hts_itr_t * itr = tbx_itr_queryi(tabix_index, tid, region.pos1, region.pos2);
        bcf1_t * bcf_record = bcf_init();
        kstring_t kstr = {0, 0, 0};
        
        // iterate over VCF lines
        while (tbx_itr_next(vcf, tabix_index, itr, &kstr) >= 0) {
            vcf_parse(&kstr, bcf_header, bcf_record);
            
            set<int> alleles;
            // init a genotype array
            int32_t* genotypes = nullptr;
            int arr_size = 0;
            // and query it
            int num_genotypes = bcf_get_genotypes(bcf_header, bcf_record, &genotypes, &arr_size);
            int ploidy = num_genotypes / bcf_hdr_nsamples(bcf_header);
            // look at the genotype of this sample
            for (int i = ploidy * sample_idx, n = ploidy * (sample_idx + 1); i < n; ++i) {
                if (genotypes[i] == bcf_int32_vector_end) {
                    // sample has lower ploidy
                    break;
                }
                if (bcf_gt_is_missing(genotypes[i])) {
                    continue;
                }
                alleles.insert(bcf_gt_allele(genotypes[i]));
            }
            free(genotypes);
            
            if (alleles.size() > 1 || (alleles.size() == 1 && !alleles.count(0))) {
                // the sample has a non-ref allele
                
                // we want to know the range of alleles that either a) exist in the sample
                // or b) are the ref, so we add the ref allele here so that we always
                // consider it
                alleles.insert(0);
                
                // make sure the lazily-unpacked metadata is unpacked through the alleles
                if (bcf_record->unpacked & BCF_UN_INFO) {
                    bcf_unpack(bcf_record, BCF_UN_INFO);
                }
                
                int min_allele_length = numeric_limits<int>::max();
                int max_allele_length = numeric_limits<int>::min();
                for (auto i : alleles) {
                    int len = strlen(bcf_record->d.allele[i]);
                    min_allele_length = min(min_allele_length, len);
                    max_allele_length = max(max_allele_length, len);
                }
                
                if (min_allele_length == 1 && max_allele_length > 1) {
                    // indel
                    indel_bps += max_allele_length - 1;
                }
                else if (min_allele_length == max_allele_length) {
                    // substitution
                    subs_bps += max_allele_length;
                }
                // else: complex variant, we ignore it for simplicity
            }
        }
        
        bcf_destroy(bcf_record);
        tbx_itr_destroy(itr);
        if (kstr.s) {
            free(kstr.s);
        }
    }
    return make_pair(subs_bps, indel_bps);
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
