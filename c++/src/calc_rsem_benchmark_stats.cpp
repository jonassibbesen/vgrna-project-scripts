
/*
calc_rsem_benchmark_stats
Calculates benchmarking statistics for rsem simulated reference mapped reads.
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <assert.h>

#include "SeqLib/BamReader.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/GenomicRegion.h"
#include "SeqLib/GenomicRegionCollection.h"

#include "utils.hpp"

using namespace SeqLib;

unordered_map<string, pair<string, BamRecord> > parseTranscriptAlignments(const string & transcript_bam_file) {

    unordered_map<string, pair<string, BamRecord> > transcript_alignments;

    BamReader bam_reader;
    bam_reader.Open(transcript_bam_file);
    assert(bam_reader.IsOpen());

    BamRecord bam_record;

    while (bam_reader.GetNextRecord(bam_record)) { 

        assert(bam_record.GetCigar().NumQueryConsumed() == bam_record.Length());
        assert(transcript_alignments.emplace(bam_record.Qname(), make_pair(bam_record.ChrName(bam_reader.Header()), bam_record)).second);
    }

    bam_reader.Close();

    return transcript_alignments;
}

vector<string> parseTranscriptIds(const string & rsem_expression_file) {

    vector<string> transcript_ids;

    ifstream rsem_istream(rsem_expression_file);
    assert(rsem_istream.is_open());

    string element;

    while (rsem_istream.good()) {

        getline(rsem_istream, element, '\t');

        if (element.empty() || element == "transcript_id") {

            rsem_istream.ignore(numeric_limits<streamsize>::max(), '\n');
            continue;
        }

        transcript_ids.emplace_back(element);
        rsem_istream.ignore(numeric_limits<streamsize>::max(), '\n');
    }

    rsem_istream.close();

    return transcript_ids;
}

int main(int argc, char* argv[]) {

    if (!(argc == 4 || argc == 5)) {

        cout << "Usage: calc_rsem_benchmark_stats <read_bam> <transcript_bam> <rsem_expression_file> (<enable_debug_output>) > statistics.txt" << endl;
        return 1;
    }

    printScriptHeader(argc, argv);

    BamReader bam_reader;
    bam_reader.Open(argv[1]);
    assert(bam_reader.IsOpen());

    auto transcript_alignments = parseTranscriptAlignments(argv[2]);
    cerr << "Number of transcript alignments: " << transcript_alignments.size() << endl;

    auto transcript_ids = parseTranscriptIds(argv[3]);
    cerr << "Number of transcript names: " << transcript_ids.size() << "\n" << endl;

    bool debug_output = (argc == 5);

    stringstream base_header; 
    base_header << "TruthAlignmentLength" << "\t" << "IsMapped" << "\t" << "MapQ" << "\t" << "Length" << "\t" << "SoftClipLength" << "\t" << "Overlap";

    if (debug_output) {

        cout << "Name" << "\t" << "Alignment" << "\t" << "TruthAlignment" << "\t" << base_header.str() << endl;
    } 

    BamRecord bam_record;

    unordered_map<string, uint32_t> benchmark_stats;

    uint32_t num_reads = 0;
    double sum_overlap = 0;

    while (bam_reader.GetNextRecord(bam_record)) { 

        if (bam_record.SecondaryFlag()) {

            continue;
        }

        num_reads++;

        auto read_name_split = splitString(bam_record.Qname(), '_');
        assert(read_name_split.size() == 5);

        auto read_transcript_id = transcript_ids.at(stoi(read_name_split.at(2)) - 1);
        uint32_t read_transcript_pos = stoi(read_name_split.at(3)) + 1;

        auto read_name_end_split = splitString(read_name_split.at(4), '/');
        assert(read_name_end_split.size() <= 2); 

        if (read_name_end_split.size() == 2) {

            if (read_name_end_split.back() == "2") {

                read_transcript_pos += stoi(read_name_end_split.front()) - bam_record.Length();
            }

        } else if (!bam_record.FirstFlag()) {

            read_transcript_pos += stoi(read_name_end_split.front()) - bam_record.Length();
        }

        auto transcript_alignments_it = transcript_alignments.find(read_transcript_id);
        assert(transcript_alignments_it != transcript_alignments.end());

        bool is_forward = false;

        if ((!transcript_alignments_it->second.second.ReverseFlag() && read_name_split.at(1) == "0") or (transcript_alignments_it->second.second.ReverseFlag() && read_name_split.at(1) =="1")) {

            is_forward = true;
        }

        if (!is_forward) {

            read_transcript_pos = transcript_alignments_it->second.second.Length() - (read_transcript_pos - 1) - (bam_record.Length() - 1);
        }

        uint32_t read_transcript_genomic_pos = transcript_alignments_it->second.second.Position();
        read_transcript_genomic_pos += read_transcript_pos - 1;

        uint32_t cur_transcript_pos = 0;

        Cigar transcript_read_cigar;

        for (auto & field: transcript_alignments_it->second.second.GetCigar()) {

            if (field.ConsumesQuery()) {

                uint32_t new_field_length = 0;

                if (cur_transcript_pos <= read_transcript_pos && read_transcript_pos <= cur_transcript_pos + field.Length() - 1) {

                    assert(transcript_read_cigar.TotalLength() == 0);

                    new_field_length = min(static_cast<uint32_t>(bam_record.Length()), cur_transcript_pos + field.Length() - (read_transcript_pos - 1));
                    assert(new_field_length > 0);

                } else if (transcript_read_cigar.size() > 0) {

                    new_field_length = min(static_cast<uint32_t>(bam_record.Length() - transcript_read_cigar.NumQueryConsumed()), field.Length());
                    assert(new_field_length > 0);
                }

                cur_transcript_pos += field.Length();

                if (new_field_length > 0) {

                    transcript_read_cigar.add(CigarField(field.Type(), new_field_length));
                
                    if (transcript_read_cigar.NumQueryConsumed() == bam_record.Length()) {

                        break;
                    }
                
                } else if (!field.ConsumesReference()) {

                    assert(field.Type() == 'I' || field.Type() == 'S');
                    read_transcript_genomic_pos -= field.Length();
                }

            } else if (transcript_read_cigar.size() > 0) {

                transcript_read_cigar.add(field);
            
            } else {

                assert(field.ConsumesReference());
                read_transcript_genomic_pos += field.Length();
            }
        }

        uint32_t soft_clip_length = 0;
        double overlap = 0;

        string read_genomic_regions_str = "";
        auto transcript_cigar_genomic_regions = cigarToGenomicRegions(transcript_read_cigar, 0, read_transcript_genomic_pos);

        if (bam_record.MappedFlag()) {

            assert(bam_record.GetCigar().NumQueryConsumed() == bam_record.Length());
            assert(bam_record.GetCigar().NumQueryConsumed() == transcript_read_cigar.NumQueryConsumed());

            soft_clip_length = cigarTypeLength(bam_record.GetCigar(), 'S');

            if (bam_record.ChrName(bam_reader.Header()) == transcript_alignments_it->second.first) {
        
                auto read_cigar_genomic_regions = cigarToGenomicRegions(bam_record.GetCigar(), 0, bam_record.Position());

                read_cigar_genomic_regions.CreateTreeMap();
                transcript_cigar_genomic_regions.CreateTreeMap();

                auto cigar_genomic_regions_intersection = transcript_cigar_genomic_regions.Intersection(read_cigar_genomic_regions, true);
                overlap = cigar_genomic_regions_intersection.TotalWidth() / static_cast<double>(transcript_cigar_genomic_regions.TotalWidth());

                if (debug_output) {

                    read_genomic_regions_str = genomicRegionsToString(read_cigar_genomic_regions);
                }
            }
        }

        stringstream benchmark_stats_ss;

        benchmark_stats_ss << transcript_cigar_genomic_regions.TotalWidth();
        benchmark_stats_ss << "\t" << bam_record.MappedFlag();
        benchmark_stats_ss << "\t" << bam_record.MapQuality();
        benchmark_stats_ss << "\t" << bam_record.Length();
        benchmark_stats_ss << "\t" << soft_clip_length;
        benchmark_stats_ss << "\t" << overlap;   

        if (debug_output) {

            cout << bam_record.Qname();
            cout << "\t" << bam_record.ChrName(bam_reader.Header()) << ":" << read_genomic_regions_str;
            cout << "\t" << transcript_alignments_it->second.first << ":" << genomicRegionsToString(transcript_cigar_genomic_regions);
            cout << "\t" << benchmark_stats_ss.str();
            cout << endl;

        } else {         

            auto benchmark_stats_it = benchmark_stats.emplace(benchmark_stats_ss.str(), 0);
            benchmark_stats_it.first->second++;  
        }

        sum_overlap += overlap;

        if (num_reads % 1000000 == 0) {

            cerr << "Number of analysed reads: " << num_reads << endl;
        }
    }

    bam_reader.Close();

    if (!debug_output) {

        cout << "Count" << "\t" << base_header.str() << endl;

        for (auto & stats: benchmark_stats) {

            cout << stats.second << "\t" << stats.first << endl;
        }
    }

    cerr << "\nTotal number of analysed reads: " << num_reads << endl;
    cerr << "Average overlap: " << sum_overlap/num_reads << endl;

	return 0;
}
