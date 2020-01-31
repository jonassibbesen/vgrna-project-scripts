
/*
calc_rsem_benchmark_stats
Calculates benchmarking statistics for rsem simulated mapped reads.
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

unordered_map<string, BamRecord> parseTranscriptAlignments(const string & transcript_bam_file) {

    unordered_map<string, BamRecord> transcript_alignments;

    BamReader bam_reader;
    bam_reader.Open(transcript_bam_file);
    assert(bam_reader.IsOpen());

    BamRecord bam_record;

    while (bam_reader.GetNextRecord(bam_record)) { 

        auto transcript_id_suffix = bam_record.Qname().substr(bam_record.Qname().size() - 2);

        if (transcript_id_suffix == "_1" || transcript_id_suffix == "_2" || transcript_id_suffix == "_3") {

            assert(transcript_alignments.emplace(bam_record.Qname().substr(0, bam_record.Qname().size() - 2), bam_record).second);
            
        } else {

            assert(transcript_alignments.emplace(bam_record.Qname(), bam_record).second);
        }
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

        auto transcript_id_suffix = element.substr(element.size() - 2);

        if (transcript_id_suffix == "_1" || transcript_id_suffix == "_2" || transcript_id_suffix == "_3") {

            transcript_ids.emplace_back(element.substr(0, element.size() - 2));
            
        } else {

            transcript_ids.emplace_back(element);
        }

        rsem_istream.ignore(numeric_limits<streamsize>::max(), '\n');
    }

    rsem_istream.close();

    return transcript_ids;
}

GRC cigarToGenomicRegionCollection(const Cigar & cigar, const uint32_t start_pos) {

    GRC cigar_genomic_regions;
    
    uint32_t match_length = 0;
    uint32_t deletion_length = 0;

    for (auto & field: cigar) {

        if (field.Type() == 'M' || field.Type() == '=' || field.Type() == 'X') {

            uint32_t cur_pos = start_pos + match_length + deletion_length;

            cigar_genomic_regions.add(GenomicRegion(0, cur_pos, cur_pos + field.Length() - 1));
            match_length += field.Length();

        } else if (field.Type() == 'D' || field.Type() == 'N') {

            deletion_length += field.Length();
        }
    }

    return cigar_genomic_regions;
}

int main(int argc, char* argv[]) {

    if (argc != 4) {

        cout << "Usage: calc_rsem_benchmark_stats <read_bam_name> <transcript_bam_name> <rsem_expression_name> > statistics.txt" << endl;
        return 1;
    }

    printScriptHeader(argc, argv);

    BamReader bam_reader;
    bam_reader.Open(argv[1]);
    assert(bam_reader.IsOpen());

    auto transcript_alignments = parseTranscriptAlignments(argv[2]);
    cerr << transcript_alignments.size() << endl;

    auto transcript_ids = parseTranscriptIds(argv[3]);
    cerr << transcript_ids.size() << endl;

    unordered_map<string, uint32_t> consensus_stats;

    BamRecord bam_record;

    uint32_t num_reads = 0;

    while (bam_reader.GetNextRecord(bam_record)) { 

        num_reads++;

        if (num_reads % 1000000 == 0) {

            cerr << num_reads << endl;
        }

        if (bam_record.SecondaryFlag()) {

            continue;
        }

        if (!bam_record.MappedFlag()) {

            assert(bam_record.MapQuality() == 0);

            auto consensus_stats_it = consensus_stats.emplace("0\t0", 0);
            consensus_stats_it.first->second++;

            continue;
        }

        auto read_name_split = splitString(bam_record.Qname(), '_');
        assert(read_name_split.size() == 6);

        auto read_transcript_id = transcript_ids.at(stoi(read_name_split.at(2)) - 1);
        uint32_t read_transcript_pos = stoi(read_name_split.at(3)) + 1;

        if (!bam_record.FirstFlag()) {

            read_transcript_pos += stoi(read_name_split.at(4)) - bam_record.Length();
        }

        auto transcript_alignments_it = transcript_alignments.find(read_transcript_id);
        assert(transcript_alignments_it != transcript_alignments.end());

        if (bam_record.ChrName() != transcript_alignments_it->second.ChrName()) {
        
            stringstream consensus_ss;
            consensus_ss << bam_record.MapQuality();
            consensus_ss << "\t0";

            auto consensus_stats_it = consensus_stats.emplace(consensus_ss.str(), 0);
            consensus_stats_it.first->second++;

            continue;
        }

        bool is_forward = false;

        if ((!transcript_alignments_it->second.ReverseFlag() && read_name_split.at(1) == "0") or (transcript_alignments_it->second.ReverseFlag() && read_name_split.at(1) =="1")) {

            is_forward = true;
        }

        if (!is_forward) {

            read_transcript_pos = transcript_alignments_it->second.Length() - (read_transcript_pos - 1) - (bam_record.Length() - 1);
        }

        if (is_forward) {

            continue;
        }

        uint32_t cur_transcript_pos = 0;
        uint32_t read_transcript_genomic_pos = transcript_alignments_it->second.Position();

        Cigar transcript_read_cigar;

        for (auto & field: transcript_alignments_it->second.GetCigar()) {

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
                } 

            } else if (transcript_read_cigar.size() > 0) {

                transcript_read_cigar.add(field);
            
            } else {

                read_transcript_genomic_pos += field.Length();
            }
        }

        read_transcript_genomic_pos += read_transcript_pos - 1; 

        auto read_cigar_genomic_regions = cigarToGenomicRegionCollection(bam_record.GetCigar(), bam_record.Position());
        auto transcript_cigar_genomic_regions = cigarToGenomicRegionCollection(transcript_read_cigar, read_transcript_genomic_pos);

        read_cigar_genomic_regions.CreateTreeMap();
        transcript_cigar_genomic_regions.CreateTreeMap();

        auto cigar_genomic_regions_intersection = read_cigar_genomic_regions.Intersection(transcript_cigar_genomic_regions, true);

        float consensus = cigar_genomic_regions_intersection.TotalWidth() * 2 / static_cast<float>(read_cigar_genomic_regions.TotalWidth() + transcript_cigar_genomic_regions.TotalWidth()); 

        stringstream consensus_ss;
        consensus_ss << bam_record.MapQuality();
        consensus_ss << "\t" << consensus;

        auto consensus_stats_it = consensus_stats.emplace(consensus_ss.str(), 0);
        consensus_stats_it.first->second++;

        // cerr << endl;
        // cerr << read_transcript_id << endl;
        // cerr << read_transcript_pos << endl;
        // cerr << bam_record;
        // cerr << transcript_read_cigar << endl;
        // cerr << read_transcript_genomic_pos << endl;
        // cerr << transcript_alignments_it->second;
        // cerr << read_cigar_genomic_regions.AsGenomicRegionVector() << endl;
        // cerr << transcript_cigar_genomic_regions.AsGenomicRegionVector() << endl;
        // cerr << cigar_genomic_regions_intersection.AsGenomicRegionVector() << endl;
        // cerr << consensus << endl;
        // cerr << endl;
    }

    bam_reader.Close();

    cout << "Count" << "\t" << "MapQ" << "\t" << "Consensus" << endl;

    for (auto & stats: consensus_stats) {

        cout << stats.second << "\t" << stats.first << endl;
    }

	return 0;
}
