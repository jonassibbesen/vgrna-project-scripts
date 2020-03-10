
/*
calc_read_transcript_overlap_stats
Calculates overlapping statistics between mapped reads and 
transcripts in bam format.
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

unordered_map<string, SeqLib::GRC> createTranscriptGenomicRegions(const string & transcript_bam_file, const uint32_t min_deletion_length) {

    unordered_map<string, SeqLib::GRC> transcript_genomic_regions;

    BamReader bam_reader;
    bam_reader.Open(transcript_bam_file);
    assert(bam_reader.IsOpen());

    BamRecord bam_record;

    while (bam_reader.GetNextRecord(bam_record)) {

        assert(bam_record.GetCigar().NumQueryConsumed() == bam_record.Length());

        auto transcript_genomic_regions_it = transcript_genomic_regions.emplace(bam_record.ChrName(bam_reader.Header()), SeqLib::GRC());
        cigarToGenomicRegions(&(transcript_genomic_regions_it.first->second), bam_record.GetCigar(), bam_record.Position(), min_deletion_length);
    }

    bam_reader.Close();

    for (auto & chrom_regions: transcript_genomic_regions) {

        chrom_regions.second.MergeOverlappingIntervals();
        chrom_regions.second.CreateTreeMap();
    }

    return transcript_genomic_regions;
}

int main(int argc, char* argv[]) {

    if (argc != 5) {

        cout << "Usage: calc_read_transcript_overlap_stats <read_bam_name> <transcript_bam_name> <min_deletion_length> <debug_output (use 0 or 1)> > statistics.txt" << endl;
        return 1;
    }

    printScriptHeader(argc, argv);

    BamReader bam_reader;
    bam_reader.Open(argv[1]);
    assert(bam_reader.IsOpen());

    uint32_t min_deletion_length = stoi(argv[3]);

    auto transcript_genomic_regions = createTranscriptGenomicRegions(argv[2], min_deletion_length);
    uint32_t total_width = 0;

    for (auto & chrom_regions: transcript_genomic_regions) {

        total_width += chrom_regions.second.TotalWidth();
    }

    cerr << "Number of transcript bases: " << total_width << "\n" << endl;

    bool debug_output = stoi(argv[4]);

    if (debug_output) {

        cout << "Name" << "\t" << "MapQ" << "\t" << "ReadLength" << "\t" << "ReadNumSoftClip" << "\t" << "Overlap" << endl;
    }

    BamRecord bam_record;

    unordered_map<string, uint32_t> overlap_stats;

    uint32_t num_reads = 0;
    float sum_overlap = 0;

    while (bam_reader.GetNextRecord(bam_record)) { 

        if (bam_record.SecondaryFlag()) {

            continue;
        }

        num_reads++;

        uint32_t read_num_soft_clip = 0;
        float overlap = 0;

        auto transcript_genomic_regions_it = transcript_genomic_regions.find(bam_record.ChrName(bam_reader.Header()));

        if (bam_record.MappedFlag()) { 

            assert(bam_record.GetCigar().NumQueryConsumed() == bam_record.Length());
            read_num_soft_clip = numSoftClippedBases(bam_record.GetCigar());

            if (transcript_genomic_regions_it != transcript_genomic_regions.end()) {

                auto read_cigar_genomic_regions = cigarToGenomicRegions(bam_record.GetCigar(), bam_record.Position(), 0);
                read_cigar_genomic_regions.CreateTreeMap();

                auto cigar_genomic_regions_intersection = transcript_genomic_regions_it->second.Intersection(read_cigar_genomic_regions, true);

                overlap = cigar_genomic_regions_intersection.TotalWidth() / static_cast<float>(bam_record.Length()); 
            }
        }

        if (debug_output) {

            cout << bam_record.Qname();
            cout << "\t" << bam_record.MapQuality();
            cout << "\t" << bam_record.Length();
            cout << "\t" << read_num_soft_clip;
            cout << "\t" << overlap;
            cout << endl;

        } else {

            stringstream overlap_stats_ss;
            overlap_stats_ss << bam_record.MapQuality();
            overlap_stats_ss << "\t" << bam_record.Length();
            overlap_stats_ss << "\t" << read_num_soft_clip;
            overlap_stats_ss << "\t" << overlap;            

            auto overlap_stats_it = overlap_stats.emplace(overlap_stats_ss.str(), 0);
            overlap_stats_it.first->second++;  
        }

        sum_overlap += overlap;

        if (num_reads % 1000000 == 0) {

            cerr << "Number of analysed reads: " << num_reads << endl;
        }        
    }

    bam_reader.Close();

    if (!debug_output) {

        cout << "Count" << "\t" << "MapQ" << "\t" << "ReadLength" << "\t" << "ReadNumSoftClip" << "\t" << "Overlap" << endl;

        for (auto & stats: overlap_stats) {

            cout << stats.second << "\t" << stats.first << endl;
        }
    }

    cerr << "\nTotal number of analysed reads: " << num_reads << endl;
    cerr << "Average overlap: " << sum_overlap/num_reads << endl;

	return 0;
}