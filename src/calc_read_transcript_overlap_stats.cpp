
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

unordered_map<string, SeqLib::GRC> createTranscriptGenomicRegions(const string & transcript_bam_file) {

    unordered_map<string, SeqLib::GRC> transcript_genomic_regions;

    BamReader bam_reader;
    bam_reader.Open(transcript_bam_file);
    assert(bam_reader.IsOpen());

    BamRecord bam_record;

    while (bam_reader.GetNextRecord(bam_record)) {

        assert(bam_record.GetCigar().NumQueryConsumed() == bam_record.Length());

        auto transcript_genomic_regions_it = transcript_genomic_regions.emplace(bam_record.ChrName(bam_reader.Header()), SeqLib::GRC());
        cigarToGenomicRegions(&(transcript_genomic_regions_it.first->second), bam_record.GetCigar(), bam_record.Position());
    }

    bam_reader.Close();

    for (auto & chrom_regions: transcript_genomic_regions) {

        chrom_regions.second.MergeOverlappingIntervals();
        chrom_regions.second.CreateTreeMap();
    }

    return transcript_genomic_regions;
}

float fractionSoftClippedBases(const Cigar & cigar) {

    float num_soft_clipped_bases = 0;

    for (auto & field: cigar) {

        if (field.Type() == 'S') {

            num_soft_clipped_bases += field.Length();
        }
    }

    return num_soft_clipped_bases / cigar.NumQueryConsumed();
}

int main(int argc, char* argv[]) {

    if (argc != 3) {

        cout << "Usage: calc_read_transcript_overlap_stats <read_bam_name> <transcript_bam_name> > statistics.txt" << endl;
        return 1;
    }

    printScriptHeader(argc, argv);

    BamReader bam_reader;
    bam_reader.Open(argv[1]);
    assert(bam_reader.IsOpen());

    auto transcript_genomic_regions = createTranscriptGenomicRegions(argv[2]);

    uint32_t total_width = 0;

    for (auto & chrom_regions: transcript_genomic_regions) {

        total_width += chrom_regions.second.TotalWidth();
    }

    cerr << total_width << endl;

    unordered_map<string, uint32_t> overlap_stats;

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

            auto overlap_stats_it = overlap_stats.emplace("0\t0\t0", 0);
            overlap_stats_it.first->second++;

            continue;
        }

        assert(bam_record.GetCigar().NumQueryConsumed() == bam_record.Length());

        auto transcript_genomic_regions_it = transcript_genomic_regions.find(bam_record.ChrName(bam_reader.Header()));
        
        if (transcript_genomic_regions_it == transcript_genomic_regions.end()) {
        
            stringstream overlap_ss;
            overlap_ss << bam_record.MapQuality();
            overlap_ss << "\t" << fractionSoftClippedBases(bam_record.GetCigar());
            overlap_ss << "\t0";

            auto overlap_stats_it = overlap_stats.emplace(overlap_ss.str(), 0);
            overlap_stats_it.first->second++;

            continue;
        }

        auto read_cigar_genomic_regions = cigarToGenomicRegions(bam_record.GetCigar(), bam_record.Position());
        read_cigar_genomic_regions.CreateTreeMap();

        auto cigar_genomic_regions_intersection = transcript_genomic_regions_it->second.Intersection(read_cigar_genomic_regions, true);

        auto overlap = cigar_genomic_regions_intersection.TotalWidth() / static_cast<float>(read_cigar_genomic_regions.TotalWidth()); 

        stringstream overlap_ss;
        overlap_ss << bam_record.MapQuality();
        overlap_ss << "\t" << fractionSoftClippedBases(bam_record.GetCigar());
        overlap_ss << "\t" << overlap;

        auto overlap_stats_it = overlap_stats.emplace(overlap_ss.str(), 0);
        overlap_stats_it.first->second++;
    }

    bam_reader.Close();

    cout << "Count" << "\t" << "MapQ" << "\t" << "FracSoftClip" << "\t" << "Overlap" << endl;

    for (auto & stats: overlap_stats) {

        cout << stats.second << "\t" << stats.first << endl;
    }

	return 0;
}
