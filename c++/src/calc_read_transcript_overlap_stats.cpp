
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

        auto chrom = bam_record.ChrName(bam_reader.Header());

        if (chrom.size() >= 4 && chrom.substr(0, 3) == "chr") {

            chrom = chrom.substr(3);
        }

        auto transcript_genomic_regions_it = transcript_genomic_regions.emplace(chrom, SeqLib::GRC());
        cigarToGenomicRegions(&(transcript_genomic_regions_it.first->second), bam_record.GetCigar(), bam_record.Position(), min_deletion_length);
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

    if (argc != 4) {

        cout << "Usage: calc_read_transcript_overlap_stats <read_bam_name> <transcript_bam_name> <min_deletion_length> > statistics.txt" << endl;
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

    BamRecord bam_record;

    cout << "Name" << "\t" << "MapQ" << "\t" << "FracSoftClip" << "\t" << "Overlap" << endl;

    uint32_t num_reads = 0;
    float sum_overlap = 0;

    while (bam_reader.GetNextRecord(bam_record)) { 

        if (bam_record.SecondaryFlag()) {

            continue;
        }

        num_reads++;

        float overlap = 0;
        float frac_soft_clip = 0;

        auto transcript_genomic_regions_it = transcript_genomic_regions.find(bam_record.ChrName(bam_reader.Header()));

        if (bam_record.MappedFlag()) { 

            assert(bam_record.GetCigar().NumQueryConsumed() == bam_record.Length());
            frac_soft_clip = fractionSoftClippedBases(bam_record.GetCigar());

            if (transcript_genomic_regions_it != transcript_genomic_regions.end()) {

                auto read_cigar_genomic_regions = cigarToGenomicRegions(bam_record.GetCigar(), bam_record.Position(), 0);
                read_cigar_genomic_regions.CreateTreeMap();

                auto cigar_genomic_regions_intersection = transcript_genomic_regions_it->second.Intersection(read_cigar_genomic_regions, true);

                overlap = cigar_genomic_regions_intersection.TotalWidth() / static_cast<float>(bam_record.Length()); 
            }
        }

        cout << bam_record.Qname();
        cout << "\t" << bam_record.MapQuality();
        cout << "\t" << frac_soft_clip;
        cout << "\t" << overlap;
        cout << endl;

        sum_overlap += overlap;

        if (num_reads % 1000000 == 0) {

            cerr << "Number of analysed reads: " << num_reads << endl;
        }        
    }

    bam_reader.Close();

    cerr << "\nTotal number of analysed reads: " << num_reads << endl;
    cerr << "Average overlap: " << sum_overlap/num_reads << endl;

	return 0;
}
