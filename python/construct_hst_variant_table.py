
'''
construct_hst_variant_table.py 
'''

import sys
import os
import subprocess

import pickle
import gzip

from utils import *


def parse_hst_info(filename):

	hst_info = {}

	hst_file = gzip.open(filename, "rb")

	for line in hst_file:

		line_split = line.split("\t")
		assert(len(line_split) == 5)

		if line_split[0] == "Name":

			continue

		if not line_split[2] in hst_info:

			hst_info[line_split[2]] = []

		hst_info[line_split[2]].append((line_split[0], [line_split[4].split(",")[0].split("_")[i] for i in [2,4]]))

	hst_file.close()

	return hst_info


printScriptHeader()

if len(sys.argv) != 4:

	print("Usage: python construct_hst_variant_table.py <hst_input_name> <variant_vcf_name> <output_fasta_name>\n")
	sys.exit(1)


hst_info = parse_hst_info(sys.argv[1])
print(len(hst_info))

variant_file = gzip.open(sys.argv[2], "rb")
out_file = open(sys.argv[3], "w")

out_file.write("Chrom\tPos\tRef\tAllele\tHSTs\n")

sample_names = {}

for line in variant_file:

	line_split = line.split("\t")
	line_split[-1] = line_split[-1].strip()

	if line_split[0] == "#CHROM":

		assert(len(line_split) >= 10)

		for i in range(9, len(line_split)):

			sample_names[line_split[i]] = i

		continue

	if line_split[0][0] == "#":

		continue

	assert(len(line_split) >= 10)

	transcripts = [x.split("=")[1].split(",") for x in line_split[7].split(";") if x.split("=")[0] == "TRANSCIPTS"]
	assert(len(transcripts) == 1)

	alt_alleles = line_split[4].split(",")
	variant_hsts = [[] for x in range(1 + len(alt_alleles))]

	for transcript in transcripts[0]:

		for hst in hst_info[transcript]:

			gt = line_split[sample_names[hst[1][0]]]
			assert(not "/" in gt)

			allele = gt.split("|")[int(hst[1][1])]

			if allele != ".":

				variant_hsts[int(allele)].append(hst[0])

	out_file.write(line_split[0] + "\t" + line_split[1] + "\t" + line_split[3] + "\t" + line_split[3] + "\t" + ",".join(variant_hsts[0]) + "\n")

	for i in range(len(alt_alleles)):

		out_file.write(line_split[0] + "\t" + line_split[1] + "\t" + line_split[3] + "\t" + alt_alleles[i] + "\t" + ",".join(variant_hsts[i + 1]) + "\n")

variant_file.close()
out_file.close()

print("Done")
