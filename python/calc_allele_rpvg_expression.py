
'''
calc_allele_rpvg_expression.py 
Calculates allele expression from rpvg expression
estimates. Note that the input variant (vcf) file 
need to annotated with transcript names 
(INFO:TRANSCIPTS tag).  
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

def parse_rpvg_expression(filename, prob_threshold):

	rpvg_exp = {}
	
	rpvg_file = gzip.open(filename, "rb")

	for line in rpvg_file:

		line_split = line.split("\t")
		assert(len(line_split) == 7)

		if line_split[0] == "Name":

			continue

		assert(not line_split[0] in rpvg_exp)

		if float(line_split[4]) >= prob_threshold:

			rpvg_exp[line_split[0]] = float(line_split[6])

	rpvg_file.close()

	return rpvg_exp


printScriptHeader()

if len(sys.argv) != 6:

	print("Usage: python calc_allele_rpvg_expression.py <variant_vcf_gz_name> <hst_input_gz_name> <rpvg_expression_gz_name> <rpvg_prob_threshold> <output_fasta_name>\n")
	sys.exit(1)


hst_info = parse_hst_info(sys.argv[2])
print(len(hst_info))

rpvg_exp = parse_rpvg_expression(sys.argv[3], float(sys.argv[4]))
print(len(rpvg_exp))

variant_file = gzip.open(sys.argv[1], "rb")
out_file = open(sys.argv[5], "w")

out_file.write("Chrom\tPos\tType\tLength\tExpression\n")

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

	alt_alleles = line_split[4].split(",")
	allele_exp = [0.0 for x in range(1 + len(alt_alleles))]

	transcripts = [x.split("=")[1].split(",") for x in line_split[7].split(";") if x.split("=")[0] == "TRANSCIPTS"]
	assert(len(transcripts) == 1)

	for transcript in transcripts[0]:

		if transcript in hst_info:

			for hst in hst_info[transcript]:

				gt = line_split[sample_names[hst[1][0]]]
				assert(not "/" in gt)

				allele = gt.split("|")[int(hst[1][1])]

				if allele != ".":

					cur_exp = 0.0

					if hst[0] in rpvg_exp:

						cur_exp = rpvg_exp[hst[0]]

					allele_exp[int(allele)] += cur_exp

	out_file.write(line_split[0] + "\t" + line_split[1] + "\tRef\t0\t" + str(allele_exp[0]) + "\n")

	for i in range(len(alt_alleles)):

		allele_type_length = getAlleleTypeLength(line_split[3], alt_alleles[i])
		assert(allele_type_length[0] != "Ref")

		out_file.write(line_split[0] + "\t" + line_split[1] + "\t" + allele_type_length[0] + "\t" + str(allele_type_length[1]) + "\t" + str(allele_exp[i + 1]) + "\n")

variant_file.close()
out_file.close()

print("Done")
