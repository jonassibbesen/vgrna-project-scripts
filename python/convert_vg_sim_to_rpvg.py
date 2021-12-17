
'''
convert_vg_sim_to_rpvg.py 
Convert vg sim path count output to rpvg 
expression output.  
'''

import sys
import os
import subprocess

import pickle
import gzip

from Bio.Seq import Seq
from Bio import SeqIO

from utils import *


def parse_path_counts(filename):

	path_counts = {}
	
	vg_sim_file = gzip.open(filename, "rb")

	for line in vg_sim_file:

		line_split = line.split("\t")
		assert(len(line_split) == 4)

		if line_split[0] == "read":

			continue

		if line_split[1] in path_counts: 

			path_counts[line_split[1]] += 1

		else:

			path_counts[line_split[1]] = 1

	vg_sim_file.close()

	return path_counts


printScriptHeader()

if len(sys.argv) != 3:

	print("Usage: python convert_vg_sim_to_rpvg.py <vg_sim_gz_name> <output_file_name>\n")
	sys.exit(1)


path_counts = parse_path_counts(sys.argv[1])
print(len(path_counts))

out_file = open(sys.argv[2], "w")
out_file.write("Name\tClusterID\tLength\tEffectiveLength\tHaplotypeProbability\tReadCount\tTPM\n")	

for path, count in path_counts.items():

	out_file.write(path + "\t0\t0\t0\t0\t" + str(count) + "\t" + str(count) + "\n")

out_file.close()

print("Done")
