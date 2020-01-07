import sys
import os
import subprocess

def write_average_exp(hap1_exp, hap2_exp, exp_out_file):

	assert(len(hap1_exp) != 0 or len(hap2_exp) != 0)

	average_tpm = 0

	if len(hap1_exp) != 0 and len(hap2_exp) != 0:

		average_tpm = (float(hap1_exp[5]) + float(hap2_exp[5])) / 2

	if len(hap1_exp) != 0:

		hap1_exp[4:] = ["0"] * 4
		hap1_exp[5] = str(average_tpm)

		exp_out_file.write("\t".join(hap1_exp) + "\n")

	if len(hap2_exp) != 0:
	
		hap2_exp[4:] = ["0"] * 4
		hap2_exp[5] = str(average_tpm)

		exp_out_file.write("\t".join(hap2_exp) + "\n")


script_dir = os.path.dirname(os.path.abspath(__file__))

git_commit_hash = subprocess.check_output(["git", "-C", script_dir, "rev-parse", "HEAD"]).strip()
git_branch = subprocess.check_output(["git", "-C", script_dir, "rev-parse", "--abbrev-ref", "HEAD"]).strip()

print(git_branch + " " + git_commit_hash)
print(sys.argv)
print

if len(sys.argv) != 3:

	print("Usage: python average_diploid_rsem_expression.py <input_name> <output_name>\n")
	sys.exit(1)

exp_in_file = open(sys.argv[1], "r")
exp_out_file = open(sys.argv[2], "w")

hap1_exp = []
hap2_exp = []

prev_transcript_id = ""

for line in exp_in_file:

	line_split = line.strip().split("\t")
	
	if line_split[0] == "transcript_id":

		exp_out_file.write(line)
		continue

	# Assumes haplotype-specefic transcriptcripts are next to each other
	assert(line_split[0][-2] == "_" or hap_id == "_")

	hap_id = line_split[0][-1]
	cur_transcript_id = line_split[0][:-2]

	if cur_transcript_id != prev_transcript_id and prev_transcript_id != "":

		write_average_exp(hap1_exp, hap2_exp, exp_out_file)

		hap1_exp = []
		hap2_exp = []

	if hap_id == "1":

		assert(len(hap1_exp) == 0)
		hap1_exp = line_split

	else:

		assert(hap_id == "2")
		assert(len(hap2_exp) == 0)
		hap2_exp = line_split

	prev_transcript_id = cur_transcript_id

write_average_exp(hap1_exp, hap2_exp, exp_out_file)

exp_in_file.close()
exp_out_file.close()

print("Done")
