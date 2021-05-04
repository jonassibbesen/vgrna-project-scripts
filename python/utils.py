
'''
utils.py
'''

import os
import sys
import subprocess

def printScriptHeader():

	script_dir = os.path.dirname(os.path.abspath(__file__))

	git_commit_hash = subprocess.check_output(["git", "-C", script_dir, "rev-parse", "HEAD"]).strip()
	git_branch = subprocess.check_output(["git", "-C", script_dir, "rev-parse", "--abbrev-ref", "HEAD"]).strip()

	sys.stderr.write(str(git_branch) + " " + str(git_commit_hash) + "\n")
	sys.stderr.write(" ".join(sys.argv) + "\n\n")

def getAlleleTypeLength(ref_allele, alt_allele):

    assert(len(ref_allele) > 0)
    assert(len(alt_allele) > 0)

    if (ref_allele == alt_allele):

        return ("Ref", 0)

    right_trim_len = 0

    for i in range(1, min(len(ref_allele), len(alt_allele)) + 1):

        if ref_allele[-1 * i] == alt_allele[-1 * i]:

            right_trim_len += 1
        
        else:

        	break
        
    if right_trim_len > 0:

    	ref_allele = ref_allele[:(-1 * right_trim_len)]
    	alt_allele = alt_allele[:(-1 * right_trim_len)]

    left_trim_len = 0

    for i in range(0, min(len(ref_allele), len(alt_allele))):

        if ref_allele[i] == alt_allele[i]:

            left_trim_len += 1
        
        else:

        	break

    if left_trim_len > 0:

    	ref_allele = ref_allele[left_trim_len:]
    	alt_allele = alt_allele[left_trim_len:]

    if len(ref_allele) == len(alt_allele) and len(ref_allele) == 1:

       	return ("SNV", 1)

    elif len(ref_allele) == 0:

       	assert(len(alt_allele) > 0)
        return ("Ins", len(alt_allele))

    elif len(alt_allele) == 0:

       	assert(len(ref_allele) > 0)
        return ("Del", -1 * len(ref_allele))    

    else:

       	assert(len(ref_allele) > 0)
       	assert(len(alt_allele) > 0)
        return ("Com", len(alt_allele) - len(ref_allele))
    
