
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
