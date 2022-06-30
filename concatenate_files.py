import glob
import sys
import os
import argparse

# find folder names 
directory = os.getcwd()
dirs = [x[0] for x in os.walk(directory)]

illumina_barcodes = set([x.split("_L00")[0] for x in dirs if "_L00" in x])
print(illumina_barcodes)

for ib in illumina_barcodes:
	# get all fastq names
	fq_names = [x.split("/")[-1] for x in glob.glob(ib + '_L001_/' + "*fastq.gz")]
	
	for fq in fq_names:
		ib_name = ib.split("/")[-1]
		if not os.path.isdir("combined_" + ib_name):
			os.makedirs("combined_" + ib_name)
		command = "cat " + ib + "_L001_/" + fq + " " + ib + "_L002_/" + fq + " > combined_" + ib_name + "/"+ fq 
		print(command)
		os.system('sbatch --wrap="' + command + '"')


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-d", "--directory", help="directory that the files are in. Can be a list", nargs = "+", required = True)
	args = parser.parse_args()

	



if __name__=="__main__":
	main()