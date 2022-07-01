import argparse
import dnaio
from rapidfuzz import fuzz
from os.path import exists
import gzip


def rev_c(seq):
	"""
	simple function that reverse complements a given sequence
	"""
	tab = str.maketrans("ACTGN", "TGACN")
	# first reverse the sequence
	seq = seq[::-1]
	# and then complement
	seq = seq.translate(tab)
	return seq


def fast_fuzz(s1, s2):
	if (s1 in s2) or (s2 in s1):
		return 100
	else:
		return fuzz.partial_ratio(s1, s2)


def find_barcode(seq, barcodes, min_score, max_ambiguity, stored_results):
	scores = []
	names = []

	if seq in stored_results.keys():
		return stored_results[seq]

	for name, bc in barcodes.items():
		score = fast_fuzz(seq, bc)
		scores.append(score)
		names.append(name)

	# Reject if similar to >1 barcode (ambiguous)
	if len([a for a in scores if a > max_ambiguity]) > 1:
		stored_results[seq] = -1

	elif max(scores) >= min_score:
		stored_results[seq] = names[scores.index(max(scores))]

	else:
		stored_results[seq] = -1
	
	return stored_results[seq], stored_results


def initialise_d(forward_primers, reverse_primers):
	d = {}
	for i in list(forward_primers.keys()) + [-1]:
		for j in list(reverse_primers.keys()) + [-1]:
			d[str(i) + "_" + str(j)] = ""
	return d


def write_out_d(d, output):
	for key, value in d.items():

		filename = output + "_" + key + ".fastq.gz"

		if exists(filename):
			mode = "ab"
		else:
			mode = "wb"

		with gzip.open(filename, mode) as file:
			file.write(value.encode())

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--fastq", required=True)
	parser.add_argument("-p", "--primers", required=True, help="This should be a headerless, three column csv with columns 1=name,2=barcode,3=f or r")
	parser.add_argument("-s", "--min_score", type=float, default=92, help="When matching barcodes")
	parser.add_argument("--max_ambiguity", type=float, default=76, help="If another barcode has this score or higher, then ignore this read because it's ambiguous")
	parser.add_argument("-l", "--length", type=int, default=35, help="length to look at start and end of reads "
	                                                                 "for barcode")
	parser.add_argument("--ignore_rc", action="store_true", default=False, help="by default it looks at the read and its reverse complement. Use this command to only look at the forward read")
	parser.add_argument("-o", "--output", required=True)
	args = parser.parse_args()

	assert args.min_score >= args.max_ambiguity

	# Read in primers
	forward_primers = {}
	reverse_primers = {}
	with open(args.primers) as file:
		for i, line in enumerate(file):
			if i == 0:
				continue

			line2 = line.rstrip().split(",")

			this_name = line2[0]
			this_bc = line2[1].upper()
			f_or_r = line2[2]

			if f_or_r == "R":
				this_bc = rev_c(this_bc)

			#bc_number = int(line2[2].replace("pTwist_nano_F", "").replace("pTwist_nano_R", ""))

			# if f_or_r == "F" and bc_number > 8:
			# 	continue

			if f_or_r == "F":
				forward_primers[this_name] = this_bc
			else:
				reverse_primers[this_name] = this_bc

	to_write_d = initialise_d(forward_primers, reverse_primers)

	print(forward_primers)
	print(reverse_primers)
	print(to_write_d)

	stored_results1 = {}
	stored_results2 = {}

	with dnaio.open(args.fastq) as file:
		since_written = 0
		for i, record in enumerate(file):
			forward_results = []
			reverse_results = []
			for rc in [True, False]:
				if rc:
					if args.ignore_rc:
						continue

					seq = rev_c(str(record.sequence))
				
				else:
					seq = str(record.sequence)

				s1 = seq[0:args.length]
				s2 = seq[-args.length:]

				bc1, stored_results1 = find_barcode(s1, forward_primers, args.min_score, args.max_ambiguity, stored_results1)
				bc2, stored_results2 = find_barcode(s2, reverse_primers, args.min_score, args.max_ambiguity, stored_results2)

				if rc:
					reverse_results = [bc1, bc2]
				else:
					forward_results = [bc1, bc2]

			if args.ignore_rc:
				if -1 in forward_results:
					continue
				bc1 = forward_results[0]
				bc2 = forward_results[1]
			else:
				if forward_results == [-1, -1] and -1 not in reverse_results:
					# reverse results are good
					bc1 = reverse_results[0]
					bc2 = reverse_results[1]
				elif -1 not in forward_results and reverse_results == [-1, -1]:
					bc1 = forward_results[0]
					bc2 = forward_results[1]



			to_write_d[str(bc1) + "_" + str(bc2)] += "\n".join(["@"+record.name, record.sequence, "+", record.qualities]) + "\n"
			since_written += 1

			if since_written == 10_000:
				write_out_d(to_write_d, args.output)
				to_write_d = initialise_d(forward_primers, reverse_primers)
				since_written = 0
				print(i)

		write_out_d(to_write_d, args.output)



if __name__ == "__main__":
	main()




