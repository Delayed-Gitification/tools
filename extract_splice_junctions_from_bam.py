import argparse
import pysam
from fuzzywuzzy import fuzz






def write_out(file, to_write):
	file.write(to_write)


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-b", "--bam", required=True)
	parser.add_argument("--min_intron_length", default=50, type=int)
	parser.add_argument("--chunk_size", default=1_000_000, type=int)
	parser.add_argument("-o", "--output", required=True)
	parser.add_argument("--early_stop", default=-1, type=int)
	args = parser.parse_args()


	with pysam.AlignmentFile(args.bam) as bam, open(args.output, 'w') as out:
		record_number = 0

		out.write("record,reference,secondary,junctions\n")
		skipped = 0
		for record in bam:
			try:
				if record.flag == 256:
					secondary = True
				else:
					secondary = False

				record_number += 1

				if record_number % 10_000 == 0:
					print(record_number)

				positions = record.get_reference_positions()

				junctions = []

				for i in range(len(positions) - 1):
					distance = positions[i + 1] - positions[i]
					if distance >= args.min_intron_length:
						junctions.append(str(positions[i]) + "-" + str(positions[i + 1]))

				to_write += ','.join([str(record_number), ref, str(secondary)])

				to_write += "," + ";".join(junctions) + "\n"


				if len(to_write) > args.chunk_size:
					write_out(out, to_write)
					to_write = ""

				if record_number > args.early_stop > 0:
					break
			except:
				error=1

		write_out(out, to_write)

	print(str(skipped) + " skipped")

if __name__ == '__main__':
	main()