import argparse
import pysam

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

	output_d = {}

	with pysam.AlignmentFile(args.bam) as bam:
		record_number = 0

		skipped = 0
		try:
			for record in bam:
				if record.flag == 256:
					secondary = True
				else:
					secondary = False

				record_number += 1

				if record_number % 10_000 == 0:
					print(record_number)

				positions = record.get_reference_positions()

				mapping_quality = record.mapping_quality

				ref = record.reference_name
				junctions = []

				for i in range(len(positions) - 1):
					distance = positions[i + 1] - positions[i]
					if distance >= args.min_intron_length:
						junctions.append(str(positions[i]) + "-" + str(positions[i + 1]))

				to_write = ','.join([ref, str(mapping_quality), str(secondary)])

				to_write += "," + ";".join(junctions)

				if to_write in output_d.keys():
					output_d[to_write] += 1
				else:
					output_d[to_write] = 1


				if record_number > args.early_stop > 0:
					break
		except:
			skipped+=1

		with open(args.output, 'w') as out:
			out.write("reference,mapping_quality,secondary,junctions,number_of_reads\n")
			for key, value in output_d.items():
				out.write(key + "," + str(value) + "\n")

	print(str(skipped) + " skipped")

if __name__ == '__main__':
	main()