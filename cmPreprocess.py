import argparse
import gzip
import time
import os
start = time.perf_counter()

parser = argparse.ArgumentParser('cmPreprocess turns scRNA fastq files (currently supports 10x v2 and 10x v3) into fasta \
files. Each fasta file is a cell, defined by the cell barcode')
parser.add_argument('read1', help='The name of the first fastq file')
parser.add_argument('read2', help='The name of the second fastq file')
parser.add_argument('output', help='Output folder')
parser.add_argument('--num-bases', type=int, default=16, help='The number of base pairs to extract from the first read to pull the cell barcode. 10xv2 is 16 (default) and 10xv3 is 18.')
parser.add_argument('--min-read-num', type=int, default=1000, help='The minimum number of reads required to retain a cell (default 1000)')
parser.add_argument('--chunk-size', type=int, default=200, help='How many reads to store in memory until write to disk. Faster to be the same \
as min-read-num, but will use more memory (default 200).')
parser.add_argument('--compress', default = False, action='store_true', help = 'Compress (gzip) the output fasta files in your output directory (default no)')
parser.add_argument('--no-compress', dest='compress', action='store_false',  help = 'Do not compress (gzip) the output fasta files in your output directory (default)')

# parse the command line arguments
args = parser.parse_args()

# set up dict to count number of reads per cell barcode
# and another to hold the sequence
cell_size_dict = dict() 
cell_fq_dict = dict()
keeper_bc = []

# make output dir
if not os.path.exists(args.output):
	os.makedirs(args.output)

# read in the paired fastq files
with gzip.open(args.read1, 'rt') as read1, gzip.open(args.read2, 'rt') as read2:
	lines1 = []
	lines2 = []
	# iterate over the lines in the files
	for line1, line2 in zip(read1, read2):
		lines1.append(line1)
		lines2.append(line2)
		if len(lines1) == 4 and len(lines2) == 4:
			# extract the specified number of base pairs from the first read
			seq = lines1[1][:args.num_bases]
			# if it is a new barcode then add to dicts
			if seq not in cell_size_dict:
				cell_size_dict[seq] = 1
				cell_fq_dict[seq] = lines2[:2]
			# if the number of seq reaches the min read read num, then
			# mark as in the keeper_bc list as a kept barcode	
			if cell_size_dict[seq] == args.min_read_num:
				keeper_bc.append(seq)
			# if less than n seq and below the chunk size, 
			# then append info and increment the counter
			if seq in cell_size_dict and \
				cell_size_dict[seq] <= args.chunk_size:
					cell_size_dict[seq] += 1
					cell_fq_dict[seq].extend(lines2[:2])
			# if you hit chunk size and it is a keeper seq,
			# then write collected info to file and reset counter
			elif cell_size_dict[seq] > args.chunk_size and seq in keeper_bc:
				cell_size_dict[seq] = 1
				f = open(args.output + '/' + seq + '.fa', 'a')
				cell_fq_dict[seq].extend(lines2[:2])
				f.write(''.join(cell_fq_dict[seq]))
				cell_fq_dict[seq] = []
				f.close()
			# if you hit chunk size but aren't certain whether this is a 
			# keeper seq, then keep writing to dict
			elif cell_size_dict[seq] > args.chunk_size:
				cell_size_dict[seq] += 1
				cell_fq_dict[seq].extend(lines2[:2])
			else:
				print('Whoops: ' + cell_size_dict[seq])
			# clear the fastq info
			lines1 = []
			lines2 = []


# one more loop to print out remainder entries
for seq in set(keeper_bc):
	f = open(args.output + '/' + seq + '.fa', 'a')
	f.write(''.join(cell_fq_dict[seq]))
	f.close()

# optional compress
if args.compress:
	import gzip
	import shutil
	for seq in set(keeper_bc):
		with open(args.output + '/' + seq + '.fa', 'rb') as f_in:
			with gzip.open(args.output + '/' + seq + '.fa.gz', 'wb') as f_out:
				shutil.copyfileobj(f_in, f_out)
		os.remove(args.output + '/' + seq + '.fa')
	

end = time.perf_counter()
elapsed = end - start
print('Execution time: {:.3f} seconds'.format(elapsed))
