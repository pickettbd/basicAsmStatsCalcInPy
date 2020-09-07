
def extractNlens(seq):

	n_lens = []

	n_len = 0
	i = 0
	while i < len(seq):
		
		while i < len(seq) and seq[i] == 'N':
			n_len += 1
			i += 1

		if n_len > 0:
			n_lens.append(n_len)
			n_len = 0

		i += 1

	return n_lens

def parseInputFasta(ifn):
	
	seq_lens = []
	n_lens = []
	seqs_with_ns_count = 0

	with open(ifn, 'r') as ifd:
		
		line = ifd.readline()

		while line != '':
			header = line
			seq = ifd.readline().rstrip('\n').upper()

			line = ifd.readline()

			while line != '' and line[0] != '>':
				seq += line.rstrip('\n').upper()
				line = ifd.readline()

			n_len = extractNlens(seq)

			if len(n_len) > 0:
				seqs_with_ns_count += 1
				n_lens.extend(n_len)
			
			seq_lens.append(len(seq) - sum(n_len))
	
	return sorted(seq_lens)[::-1], sorted(n_lens)[::-1], seqs_with_ns_count

def determineGcutoff(seq_total, ngx_cutoffs):
	
	for i in range(len(ngx_cutoffs) - 1, -1, -1):
		if seq_total >= ngx_cutoffs[i]:
			return i + 1
	
	return 0


if __name__ == "__main__":
	import sys
	import statistics as stats

	if len(sys.argv) != 4:
		sys.stderr.write("ifn, ofn, genome_size\n\n")
		sys.exit()

	ifn = sys.argv[1]
	ofn = sys.argv[2]
	genome_size = int(sys.argv[3])

	seq_lens, n_lens, seqs_with_ns_count = parseInputFasta(ifn)
	
	seq_total = sum(seq_lens)

	percents = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 ]
	nx_cutoffs = [ int( p * seq_total ) for p in percents ]
	ngx_cutoffs = [ int( p * genome_size ) for p in percents ]
	
	g_cutoff = determineGcutoff(seq_total, ngx_cutoffs)
	
	nx = []
	ngx = []
	lx = []
	lgx = []
	nl_nucs = []
	nlg_nucs = []

	nx_index = 0
	ngx_index = 0
	sl_index = 0
	sl_count = 0

	while sl_index < len(seq_lens) and (nx_index < 10 or ngx_index < 10):
		
		seq_len = seq_lens[sl_index]
		sl_count += seq_len
		
		while nx_index < 10 and sl_count >= nx_cutoffs[nx_index]:
			nx_index += 1
			nx.append(seq_len)
			lx.append(sl_index + 1)
			nl_nucs.append(sl_count)

		while ngx_index < 10 and sl_count >= ngx_cutoffs[ngx_index]:
			ngx_index += 1
			ngx.append(seq_len)
			lgx.append(sl_index + 1)
			nlg_nucs.append(sl_count)
			
		sl_index += 1

	if len(nx) != 10 or len(ngx) != g_cutoff or len(lx) != 10 or len(lgx) != g_cutoff or len(nl_nucs) != 10 or len(nlg_nucs) != g_cutoff:
		sys.stderr.write("Houston, we have a problem..\n")
		sys.stderr.write("nx len: " + str(len(nx)) + '\n')
		sys.stderr.write("ngx len: " + str(len(ngx)) + '\n')
		sys.stderr.write("lx len: " + str(len(lx)) + '\n')
		sys.stderr.write("lgx len: " + str(len(lgx)) + '\n')
		sys.stderr.write("nl_nucs len: " + str(len(nl_nucs)) + '\n')
		sys.stderr.write("nlg_nucs len: " + str(len(nlg_nucs)) + '\n')
		sys.exit(1)
	
	c1 = ["x", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"]
	c2 = ["N"] + list(map(str, nx))
	c3 = ["L"] + list(map(str, lx))
	c4 = ["N/L nt"] + list(map(str, nl_nucs))
	c5 = ["NG"] + list(map(str, ngx)) + ["NA"] * (10 - len(ngx))
	c6 = ["LG"] + list(map(str, lgx)) + ["NA"] * (10 - len(lgx))
	c7 = ["NG/LG nt"] + list(map(str, nlg_nucs)) + ["NA"] * (10 - len(nlg_nucs))

	c1w = len(sorted(c1, key=lambda x: len(x))[-1])
	c2w = len(sorted(c2, key=lambda x: len(x))[-1])
	c3w = len(sorted(c3, key=lambda x: len(x))[-1])
	c4w = len(sorted(c4, key=lambda x: len(x))[-1])
	c5w = len(sorted(c5, key=lambda x: len(x))[-1])
	c6w = len(sorted(c6, key=lambda x: len(x))[-1])
	c7w = len(sorted(c7, key=lambda x: len(x))[-1])

	seq_lens_mean = stats.mean(seq_lens) if len(seq_lens) > 0 else 0
	seq_lens_stdev = stats.stdev(seq_lens, xbar=seq_lens_mean) if len(seq_lens) > 1 else 0
	n_lens_mean = stats.mean(n_lens) if len(n_lens) > 0 else 0
	n_lens_stdev = stats.stdev(n_lens, xbar=n_lens_mean) if len(n_lens) > 1 else 0

	with open(ofn, 'w') as ofd:
		ofd.write("                       Assumed Genome Size: " + str(genome_size) + '\n')
		ofd.write("                 Total Number of Sequences: " + str(len(seq_lens)) + '\n')
		ofd.write("Total Number of Nucleotides (excluding Ns): " + str(seq_total) + '\n')
		ofd.write("         Average (st. dev) Sequence Length: " + str(seq_lens_mean) + " (" + str(seq_lens_stdev) + ")\n")
		ofd.write("             Number of Sequences with Gaps: " + str(seqs_with_ns_count) + '\n')
		ofd.write("                            Number of Gaps: " + str(len(n_lens)) + '\n')
		ofd.write("                        Total Size of Gaps: " + str(sum(n_lens)) + '\n')
		ofd.write("              Average (st. dev) Gap Length: " + str(n_lens_mean) + " (" + str(n_lens_stdev) + ")\n")
		ofd.write('\n')

		for i in range(0, 11, 1):
			ofd.write(" " * (c1w - len(c1[i])))
			ofd.write(c1[i] + "  ")
			ofd.write(" " * (c2w - len(c2[i])))
			ofd.write(c2[i] + "  ")
			ofd.write(" " * (c3w - len(c3[i])))
			ofd.write(c3[i] + "  ")
			ofd.write(" " * (c4w - len(c4[i])))
			ofd.write(c4[i] + "  ")
			ofd.write(" " * (c5w - len(c5[i])))
			ofd.write(c5[i] + "  ")
			ofd.write(" " * (c6w - len(c6[i])))
			ofd.write(c6[i] + "  ")
			ofd.write(" " * (c7w - len(c7[i])))
			ofd.write(c7[i] + '\n')
			
