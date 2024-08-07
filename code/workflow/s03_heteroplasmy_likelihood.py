from Bio import SeqIO
import math
import statistics
import sys
import annotate

# ref and read must have same length
def diff(ref, ref_pos, read, read_pos):
	return [ (i+read_pos, read[i], i+ref_pos, ref[i]) for (i,base) in enumerate(ref) if base!=read[i] ]
	

def reverse_complement(s):
	rc = list(s)
	left, right = 0, len(rc)-1
	complement = {'A':'T','T':'A','C':'G','G':'C'}
	while left < right:
		rc[left], rc[right] = complement.get(rc[right],rc[right]), complement.get(rc[left],rc[left])
		left, right = left+1, right-1
	return ''.join(rc)

# -----------------------------------------------------------
def collect_info_alleles(cigar, ref, read, reference_pos, phred, profile={}):
	break_points = [ i[0] for i in enumerate(cigar) if i[1].isalpha() ]
	start, ref_pos, read_pos = 0, 0, 0
	count = 0
	# print(break_points)
	for b in break_points:
		op, length = cigar[b], int(cigar[start:b])
		# print(">", op, length) 
		# print("ref_pos = ", ref_pos, "read_pos =", read_pos)
		start = int(b)+1
		if op == 'S':
			read_pos += length
		elif op == 'M':
			ref_substr = ref[ref_pos:ref_pos+length]
			read_substr = read[read_pos:read_pos+length]
			differences = diff(ref_substr, ref_pos, read_substr, read_pos)
			count += len(differences)
			read_pos += length
			ref_pos += length
			for d in differences:
				if reference_pos + d[2] not in profile:
					profile[reference_pos+d[2]] = {'ref':d[3], 'err': {}}
					profile[reference_pos+d[2]]['err'][d[3]] = []
				profile[reference_pos+d[2]]['err'].setdefault(d[1], [])
				profile[reference_pos+d[2]]['err'][d[1]].append(10**(-(ord(phred[d[0]])-33)/10))

				if d[3] != profile[reference_pos+d[2]]['ref']:
					print("Something is wrong", d[3], profile[reference_pos+d[2]])
					raise Exception("QUIT")
		elif op == 'I':
			read_substr = read[read_pos:read_pos+length]

			if reference_pos + ref_pos -1 not in profile:
				profile[reference_pos+ref_pos-1] = {'ref':ref[ref_pos-1], 'err':{}}
				profile[reference_pos+ref_pos-1]['err'][ref[ref_pos-1]] = []
			

			profile[reference_pos+ref_pos-1]['err'].setdefault('I',[])

			e = 0
			for i in range(length):
				e = e + 10**(-(ord(phred[read_pos+i])-33)/10)
			profile[reference_pos+ref_pos-1]['err']['I'].append(float(e/length))
			
			# profile[reference_pos+ref_pos]['err']['I'].append(10**(-(ord(phred[ref_pos])-33)/10))
			
			read_pos += length 			

		elif op == 'D':
			ref_substr = ref[ref_pos:ref_pos+length]
			
			for i in range(length):
				
				if reference_pos+ref_pos not in profile:
					profile[reference_pos+ref_pos] = {'ref':ref_substr[i], 'err': {}}
					profile[reference_pos+ref_pos]['err'][ref_substr[i]] = []


				profile[reference_pos+ref_pos]['err'].setdefault('D', [])
				profile[reference_pos+ref_pos]['err']['D'].append(0)
				# profile[reference_pos+ref_pos]['err']['D'].append(10**(-(ord(phred[ref_pos])-33)/10))
				ref_pos += 1

		elif op == 'H':
			pass
	# print(count, "potential SNPs found.")

def extract_ref_seq(sequence, cigar, ref_pos, sam_flag):
	break_points = [ i[0] for i in enumerate(cigar) if i[1].isalpha() ]
	L, start = 0, 0
	for b in break_points:
		op, length = cigar[b], int(cigar[start:b])
		start = int(b)+1
		if op in ['M','D']:
			L += length
	# if (sam_flag & 82) == 82 or (sam_flag & 162) == 162:
	# 	return reverse_complement(str(sequence[ref_pos : ref_pos+L]))
	return str(sequence[ref_pos : ref_pos+L])


# -----------------------------------------------------------
def collect_info_ref(cigar, ref, read, reference_pos, phred, profile):
	break_points = [ i[0] for i in enumerate(cigar) if i[1].isalpha() ]
	start, ref_pos, read_pos = 0, 0, 0
	for b in break_points:
		op, length = cigar[b], int(cigar[start:b])
		start = int(b)+1
		if op == 'S':
			read_pos += length
		elif op == 'M':
			ref_substr = ref[ref_pos:ref_pos+length]
			read_substr = read[read_pos:read_pos+length]
			for i,c in enumerate(read_substr):
				cur_pos = reference_pos + ref_pos + i
				if cur_pos in profile:
					if c == profile[cur_pos]['ref']:
						profile[cur_pos]['err'][c].append(10**(-(ord(phred[read_pos+i])-33)/10))
			read_pos += length
			ref_pos += length
		elif op == 'I':
			read_pos += length
		elif op == 'D':
			ref_pos += length
		elif op == 'H':
			pass


# -----------------------------------------------------------
# heteroplasmy, Ye 2014
# -----------------------------------------------------------
def log_likelihood(major_allele, major_f, errors):
	value = 0
	for allele in errors:
		if allele == major_allele:
			for err in errors[allele]:
				v = (1-major_f)*err + major_f*(1-err)
				if v>0:
					value += math.log10(v)
		else:
			for err in errors[allele]:
				v = (1-major_f)*(1-err) + major_f*err
				if v>0:
					value += math.log10(v)
	return value


# -----------------------------------------------------------
def analyze_positions(profile, ANNOTATION_FILE, out_csv):
	with open(out_csv, 'w') as f:
		# print("Pos,Ref,A,C,G,T,D,I,Total,Ea,Ec,Eg,Et,Ed,Ei,Score,GeneProduct")
		f.write("Pos,Ref,A,C,G,T,D,I,Total,Percentage,Ea,Ec,Eg,Et,Ed,Ei,Score,GeneProduct\n")
		for pos, p in sorted(profile.items()):
			# biopython is off by one position
			pos = pos + 1 

			# print('%s,%s' % (pos,p['ref']), end=",")
			f.write(str(pos) + ',' + str(p['ref']) + ',')

			# determine major allele freq
			major_allele, major_f, count = None, None, 0
			for a in p['err']:
				if a != 'I':
					if (major_f is None) or (major_f < len(p['err'][a])):
						major_allele, major_f = a, len(p['err'][a])
					count += len(p['err'][a])
			major_f = major_f / count

			for base in ['A','C','G','T','D','I']:
				if base not in p['err']:
					# print(0, end=',')
					f.write(str(0) + ',')
				else:
					# print(len(p['err'][base]), end=',')
					f.write(str(len(p['err'][base])) + ',')
					
			# print(count,end=',')
			f.write(str(count) + ',')

			remaining_sum = 1 - major_f
			f.write(str(remaining_sum) + ',')
			
			for base in ['A','C','G','T','D','I']:
				if (base not in p['err']) or (p['err'][base]==[]):
					# print(0, end=',')
					f.write(str(0) + ',')
				else:			
					# print(statistics.mean(p['err'][base]), end=',')
					f.write(str(statistics.mean(p['err'][base])) + ',')
					

			# print(log_likelihood(major_allele, major_f, p['err']) - log_likelihood(major_allele, 1, p['err']), end=',')
			f.write(str(log_likelihood(major_allele, major_f, p['err']) - log_likelihood(major_allele, 1, p['err'])) + ',')

			annotations = annotate.get_entries(ANNOTATION_FILE)
			loc = annotate.search(annotations, pos)
			if loc == -1:
				# print("non-coding")
				f.write("non-coding\n")
			else:
				gene_product = annotations[loc]
				# print(gene_product[2])
				f.write(gene_product[2] + '\n')

# -----------------------------------------------------------
# Detect heteroplasmy candidates
def parse_sam_file(filename, ref_seq, annotation_file, out_csv):
	profile = {}
	count = 0
	# 1. Detect heteroplasmy sites.  Update profile of alleles that are different from ref.
	with open(filename, "rU") as file:
		for line in file:
			if line[0] != '@':
				count += 1
				items = line.split('\t')
				sam_flag, ref_pos, cigar, read, phred = int(items[1]), int(items[3])-1, items[5], items[9], items[10]
				seq = extract_ref_seq(ref_seq, cigar, ref_pos, sam_flag)
				collect_info_alleles(cigar, seq, read, ref_pos, phred, profile)

	# 2. Update profiles of alleles that are the same as ref.
	# Note: have to do this because reads with bases matching reference bases might be missed
	# if we go through only one pass (step 1)
	with open(filename, "rU") as file:
		for line in file:
			if line[0] != '@':
				items = line.split('\t')
				sam_flag, ref_pos, cigar, read, phred = int(items[1]), int(items[3])-1, items[5], items[9], items[10]
				seq = extract_ref_seq(ref_seq, cigar, ref_pos, sam_flag)
				collect_info_ref(cigar, seq, read, ref_pos, phred, profile)

	analyze_positions(profile, annotation_file, out_csv)

	# print(profile[63754])

def process(**params):
	ref = params['ref']
	record = SeqIO.read(ref, "fasta")
	alignment_output = params['out_filtered_sam']
	annotation_file = params['annotation']
	out_csv = params['out_csv']
	parse_sam_file(alignment_output, record.seq, annotation_file, out_csv)


if __name__ == '__main__':
	if len(sys.argv) != 5:
		print("Usage:", sys.argv[0] + "  reference.fasta  alignment_output.sam  annotation_file output_csv")
	else:
		# record = SeqIO.read(sys.argv[1], "fasta")
		params = {
			'ref' : sys.argv[1],
			'out_filtered_sam' : sys.argv[2],
			'annotation' : sys.argv[3],
			'out_csv' : sys.argv[4],
		}
		
		process(**params)
		# parse_sam_file(alignment_output, record.seq, annotation_file, out_csv)
