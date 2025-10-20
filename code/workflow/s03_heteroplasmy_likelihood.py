from Bio import SeqIO
import math
import statistics
from scipy.stats import chi2
from statsmodels.stats.multitest import fdrcorrection
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

def calculate_required_coverage(p_min, power_levels=[0.95, 0.99, 0.999], k=1):
	"""
	计算在给定统计功效下，至少观察到k次次要等位基因所需的最低覆盖度N。
	该函数目前实现了k=1时的计算。
	公式: N >= log(1 - power) / log(1 - p_min)
	
	Args:
		p_min (float): 用于计算的次要等位基因频率档位。
		power_levels (list): 需要计算的统计功效水平列表。
		k (int): 需要观察到的最少次数。

	Returns:
		dict: 包含不同功效水平下所需N值的字典。
	"""
	required_n = {}
	# 如果频率为0或负数，无法检测，返回'NA'
	if p_min <= 0:
		return {'N_req_95': 'NA', 'N_req_99': 'NA', 'N_req_999': 'NA'}

	# 此处计算仅适用于 k=1 的情况
	log_1_minus_p = math.log(1 - p_min)

	# 95% 功效
	n_95 = math.log(1 - 0.95) / log_1_minus_p
	required_n['N_req_95'] = math.ceil(n_95)

	# 99% 功效
	n_99 = math.log(1 - 0.99) / log_1_minus_p
	required_n['N_req_99'] = math.ceil(n_99)

	# 99.9% 功效
	n_999 = math.log(1 - 0.999) / log_1_minus_p
	required_n['N_req_999'] = math.ceil(n_999)
	
	return required_n

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
	"""
	Analyzes genomic positions to call heteroplasmy, incorporating FDR correction.
	This is a two-pass process: first collect all results, then calculate FDR and write.
	"""
	
	# --- PASS 1: Calculate scores and p-values for all positions ---
	results = []
	for pos, p in sorted(profile.items()):
		# biopython is off by one position
		pos_corrected = pos + 1

		# Determine major allele and its frequency
		major_allele, major_count, count = None, 0, 0
		for a, err_list in p['err'].items():
			if a != 'I': # Exclude insertions from frequency calculation
				num_bases = len(err_list)
				if major_allele is None or num_bases > major_count:
					major_allele, major_count = a, num_bases
				count += num_bases
		
		# Avoid division by zero for positions with no coverage
		if count == 0:
			continue

		major_f = major_count / count
		minor_f = 1 - major_f

		# --- Calculate required coverage based on observed minor allele frequency ---
		p_min_tiers = [0.01, 0.02, 0.03, 0.04, 0.05]
		# 当 minor_f 为 0 时 (纯合位点)，计算无意义
		if minor_f > 0:
			# 如果观测频率高于最高档位，则使用最高档位
			if minor_f > p_min_tiers[-1]:
				selected_p_min = p_min_tiers[-1]
			else:
				# 否则，找到第一个大于等于观测频率的档位作为计算标准
				for tier in p_min_tiers:
					if minor_f <= tier:
						selected_p_min = tier
						break
			required_coverages = calculate_required_coverage(p_min=selected_p_min, k=1)
		else:
			selected_p_min = 0
			required_coverages = {'N_req_95': 0, 'N_req_99': 0, 'N_req_999': 0}


		# Calculate LLR score (log10 likelihood ratio)
		# H1: Heteroplasmy model (frequency = major_f)
		# H0: Homoplasmy model (frequency = 1.0)
		llr_score = log_likelihood(major_allele, major_f, p['err']) - log_likelihood(major_allele, 1, p['err'])

		# Convert LLR score to p-value using the chi-squared distribution
		# Test statistic T = 2 * ln(LR) = 2 * ln(10) * log10(LR) ~= 4.605 * LLR
		# T follows a chi-squared distribution with df=1
		if llr_score > 0:
			test_statistic = 4.60517 * llr_score
			p_value = chi2.sf(test_statistic, df=1) # sf is the survival function (1 - cdf)
		else:
			p_value = 1.0

		# Store results for this position
		position_data = {
			'Pos': pos_corrected,
			'Ref': p['ref'],
			'A': len(p['err'].get('A', [])),
			'C': len(p['err'].get('C', [])),
			'G': len(p['err'].get('G', [])),
			'T': len(p['err'].get('T', [])),
			'D': len(p['err'].get('D', [])),
			'I': len(p['err'].get('I', [])),
			'Total': count,
			'Percentage': minor_f,
			'Ea': statistics.mean(p['err']['A']) if p['err'].get('A') else 0,
			'Ec': statistics.mean(p['err']['C']) if p['err'].get('C') else 0,
			'Eg': statistics.mean(p['err']['G']) if p['err'].get('G') else 0,
			'Et': statistics.mean(p['err']['T']) if p['err'].get('T') else 0,
			'Ed': statistics.mean(p['err']['D']) if p['err'].get('D') else 0,
			'Ei': statistics.mean(p['err']['I']) if p['err'].get('I') else 0,
			'Score': llr_score,
			'Score_p_value': p_value,
			'Selected_Pmin_Tier': selected_p_min,
			'N_req_95': required_coverages['N_req_95'],
			'N_req_99': required_coverages['N_req_99'],
			'N_req_999': required_coverages['N_req_999'],
		}
		results.append(position_data)

	# --- FDR CALCULATION ---
	if not results:
		print("No processable results, program exiting.")
		return

	p_values = [res['Score_p_value'] for res in results]
	
	# Use fdrcorrection from statsmodels. It returns a tuple of (rejected_flags, corrected_pvals)
	# The corrected_pvals are the q-values.
	rejected, q_values = fdrcorrection(p_values, alpha=0.05, method='indep')

	# Add q_values to the results with the new prefixed key
	for i, res in enumerate(results):
		res['Score_q_value'] = q_values[i]
	
	# --- PASS 2: Write all results to the output file ---
	with open(out_csv, 'w') as f:
		# Define header, now including the new columns for required coverage
		header = [
			"Pos", "Ref", "A", "C", "G", "T", "D", "I", "Total", 
			"Percentage", "Ea", "Ec", "Eg", "Et", "Ed", "Ei", 
			"Score", "Score_p_value", "Score_q_value", "Selected_Pmin_Tier",
			"N_req_95", "N_req_99", "N_req_999", "GeneProduct"
		]
		f.write(','.join(header) + '\n')
		
		# Get all annotations at once to avoid re-reading the file
		annotations = annotate.get_entries(ANNOTATION_FILE)

		for res in results:
			# Find annotation for the position
			loc = annotate.search(annotations, res['Pos'])
			gene_product = "non-coding" if loc == -1 else annotations[loc][2]
			
			# Write data row in the same order as the header
			row_data = [
				res['Pos'], res['Ref'], res['A'], res['C'], res['G'], res['T'], 
				res['D'], res['I'], res['Total'], res['Percentage'], res['Ea'], 
				res['Ec'], res['Eg'], res['Et'], res['Ed'], res['Ei'],
				res['Score'], res['Score_p_value'], res['Score_q_value'],
				res['Selected_Pmin_Tier'], res['N_req_95'], res['N_req_99'],
				res['N_req_999'], gene_product
			]
			f.write(','.join(map(str, row_data)) + '\n')
	
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
