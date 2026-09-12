import csv
import os
import sys
import pandas as pd

#------------------------------------------------------------------------------
def get_csvfiles(dir):
	return [ os.path.join(dir, f) for f in os.listdir(dir) if f.endswith('.csv') and not f.endswith('_links.csv') ]

#------------------------------------------------------------------------------
def get_csvfiles_by_nameList(dir, namelist):

	def getFile(f_dir, name):
		for i in f_dir:
			if name in i:
				return i

		return None

	f = []
	f_dir = [ os.path.join(dir, f) for f in os.listdir(dir) if f.endswith('.csv') and not f.endswith('_links.csv') ]
	with open(namelist,'r') as nl:
		for line in nl:
			if "," in line:
				name, organism = line.strip().split(",")
			else:
				name = line.strip()

			if getFile(f_dir, name):
					f.append(getFile(f_dir, name))

	return f

#------------------------------------------------------------------------------
def get_individual_id(csvfile):
	filename = csvfile.split('/')[-1]
	prefix = filename.split('.')[0]
	# name = prefix.split('_')[0]+"_"+prefix.split('_')[1]
	name = prefix.split('_')[0]
	return name

#------------------------------------------------------------------------------
def load_links(link_file):
	"""
	Load linkage data with allele pairs.
	Supports:
	  - legacy: pos1,pos2,allele1,allele2,count
	  - extended (indel events): adds allele1_type, allele1_len, allele1_seq, allele2_type, allele2_len, allele2_seq
	Returns: dict mapping (pos1, pos2) -> dict of {(allele1, allele2): info_dict}
	"""
	links = {}
	try:
		with open(link_file, newline='') as f:
			reader = csv.DictReader(f)
			for row in reader:
				try:
					a = int(row['pos1'])
					b = int(row['pos2'])
				except (KeyError, ValueError):
					continue

				allele1 = (row.get('allele1') or '').strip()
				allele2 = (row.get('allele2') or '').strip()
				if not allele1 or not allele2:
					continue

				try:
					count = int(float(row.get('count', 0)))
				except (TypeError, ValueError):
					continue
				if count <= 0:
					continue

				info = {'count': count}
				# Optional metadata columns (present in the new indel-event links format)
				for k in (
					'allele1_type', 'allele1_len', 'allele1_seq',
					'allele2_type', 'allele2_len', 'allele2_seq',
				):
					if k in row and row[k] not in (None, ''):
						info[k] = row[k]

				links.setdefault((a, b), {})
				links[(a, b)][(allele1, allele2)] = info
	except FileNotFoundError:
		return None
	return links

#------------------------------------------------------------------------------
def majority_allele(profile_tuple):
	# profile_tuple indices: pos, score, pct, gp, freqA, freqC, freqG, freqT, freqD, freqI,
	# Ea, Ec, Eg, Et, Ed, Ei, countA, countC, countG, countT, countD, countI, total, q, n_req
	counts = {
		'A': profile_tuple[16],
		'C': profile_tuple[17],
		'G': profile_tuple[18],
		'T': profile_tuple[19],
	}
	maj = max(counts, key=lambda k: counts[k])
	return maj

#------------------------------------------------------------------------------
def filter_csvfile(csvfile, q_value_threshold, percentage_threshold, count_threshold_suffix):
	"""
	Filters a single CSV file based on q-value, percentage, and a dynamically chosen required coverage column.
	Returns (selected_positions, all_positions) where selected passed filters.
	"""
	selected_positions = {}
	all_positions = {}
	coverage_col = f"N_req_{count_threshold_suffix}"
	
	with open(csvfile) as f:
		reader = csv.DictReader(f)
		for row in reader:
			required_cols = ['Pos', 'Score', 'A', 'C', 'G', 'T', 'D', 'I', 'Total', 
							 'GeneProduct', 'Score_q_value', 'Percentage',
							 'N_req_95', 'N_req_99', 'N_req_999',
							 'Ea', 'Ec', 'Eg', 'Et', 'Ed', 'Ei']
			if not all(col in row and row[col] and row[col] != 'NA' for col in required_cols):
				continue # Skip rows with missing or invalid data

			pos = int(row['Pos'])
			total = float(row['Total'])
			selected_n_req_value = int(float(row[coverage_col]))
			profile_tuple = (
				pos, float(row['Score']), float(row['Percentage']), row['GeneProduct'],
				float(row['A'])/total, float(row['C'])/total,
				float(row['G'])/total, float(row['T'])/total,
				float(row['D'])/total, float(row['I'])/total,
				float(row['Ea']), float(row['Ec']), float(row['Eg']),
				float(row['Et']), float(row['Ed']), float(row['Ei']),
				int(float(row['A'])), int(float(row['C'])), int(float(row['G'])), int(float(row['T'])),
				int(float(row['D'])), int(float(row['I'])), int(total),
				float(row['Score_q_value']),
				selected_n_req_value,
				(row.get('DEvents') or '').strip(),
				(row.get('IEvents') or '').strip(),
			)
			all_positions[pos] = profile_tuple

			passes_q_value = float(row['Score_q_value']) < q_value_threshold
			passes_coverage = float(row['Total']) > float(row[coverage_col])
			passes_percentage = float(row['Percentage']) >= percentage_threshold
			
			if passes_q_value and passes_coverage and passes_percentage:
				selected_positions[pos] = profile_tuple
				
	return selected_positions, all_positions

#------------------------------------------------------------------------------
def intersect(csvfiles, score_threshold, percentage_threshold, count_threshold, d_threshold, output_file, output_dir, result_dir):
	"""
	Finds intersecting positions across multiple CSV files.
	"""
	positions = {}
	# --- Validate and sanitize count_threshold ---
	valid_thresholds = ['95', '99', '999']
	sanitized_threshold_suffix = str(count_threshold).replace('.0', '')
	if sanitized_threshold_suffix not in valid_thresholds:
		print(f"Warning: Invalid count_threshold '{count_threshold}'. Only '95', '99', or '999' are allowed. Defaulting to '95'.")
		sanitized_threshold_suffix = '95'
	
	for f in csvfiles:
		person_id = get_individual_id(f)
		out_threshold_file = os.path.join(output_dir, person_id+'_threshold.txt')
		
		# Default to the global percentage threshold
		local_percentage_threshold = percentage_threshold
		q_value_threshold = score_threshold
		
		if os.path.exists(out_threshold_file):
			# Read the values from the output file
			with open(out_threshold_file, 'r') as file:
				lines = file.readlines()
				# Still allow for a per-sample percentage threshold if provided
				local_percentage_threshold = float(lines[0].split(":")[1].strip())
		else:
			print(f"Warning: Threshold file not found for {person_id}. Using default percentage threshold.")
			
		pos_selected, pos_all = filter_csvfile(f, q_value_threshold, local_percentage_threshold, sanitized_threshold_suffix)

		# linkage source for passing sites: show allele-pair linkage between heteroplasmy sites
		link_file = f.replace('.csv','_links.csv')
		links = load_links(link_file)
		linkage_source_map = {}
		if links:
			for p in pos_selected.keys():
				linkage_parts = []
				# Check all position pairs involving p
				for (pos1, pos2), allele_counts in links.items():
					other_pos = None
					if pos1 == p:
						other_pos = pos2
						allele_idx = 1  # p is pos1, other is pos2
					elif pos2 == p:
						other_pos = pos1
						allele_idx = 0  # p is pos2, other is pos1

					# Only include if the other position is also a selected heteroplasmy site
					if other_pos is not None and other_pos in pos_selected:
						for (allele1, allele2), info in sorted(allele_counts.items(), key=lambda x: -x[1]['count']):
							count = info['count']
							if allele_idx == 1:
								linkage_parts.append(f"{p}:{allele1}-{other_pos}:{allele2}({count})")
							else:
								linkage_parts.append(f"{p}:{allele2}-{other_pos}:{allele1}({count})")
				if linkage_parts:
					linkage_source_map[p] = ";".join(linkage_parts)

		for p,profile in pos_selected.items():
			if p not in positions:
				positions[p] = ([], [])
			names_list, prof_list = positions[p]
			rescued_flag = False
			source = linkage_source_map.get(p, "")
			names_list.append(person_id)
			prof_list.append((profile, rescued_flag, source))

	scatter_plot([get_individual_id(f) for f in csvfiles], positions, output_file, d_threshold, sanitized_threshold_suffix)

#------------------------------------------------------------------------------
def apply_neighbor_distance_filter(items, d_threshold):
	"""
	Populate the auxiliary nearest-neighbor distance and apply an optional filter.

	The distance is undefined for a singleton callset. Its output value remains 0
	as a plotting/QC sentinel, but the singleton must be retained. A threshold of
	0 disables this legacy post-selection filter.
	"""
	items.sort()
	if len(items) <= 1:
		return items

	items[0][-1] = items[1][0] - items[0][0]
	for j in range(1, len(items) - 1):
		items[j][-1] = min(items[j][0] - items[j - 1][0], items[j + 1][0] - items[j][0])
	items[-1][-1] = items[-1][0] - items[-2][0]

	if d_threshold <= 0:
		return items
	return [x for x in items if x[-1] >= d_threshold]

#------------------------------------------------------------------------------
def scatter_plot(ids, positions, output_file, d_threshold, threshold_suffix):
	"""
	Generates the final output CSV, now including all N_req columns.
	"""
	points = []
	with open(output_file, 'w') as f:
		n_req_column_name = f"N_req_{threshold_suffix}"
		# Keep frequency + count columns, but drop per-allele error-rate columns (Ea/Ec/Eg/Et/Ed/Ei) in the final output.
		# Keep the auxiliary neighbor-distance column "d" for plotting/QC.
		header = (f'Coordinate,Sample,Name,GP,A,C,G,T,D,I,CountA,CountC,CountG,CountT,CountD,CountI,Total,'
				  f'Score,Percentage,Score_q_value,{n_req_column_name},DEvents,IEvents,LinkageSource,d\n')
		f.write(header)

		for i,cur_id in enumerate(ids):
			items = []
			for pos,profile in positions.items():
				if cur_id in profile[0]:
					idx = profile[0].index(cur_id)
					item, rescued, source = profile[1][idx]
					items.append([
						pos, i + 1, cur_id, item[3],  # GP
						item[4], item[5], item[6], item[7], item[8], item[9],  # frequencies A,C,G,T,D,I
						item[10], item[11], item[12], item[13], item[14], item[15],  # error rates
						item[16], item[17], item[18], item[19], item[20], item[21], item[22],  # counts + total
						item[1], item[2], item[23], item[24],  # score, percentage, q_value, n_req
						item[25], item[26], source, 0  # DEvents, IEvents, linkage source, distance placeholder
					])

			selected = apply_neighbor_distance_filter(items, d_threshold)
			
			output_format = ('%d,%d,%s,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,'
							 '%d,%d,%d,%d,%d,%d,%d,%.4f,%.4f,%.4f,%d,%s,%s,%s,%d\n')
			for x in selected:
				# x layout includes [Ea..Ei] at indices 10..15; drop them from output.
				out_x = x[:10] + x[16:]
				f.write(output_format % tuple(out_x))

	print("Finish selecting sites.\n")

def process(params):
	print(params)
	csv_dir = params['csv_dir']
	score_threshold = float(params['score_threshold'])
	percentage_threshold = float(params['percentage_threshold'])
	count_threshold = float(params['count_threshold'])
	d_threshold = float(params['d_threshold'])
	name_list = params['name_list']
	organellar_type = params['organellar_type']
	result_dir = params['result_dir']
	output_dir = params['output_dir']

	if name_list:
		files = get_csvfiles_by_nameList(csv_dir, name_list)
	else:
		files = get_csvfiles(csv_dir)

	if not os.path.exists(result_dir):
		os.makedirs(result_dir)

	output_file = os.path.join(result_dir, organellar_type+"_heteroplasmy.csv")

	intersect(files, score_threshold, percentage_threshold, count_threshold, d_threshold, output_file, output_dir, result_dir)

	return output_file
	

#------------------------------------------------------------------------------
if __name__ == '__main__':
	# organellar type is either chloroplast or mitochondria
	if len(sys.argv) != 9 and len(sys.argv) != 10:
		print("USAGE: ", sys.argv[0], "  csv_dir score_threshold percentage_threshold count_threshold d_threshold organellar_type result_dir output_dir")
		print("or")
		print("USAGE: ", sys.argv[0], "  csv_dir score_threshold percentage_threshold count_threshold d_threshold name_list.csv organellar_type result_dir output_dir")
		sys.exit(0)

	if len(sys.argv) == 9:
		name_list = None
		organellar_type = sys.argv[6]
		result_dir = sys.argv[7]
		output_dir = sys.argv[8]
	else:
		name_list = sys.argv[6]
		organellar_type = sys.argv[7]
		result_dir = sys.argv[8]
		output_dir = sys.argv[9]

	params = {
		'csv_dir': sys.argv[1],
		'score_threshold': int(sys.argv[2]),
		'percentage_threshold': float(sys.argv[3]),
		'count_threshold': float(sys.argv[4]),
		'd_threshold': float(sys.argv[5]),
		'name_list': name_list,
		'organellar_type': organellar_type,
		'result_dir': result_dir,
		'output_dir': output_dir
	}

	# if len(sys.argv) == 4:
	# 	files = get_csvfiles(sys.argv[1])
	# else:
	# 	files = get_csvfiles_by_nameList(sys.argv[1], sys.argv[4])
	
	# score_threshold = int(sys.argv[2])
	# percentage_threshold = float(sys.argv[3])
	# # filter_csvfile(files[0], 'SRR2147184')
	# intersect(files, score_threshold, percentage_threshold)
	process(params)
