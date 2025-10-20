import csv
import os
import sys
import pandas as pd

#------------------------------------------------------------------------------
def get_csvfiles(dir):
	return [ os.path.join(dir, f) for f in os.listdir(dir) if f.endswith('.csv') ]

#------------------------------------------------------------------------------
def get_csvfiles_by_nameList(dir, namelist):

	def getFile(f_dir, name):
		for i in f_dir:
			if name in i:
				return i

		return None

	f = []
	f_dir = [ os.path.join(dir, f) for f in os.listdir(dir) if f.endswith('.csv') ]
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
def filter_csvfile(csvfile, q_value_threshold, percentage_threshold, count_threshold_suffix):
	"""
	Filters a single CSV file based on q-value, percentage, and a dynamically chosen required coverage column.
	"""
	selected_positions = {}
	# Construct the column name directly from the validated suffix
	coverage_col = f"N_req_{count_threshold_suffix}"
	
	with open(csvfile) as f:
		reader = csv.DictReader(f)
		for row in reader:
			# Check if all required columns exist and are not empty/NA
			required_cols = ['Pos', 'Score', 'A', 'C', 'G', 'T', 'D', 'I', 'Total', 
							 'GeneProduct', 'Score_q_value', 'Percentage',
							 'N_req_95', 'N_req_99', 'N_req_999']
			if not all(col in row and row[col] and row[col] != 'NA' for col in required_cols):
				continue # Skip rows with missing or invalid data
				
			# --- DYNAMIC FILTERING LOGIC ---
			passes_q_value = float(row['Score_q_value']) < q_value_threshold
			passes_coverage = float(row['Total']) > float(row[coverage_col])
			passes_percentage = float(row['Percentage']) >= percentage_threshold
			
			if passes_q_value and passes_coverage and passes_percentage:
				pos = int(row['Pos'])
				total = float(row['Total'])
				selected_n_req_value = int(float(row[coverage_col]))
				
				# Reconstruct the data tuple, now including all three N_req values for output
				profile_tuple = (
					pos, float(row['Score']), float(row['Percentage']), row['GeneProduct'],
					float(row['A'])/total, float(row['C'])/total,
					float(row['G'])/total, float(row['T'])/total,
					float(row['D'])/total, float(row['I'])/total,
					int(float(row['A'])), int(float(row['C'])), int(float(row['G'])), int(float(row['T'])),
					int(float(row['D'])), int(float(row['I'])), int(total),
					float(row['Score_q_value']), 
					selected_n_req_value
				)
				selected_positions[pos] = profile_tuple
				
	return selected_positions

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
	coverage_logic_str = f"Total > N_req_{sanitized_threshold_suffix}"
	
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
			
		# Log the thresholds being used for this sample
		thresholdout_file = os.path.join(result_dir, person_id + '_thresholdout.txt')
		with open(thresholdout_file, 'w') as file:
			file.write(f"q_value_threshold: {q_value_threshold}\n")
			file.write(f"percentage_threshold: {local_percentage_threshold}\n")
			file.write(f"count_logic: {coverage_logic_str}\n")
	
		pos = filter_csvfile(f, q_value_threshold, local_percentage_threshold, sanitized_threshold_suffix)
		for p,profile in pos.items():
			if p not in positions:
				positions[p] = ([], [])
			positions[p][0].append(person_id)
			positions[p][1].append(profile)

	scatter_plot([get_individual_id(f) for f in csvfiles], positions, output_file, d_threshold, sanitized_threshold_suffix)

#------------------------------------------------------------------------------
def scatter_plot(ids, positions, output_file, d_threshold, threshold_suffix):
	"""
	Generates the final output CSV, now including all N_req columns.
	"""
	points = []
	with open(output_file, 'w') as f:
		n_req_column_name = f"N_req_{threshold_suffix}"
		header = (f'Coordinate,Sample,Name,GP,A,C,G,T,D,I,CountA,CountC,CountG,CountT,CountD,CountI,Total,'
				  f'Score,Percentage,Score_q_value,{n_req_column_name},d\n')
		f.write(header)

		for i,cur_id in enumerate(ids):
			items = []
			for pos,profile in positions.items():
				if cur_id in profile[0]:
					idx = profile[0].index(cur_id)
					item = profile[1][idx]
					# Unpack the simplified tuple (item[18] is the single N_req value)
					items.append([
						pos, i + 1, cur_id, item[3], item[4], item[5], item[6], item[7], 
						item[8], item[9], item[10], item[11], item[12], item[13], 
						item[14], item[15], item[16], item[1], item[2], item[17], 
						item[18], 0  # 0 is placeholder for distance 'd'
					])

			# Compute the distance to the nearest neighbors
			items.sort()
			if len(items) > 1:
				items[0][-1] = items[1][0] - items[0][0]
				for j in range(1,len(items)-1):
					items[j][-1] = min(items[j][0]-items[j-1][0],  items[j+1][0]-items[j][0])
				items[-1][-1] = items[-1][0] - items[-2][0]

			selected = [ x for x in items if x[-1] >= d_threshold ]
			
			output_format = ('%d,%d,%s,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%d,%d,%d,%d,%d,%d,'
							 '%.4f,%.4f,%.4f,%d,%d\n')
			for x in selected:
				f.write(output_format % tuple(x))

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
