from bokeh.models import HoverTool, NumeralTickFormatter, FuncTickFormatter, FixedTicker
from bokeh.models import ColumnDataSource, LabelSet, TapTool, Spacer
from bokeh.models import LinearColorMapper, Range1d, Circle
from bokeh.models.widgets import AutocompleteInput, Button, TextInput, Slider, CheckboxGroup
from bokeh.models.callbacks import CustomJS
from bokeh.plotting import figure, show, output_file
from bokeh.palettes import Blues9, OrRd9, YlOrRd9, Accent8, BuGn9, Set1
from bokeh.layouts import row, column, widgetbox
from bokeh import events

import pandas
import numpy
import sys
import argparse
import csv
import locale
import numpy as np
import math

#------------------------------------------------------------------------------
# Used to set y range of heteroplasmy plot
# Plots (1,0) and (1,1)
#------------------------------------------------------------------------------
VISIBLE_SAMPLE_RANGE = (0,36)
MAX_X = 1

#------------------------------------------------------------------------------
# The entire figure is a 2x3 grid
# DIM[row,column] = (width,height) of plot at location (row,column)
#------------------------------------------------------------------------------
DIM = {
	(0,0) : (1050, 70),
	(0,1) : ( 120, 70),
	(1,0) : (1050,550),
	(1,1) : ( 120,550),
	(2,0) : (1050, 90),
	(2,1) : ( 120, 90),
}

#------------------------------------------------------------------------------
# GENE_INTERVAL[gene_symbol] = [min, max, min_0, max_0, min_1, max_1, ... ]
# Used to label gene products returned by gene zooming search
# Plot (2,0)
#------------------------------------------------------------------------------
GENE_INTERVAL = {}

#------------------------------------------------------------------------------
# HETEROPLASMY_PROBABILITIES[coord] = dict of [(sample_id, probs), .... ]
# Used to plot hbar of prob distributions of heteroplasmy of samples
# Plot (1,1)
#------------------------------------------------------------------------------
HETEROPLASMY_PROBABILITIES = {}

#------------------------------------------------------------------------------
# Command line arguments
#------------------------------------------------------------------------------
ARGS = None

#------------------------------------------------------------------------------
def get_cmd_args():
	parser = argparse.ArgumentParser(description='Create heteroplasmy plot')
	parser.add_argument('genome_name', help='Name of genome (a string)')
	parser.add_argument('genome_annotations', help='Annotations of gene products (csv)')
	parser.add_argument('heteroplasmies', help='Heteroplasmies file (csv)')
	parser.add_argument('conserved_scores', help='conserved scores file (csv)')
	parser.add_argument('output', help='Output file')
	return parser.parse_args()

#------------------------------------------------------------------------------
# MAIN
#------------------------------------------------------------------------------
def main():
	global ARGS
	ARGS = get_cmd_args()

	# Plot heteroplasmic sites
	plasmy_fig, plasmy_source = plot_heteroplasmies()

	# Exit if the input file is empty, which would cause a crash.
	if plasmy_source is None:
		print(f"Error: Input file '{ARGS.heteroplasmies}' is empty or contains no data. Cannot generate plot.", file=sys.stderr)
		# Create an empty html file to satisfy workflow dependencies if needed
		output_file(ARGS.output, title='No Data')
		# An empty show() call might be needed to write the file.
		show(column(figure(title="No Data Available in Input File")))
		return

	annotation_fig = plot_genome_annotations(plasmy_fig)
	label_source, line_source = plot_search_result_annotations(annotation_fig)

	# Plot heteroplasmy probability figure and conservation annotation (in main figure)
	prob_fig, prob_data_source = build_prob_figure(plasmy_fig)
	conservation_fig = plot_conservation_annotations(plasmy_fig, prob_data_source)

	# coverage filter
	coverage_filter1 = build_coverage_filter(plasmy_source)

	# Build widgets
	DI_box = build_DI_optional_box(plasmy_source)
	coor_input = build_search_coordinate(plasmy_fig, line_source)
	search_input = build_search_input(plasmy_fig, label_source, line_source)
	clear_button = build_clear_button(label_source, line_source)

	# Layout figures and widgets
	layout_plots(
		plasmy_fig, 
		conservation_fig, 
		annotation_fig, 
		prob_fig, 
		coverage_filter1, 
		coor_input,
		search_input, 
		clear_button,
		DI_box
	)

#------------------------------------------------------------------------------
def acgt_color(base):
	color = dict(A=Set1[7][0], C=Set1[7][1], G=Set1[7][2], T=Set1[7][3], D=Set1[7][4])
	return color[base]

def plasmy_color(row):
	if row['A']>row['C'] and row['A']>row['G'] and row['A']>row['T'] and row['A']>row['D']:
		return acgt_color('A')
	if row['C']>row['A'] and row['C']>row['G'] and row['C']>row['T'] and row['C']>row['D']:
		return acgt_color('C')
	if row['G']>row['A'] and row['G']>row['C'] and row['G']>row['T'] and row['G']>row['D']:
		return acgt_color('G')
	if row['T']>row['A'] and row['T']>row['C'] and row['T']>row['G'] and row['T']>row['D']:
		return acgt_color('T')
	return acgt_color('D')

#------------------------------------------------------------------------------
def certainty(p):
	# Using np.log2 for vectorized operation if possible, but list comprehension is fine
	return 2.0 - sum([ -q*np.log2(q) for q in p if q>0] )

def plasmy_alpha(row):
	certainty_int = [certainty([0,0,0.5,0.5]), certainty([0,0,0.05,0.95])]   
	alpha_int = [0.4,1]
	min_alpha = 0.1
	h = certainty([row['A'],row['C'],row['G'],row['T'],row['D']])
	return np.interp(h, certainty_int, alpha_int, left=min_alpha, right=1)

#------------------------------------------------------------------------------
# LAYOUT FIGURES AND WIDGETS
#------------------------------------------------------------------------------
def layout_plots(plasmy_fig, conservation_fig, annotation_fig, prob_fig, coverage_filter1, coor_input, search_input, clear_button,DI_box):
	acgt = figure(
		plot_width = DIM[0,1][0],
		plot_height = DIM[0,1][1],
		x_range = (0,6),
		y_range = (0.3,3),
		toolbar_location=None,
	)
	acgt.xgrid.grid_line_color = None
	acgt.ygrid.grid_line_color = None
	acgt.xaxis.visible = False
	acgt.xgrid.grid_line_color = None
	acgt.yaxis.visible = False
	acgt.ygrid.grid_line_color = None
	acgt.outline_line_width = 1
	acgt.outline_line_alpha = 0.5
	acgt.outline_line_color = 'gray'

	source_A = ColumnDataSource(data=dict(
		x=[1,2,3,4,5],
		y=[1,1,1,1,1],
		text=['A','C','G','T','D'],
		text_color=[acgt_color('A'), acgt_color('C'), acgt_color('G'), acgt_color('T'), acgt_color('D')],
	))
	lab_A = LabelSet(
		x='x',y='y',text='text',text_color='text_color',text_align='center',
		text_font_style = 'bold',
		source=source_A, level='glyph', render_mode='canvas')
	acgt.add_layout(lab_A)
	
	# --- Add a dummy, invisible renderer to prevent W-1000 warning ---
	acgt.circle(x=[], y=[], alpha=0)

	layout = column(
		row(
			column(plasmy_fig, conservation_fig, annotation_fig),
			column(prob_fig, acgt, widgetbox(clear_button, width=70)),
			column(widgetbox(coverage_filter1,DI_box,coor_input,search_input, width=200)),
		),
	)
	print('Saved to', ARGS.output)
	output_file(ARGS.output, mode='inline', title='Heteroplasmy in %s' % ARGS.genome_name)
	show(layout)

#------------------------------------------------------------------------------
# PLOT HETEROPLASMIC SITES
#------------------------------------------------------------------------------
def plot_heteroplasmies():
	global MAX_X, VISIBLE_SAMPLE_RANGE

	try:
		plasmy_df = pandas.read_csv(ARGS.heteroplasmies)
		if plasmy_df.empty:
			return None, None # Return None if dataframe is empty
	except (FileNotFoundError, pandas.errors.EmptyDataError):
		return None, None # Return None if file not found or is empty

	plasmy_df['color'] = [ plasmy_color(r[1]) for r in plasmy_df.iterrows() ]
	plasmy_df['alpha'] = [ plasmy_alpha(r[1]) for r in plasmy_df.iterrows() ]
	plasmy_df['alpha_original'] = plasmy_df['alpha'].copy()
		
	plasmy_source = ColumnDataSource(data=plasmy_df.to_dict('list'))
	if not plasmy_df.empty:
		if plasmy_df.max()['Coordinate'] > MAX_X:
			MAX_X = plasmy_df.max()['Coordinate']

		if VISIBLE_SAMPLE_RANGE[1] > plasmy_df['Sample'].max() + 1:
			VISIBLE_SAMPLE_RANGE = (VISIBLE_SAMPLE_RANGE[0] , plasmy_df['Sample'].max() + 1)

	p_hover = HoverTool(
		tooltips = [
			('Sample', '@Name'),
			('Coordinate', '@Coordinate'),
			('Gene Product', '@GP'),
			('A', '@A{1.1111}'),
			('C', '@C{1.1111}'),
			('G', '@G{1.1111}'),
			('T', '@T{1.1111}'),
			('D', '@D{1.1111}'),
			('I', '@I{1.1111}'),
			('Coverage', '@Total'),
			('NN distance', '@d'),
		],
		names = [ 'plasmy' ],
	)

	# --- Simplified tool string to prevent UserWarning ---
	fig = figure(
		title='Heteroplasmy in %s' % ARGS.genome_name,
		plot_width = DIM[1,0][0],
		plot_height = DIM[1,0][1],
		tools=["pan,wheel_zoom,box_zoom,undo,reset,save",p_hover],
		active_scroll="wheel_zoom",
		y_range = VISIBLE_SAMPLE_RANGE,
		output_backend="webgl",
		toolbar_location="above",
	)
	fig.xgrid.grid_line_color = None
	fig.xaxis.visible = False
	fig.ygrid.grid_line_color = None

	person_id = plasmy_df['Sample']
	person_name = plasmy_df['Name']
	y_ticks_labels = { person_id[i] : person_name[i] for i in range(len(person_id)) }
	fig.axis.ticker = FixedTicker(ticks=person_id)
	fig.yaxis.formatter = FuncTickFormatter(code="""
    var labels = %s;
    if (tick in labels) {
        return labels[tick];
    }
    return '';
""" % y_ticks_labels )

	fig.outline_line_width = 1
	fig.outline_line_alpha = 0.5
	fig.outline_line_color = 'gray'
	fig.circle(
		x = 'Coordinate',
		y = 'Sample',
		color = 'color',
		alpha = 'alpha',
		size = 6,
		name = 'plasmy',
		source = plasmy_source,
	)

	global HETEROPLASMY_PROBABILITIES
	g = plasmy_df[['Coordinate','Sample','A','C','G','T','D','I']].groupby('Coordinate')
	for gid in g.groups:
		rows = g.get_group(gid).iterrows()
		HETEROPLASMY_PROBABILITIES[gid] = [
			[r[1]['Sample'],r[1]['A'],r[1]['C'],r[1]['G'],r[1]['T'],r[1]['D'],r[1]['I']] for r in rows
		]

	return fig, plasmy_source

#------------------------------------------------------------------------------
# PLOT GENOME ANNOTATIONS
#------------------------------------------------------------------------------
def plot_genome_annotations(main_fig):
	global GENE_INTERVAL, MAX_X
	color_scheme = Accent8
	GENEPRODUCT_COLOR_SCHEME = dict(
		gene = color_scheme[0],
		exon = color_scheme[0],
		CDS = color_scheme[3],
		rRNA = color_scheme[1],
		tRNA = color_scheme[2],
		repeat_region = color_scheme[7],
	)
	entries = []
	locale.setlocale( locale.LC_ALL, 'en_US.UTF-8' )
	fill_alpha = dict(rRNA=1,tRNA=1,exon=0.25,gene=0.9,CDS=0.9,repeat_region=0.3)
	with open(ARGS.genome_annotations) as file:
		reader = csv.DictReader(file)
		for row in reader:
			if row['Direction'] in ['forward','reverse'] and row['Type'] in GENEPRODUCT_COLOR_SCHEME:
				if locale.atoi(row['Maximum']) > MAX_X:
					MAX_X = locale.atoi(row['Maximum'])
				y_coord = 1 if row['Direction']=='forward' else 0
				height =	0.2 if row['Type']=='repeat_region' else 1
				entries.append((
					locale.atoi(row['Minimum']),
					locale.atoi(row['Maximum']),
					row['Name'].strip(),
					row['Type'].strip(),
					row['# Intervals'],
					row['Direction'],
					y_coord,
					GENEPRODUCT_COLOR_SCHEME[row['Type']],
					height,
					fill_alpha[row['Type']],
					2 if row['Type']=='exon' else 1, # linewidth
				))
				if row['Type'] != 'exon':
					int_min, int_max = locale.atoi(row['Minimum']), locale.atoi(row['Maximum'])
					name = row['Name'].split(' ')[0].strip()
					if name not in GENE_INTERVAL:
						GENE_INTERVAL[name.lower()] = [name, int_min, int_max, int_min, int_max]
					else:
						GENE_INTERVAL[name.lower()].extend([int_min, int_max])
						GENE_INTERVAL[name.lower()][1] = min(int_min, GENE_INTERVAL[name.lower()][1])
						GENE_INTERVAL[name.lower()][2] = max(int_max, GENE_INTERVAL[name.lower()][2])
	entries = sorted(entries)

	source = ColumnDataSource(data=dict(
		x = [ (a[0]+a[1])//2 for a in entries ],
		y = [ a[6] for a in entries ],
		color = [ a[7] for a in entries ],
		height = [ a[8] for a in entries ],
		width = [ a[1]-a[0]+1 for a in entries ],
		lab = [ a[2] for a in entries ],
		min = [ a[0] for a in entries ],
		max = [ a[1] for a in entries ],
		dir = [ a[5] for a in entries ],
		fill_alpha = [ a[9] for a in entries ],
		line_alpha = [ a[9] for a in entries ],
		line_width = [ a[10] for a in entries ],
	))
	a_hover = HoverTool(
		tooltips = [
			('Name', '@lab'),
			('Location', '(@min, @max)'),
			('Direction', '@dir')
		],
		names = [ 'gene_product' ],
	)
	if MAX_X > 1000000:
		MAX_X = 1000000
	main_fig.x_range = Range1d(-0.05*MAX_X, MAX_X*1.05)

	fig = figure(
		plot_width = DIM[2,0][0],
		plot_height= DIM[2,0][1],
		x_range = main_fig.x_range,
		y_range = (-3.5,2),
		tools=['reset,tap,xwheel_zoom',a_hover],
		toolbar_location=None,
		active_scroll="xwheel_zoom",
	)
	fig.xgrid.grid_line_color = None
	fig.xaxis.axis_label_text_font_style = "normal"
	fig.xaxis.axis_label_text_font_size = "14pt"
	fig.xaxis[0].formatter = NumeralTickFormatter(format="0")
	fig.ygrid.grid_line_color = None
	fig.yaxis.visible = False
	fig.outline_line_width = 1
	fig.outline_line_alpha = 0.5
	fig.outline_line_color = 'gray'

	fig.rect(
		x = 'x',
		y = 'y',
		width = 'width',
		height = 'height',
		color = 'color',
		alpha = 'alpha',
		fill_alpha = 'fill_alpha',
		line_alpha = 'line_alpha',
		line_width = 'line_width',
		nonselection_color = 'color',
		width_units = 'data',
		height_units = 'data',
		name = 'gene_product',
		source = source
	)
	return fig

#------------------------------------------------------------------------------
# PLOT LABELLED RESULTS OF SEARCH
#------------------------------------------------------------------------------
def plot_search_result_annotations(annotation_fig):
	label_source = ColumnDataSource(data=dict(x=[],y=[],text=[]))
	line_source = ColumnDataSource(data=dict(xs=[],ys=[]))

	annotation_fig.multi_line(
		xs = 'xs', ys = 'ys', line_color = 'navy', line_dash = 'dotted',
		line_alpha = 0.5, line_width = 2, source = line_source,
	)
	labels = LabelSet(
		x = 'x', y = 'y', text = 'text', text_align = 'center',
		level = 'glyph', render_mode = 'canvas', source = label_source,
	)
	annotation_fig.add_layout(labels)
	return label_source, line_source

#------------------------------------------------------------------------------
# PLOT CONSERVATION ANNOTATIONS
#------------------------------------------------------------------------------
def plot_conservation_annotations(main_fig, targeted_source):
	df = pandas.read_csv(ARGS.conserved_scores)
	source = ColumnDataSource(data=dict(
		y = [0] * len(df),
		Coordinate = df['Coordinate'],
		Score = df['Score'],
	))
	
	# This JS code will be used by both tap and box_select callbacks
	callback_code = """
		const inds = source.selected.indices;
		const selected_data = source.data;
		const targeted_data = targeted_source.data;
		
		targeted_data['y'] = [];
		targeted_data['left'] = [];
		targeted_data['right'] = [];
		targeted_data['height'] = [];
		targeted_data['color'] = [];

		const samples = {};
		for (const j of inds) {
			const coord = selected_data['Coordinate'][j];
			const items = prob[coord];
			if (!items) continue;
			for (const v of items) {
				const sample_id = v[0];
				if (sample_id in samples) {
					samples[sample_id][0] += 1;
					for (let k=1; k<v.length; k++) {
						samples[sample_id][k] += v[k];
					}
				} else {
					samples[sample_id] = [1, ...v.slice(1)];
				}
			}
		}

		for (const s in samples) {
			if (Object.prototype.hasOwnProperty.call(samples, s)) {
				const u = samples[s];
				const count = u[0];
				const v = [u[1]/count, u[2]/count, u[3]/count, u[4]/count, u[5]/count];
				const y = parseInt(s);
				targeted_data['y'].push(y,y,y,y,y);
				targeted_data['left'].push(0,v[0],v[0]+v[1],v[0]+v[1]+v[2],v[0]+v[1]+v[2]+v[3]);
				targeted_data['right'].push(v[0],v[0]+v[1],v[0]+v[1]+v[2],v[0]+v[1]+v[2]+v[3],1);
				targeted_data['height'].push(0.9,0.9,0.9,0.9,0.9);
				targeted_data['color'].push('%s','%s','%s','%s','%s');
			}
		}
		targeted_source.change.emit();
	""" % (acgt_color('A'),acgt_color('C'),acgt_color('G'),acgt_color('T'),acgt_color('D'))

	callback_args = dict(
		source=source,
		targeted_source=targeted_source,
		prob=HETEROPLASMY_PROBABILITIES
	)

	# This callback responds to Tapping
	source.selected.js_on_change('indices', CustomJS(args=callback_args, code=callback_code))

	c_hover = HoverTool(
		tooltips = [
			('Coordinate', '@Coordinate'),
			('Conservation', '@Score{0,0.0000}'),
		],
		names = [ 'conserved' ],
	)

	fig = figure(
		plot_width = DIM[0,0][0],
		plot_height= DIM[0,0][1],
		x_range = main_fig.x_range,
		tools=['tap,box_select,xwheel_zoom',c_hover],
		toolbar_location=None,
		active_scroll='xwheel_zoom',
		active_tap='tap',
		active_drag='box_select',
		output_backend="webgl",
	)

	# This callback responds to Box Selecting
	fig.js_on_event(events.SelectionGeometry, CustomJS(args=callback_args, code=callback_code))

	fig.xgrid.grid_line_color = None
	fig.ygrid.grid_line_color = None
	fig.outline_line_width = 1
	fig.outline_line_alpha = 0.5
	fig.outline_line_color = 'gray'
	fig.xaxis.visible = False
	fig.yaxis.visible = False

	PALETTE_CONSERVATION_SCORE = YlOrRd9[::-1][2:]
	c_mapper = LinearColorMapper(PALETTE_CONSERVATION_SCORE, low=0, high=1)
	fig.square(
		x = 'Coordinate',
		y = 'y',
		color = {'field':'Score','transform':c_mapper},
		alpha = 1,
		size = 6,
		name = 'conserved',
		source = source,
	)
	return fig

#------------------------------------------------------------------------------
# This figure provides annotation of heteroplasmy probabilites across samples.
#------------------------------------------------------------------------------
def build_prob_figure(main_fig):
	fig = figure(
		plot_width =  DIM[1,1][0],
		plot_height = DIM[1,1][1],
		x_range = (0,1),
		y_range = main_fig.y_range,
		tools=[],
		toolbar_location=None,
	)
	fig.outline_line_width = 1
	fig.outline_line_alpha = 0.5
	fig.outline_line_color = 'gray'
	fig.xaxis.visible = False
	fig.xgrid.grid_line_color = None
	fig.yaxis.visible = False
	fig.ygrid.grid_line_color = None

	prob_source = ColumnDataSource(data=dict(y=[],left=[],right=[],height=[],color=[]))
	
	fig.hbar(y='y',left='left',right='right',color='color',height='height',source=prob_source)
	return fig, prob_source

#------------------------------------------------------------------------------
# Search provides zooming into a gene
#------------------------------------------------------------------------------
def build_search_input(main_fig, label_source, line_source):
	text = AutocompleteInput(
		title = 'Locate gene',
		value = '',
		placeholder = 'Gene symbol',
		completions = [ v[0] for k,v in GENE_INTERVAL.items() ],
	)
	text.js_on_change('value', CustomJS(
		args = dict(
			x_range = main_fig.x_range,
			label_source = label_source,
			line_source = line_source,
			interval = GENE_INTERVAL,
			y = -3
		),
		code="""
		const gene_symbol = this.value.toLowerCase();
		const data = label_source.data;
		const data_line = line_source.data;
		
		// --- Clear previous labels before drawing new ones ---
		data['x'] = [];
		data['y'] = [];
		data['text'] = [];
		data_line['xs'] = [];
		data_line['ys'] = [];
		
		if (gene_symbol.length > 0 && gene_symbol in interval) {
			const name = interval[gene_symbol][0];
			x_range.start = interval[gene_symbol][1] * 0.99;
			x_range.end = interval[gene_symbol][2] * 1.01;
			for (let i=3; i<interval[gene_symbol].length; i += 2){
				const start = interval[gene_symbol][i];
				const end = interval[gene_symbol][i+1];
				data['x'].push((start+end)*0.5);
				data['y'].push(y);
				data['text'].push(name);
				data_line['xs'].push([start, end]);
				data_line['ys'].push([y,y]);
			}
		}
		// Always emit a change to clear the plot if search is empty
		label_source.change.emit();
		line_source.change.emit();
		"""))
	return text

#------------------------------------------------------------------------------
# Search provides zooming into a specific coordinate
#------------------------------------------------------------------------------
def build_search_coordinate(main_fig, line_source):
	coor_input = TextInput(value = '', title='Locate coordinate', placeholder = 'Coordinate')
	coor_input.js_on_change('value', CustomJS(
		args = dict(
			x_range = main_fig.x_range,
			line_source = line_source,
		),
		code="""
		const coor = parseInt(this.value, 10);	
		if (!isNaN(coor) && coor > 0) {	
			x_range.start = coor - 1000;
			x_range.end = coor + 1000;
			// No need to emit change on x_range, it happens automatically.
		}
		"""))

	return coor_input

#------------------------------------------------------------------------------
# THIS CLEARS GENE NAMES LABELED BY SEARCH
#------------------------------------------------------------------------------
def build_clear_button(label_source, line_source):
	button = Button(label='Clear gene labels', width=70)
	button.js_on_event(events.ButtonClick, CustomJS(
		args = dict(
			label_source = label_source,
			line_source = line_source,
		),
		code="""
        // --- Use the most robust method to clear the data sources ---
		label_source.data['x'] = [];
		label_source.data['y'] = [];
		label_source.data['text'] = [];
		label_source.change.emit();

		line_source.data['xs'] = [];
		line_source.data['ys'] = [];
		line_source.change.emit();
	"""))
	return button

#------------------------------------------------------------------------------
# COVERAGE FILTER
#------------------------------------------------------------------------------
def build_coverage_filter(plasmy_source):
	def roundup(x):
		return int(math.ceil(x / 100.0)) * 100

	max_coverage = 0
	if plasmy_source.data['Total']: # Check if list is not empty
		max_coverage = max(plasmy_source.data['Total'])
	
	slider_callback = CustomJS(args=dict(source=plasmy_source), code="""
		const data = source.data;
		const slider_value = this.value;
		const total = data['Total'];
		const alpha = data['alpha'];
		const alpha_original = data['alpha_original'];
		for (let i = 0; i < total.length; i++) {
			if (total[i] < slider_value) {
				alpha[i] = 0;
			} else {
				alpha[i] = alpha_original[i];
			}
		}
		source.change.emit();
	""")

	slider1 = Slider(start=0, end=roundup(max_coverage), value=0, step=100, title="Coverage", width=200)
	slider1.js_on_change('value', slider_callback)
	
	return slider1

#------------------------------------------------------------------------------
# Deletion and Insertion box
#------------------------------------------------------------------------------
def build_DI_optional_box(plasmy_source):
	checkbox_callback = CustomJS(args=dict(source=plasmy_source), code="""
		const data = source.data;
		const active = this.active; // active is a list of indices, e.g. [0, 2]

		const alpha = data['alpha'];
		const alpha_original = data['alpha_original'];
		const delet = data['CountD'];
		const ins = data['CountI'];
		
		const show_del = active.includes(0);
		const show_ins = active.includes(1);
		const show_sub = active.includes(2);

		for (let i = 0; i < alpha.length; i++) {
			const is_del = delet[i] > 0;
			const is_ins = ins[i] > 0;
			// Substitution is defined as not being a deletion or insertion
			const is_sub = !is_del && !is_ins;

			let is_visible = false;
			if (is_sub && show_sub) {
				is_visible = true;
			}
			if (is_del && show_del) {
				is_visible = true;
			}
			if (is_ins && show_ins) {
				is_visible = true;
			}
			
			alpha[i] = is_visible ? alpha_original[i] : 0;
		}
		source.change.emit();
	""")

	checkbox = CheckboxGroup(labels=["Deletion sites", "Insertion sites", "Substitution"], active=[0, 1, 2])
	checkbox.js_on_change('active', checkbox_callback)
	
	return checkbox

#------------------------------------------------------------------------------

main()

