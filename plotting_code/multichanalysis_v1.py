# this is python code was developed by mduce so multiple pmt datasets can be
# compared after they have been run through analysis.py
# input: four .root files outputted from analysis.py
# output: layered histogram of four datasets for desired parameter
# 
# NOTE: code requires input of a directory with the to-be-compared .root files

import ROOT # import the CERN ROOT module
import sys
from matplotlib import pyplot as plt
import numpy as np

# jordan 
import glob

def open_tree(tree):
	
	# print total number of recorded events
	print "Entries:", tree.GetEntries() 
	
	# initialize lists for data storage
	charge_c=[] # charge, measurment channel
	deltat_c=[] # the time btwn measurement PMT pulse and trigger PMT pulse

	# for 1:number of tot recorded events in dataset
	for i in range(tree.GetEntries()):
		
		# Get current entry in root try
		tree.GetEntry(i)

		# check that measurement pmt voltage exceeded 5mV
		if (tree.peak_voltage > -5.0): continue	

		#print "Deltat, Charge:", tree.deltat, tree.charge
		charge_c.append(tree.charge - tree.charge_empty)
		deltat_c.append(tree.deltat)
	print "Extracting data from root file..."
	
	return deltat_c, charge_c

def open_root(fi1):
	# open ROOT file
	f1 = ROOT.TFile.Open(fi1)
	# open ROOT tree
	tree1 = f1.Get("output")
	# Get the data from the ROOT tree using open_tree function
	ht1 = open_tree(tree1)
	return ht1

# takes in folder of .root files
# outputs [(dt1, ch1), (dt2,ch2), ... , (dtn, chn)]
# list[1][1] = dt1, list[n] = dtn, chn
def read_dir(dir_name):
	value_list = []
	files = glob.glob(dir_name + "/*.root")
	print "Opening directory..."
	for file in files:
		print "opening file ", file
		dt, ch = open_root(file)
		value_list.append((dt,ch))
		print "Data extracted from file"
	return value_list

def normalize(data, crop_min, crop_max, bw):
	# crop dataset so that all datasets can be plotted together
	# limit deltat values to those between 0 and X ns
	y1 = [x for x in data if x>=crop_min and x<=crop_max]
	# slight variations in the datasets mean t doesnt exactly
	# equal zero anywhere, so we'll force the min and maxes
	# of each dataset to be 0 and X ns - this is to make
	# the bin widths consistent across all plots
	y1[0] = crop_min
	y1[-1] = crop_max
	#print float(min(y1)), float(max(y1))

	# declare desired bin width and build array
	# note +val and ,val must be equal
	b = np.arange(float(min(y1)),float(max(y1)) + bw, bw)

	# build histogram of data (y1) and desired bins
	bin_height, bin_boundary = np.histogram(y1,bins=len(b))
	
	# obtain bin width, print
	width = bin_boundary[1]-bin_boundary[0]
	
	# normalize bin heights to max bin height = 1
	bin_height = bin_height/float(max(bin_height))

	# quick check of values for eased troubleshooting
	#print "the bin width is: ", width
	#print "the length of b is ", len(b)
	#print "the length of binheight is ", len(bin_height)

	return bin_boundary, bin_height, width

#def stack_peaks(x): this will be a future implementation
#	index_min = max(range(len(values)), key=values.__getitem__)

def get_peak_index(x):
	pk_val_indxs = np.where(x == 1)
	# pk is a tuple of arrays containing the indices at which x == 1
	# pk = (array([indx]),)
	pk_indx = pk_val_indxs[0][0]
	return pk_indx

def get_indexes(x):
	# input: list of bin heights

	# output of get_peak = index at binheight = 1
	pk_indx = get_peak_index(x)
	array = np.asarray(x)
	ar1 = array[:pk_indx]
	ar2 = array[pk_indx:]
	
	# get indices of ar1,ar2 where ar1[i] ~= .5
	fwhm_indx_l = (np.abs(ar1 - .5)).argmin()
	fwhm_indx_r = (np.abs(ar2 - .5)).argmin()

	# all return vals of type numpy.int64
	return fwhm_indx_l, fwhm_indx_r, pk_indx

def plot(x):
	# x is of type list
	# x[0] is of type tuple
	# x[0][0] is of type list

	# initialize lists for hist bins locations, heights, widths
	t_bin_boundary = []
	t_bin_height = []
	t_bin_width = []

	ch_bin_boundary = []
	ch_bin_height = []
	ch_bin_width = []
	
	# user plug in crop values if non-standard
	#t_min = raw_input("Lower value of time: ")
	#t_max = raw_input("Upper value of time: ")
	#ch_min = raw_input("Lower value of charge: ")
	#ch_max = raw_input("Upper value of charge: ")

	# standard lower and upper lims of histograms
	t_min = 0
	t_max = 0
	ch_min = 20
	ch_max = 100
	
	bin_width = .1

	print "normalizing datasets..."
	# crop and normalize datasets, allocate tuple values
	for e in range(len(x)):
		# normalize and crop data
		t_out = normalize(x[e][0], t_min, t_max, bin_width)
		ch_out = normalize(x[e][1], ch_min, ch_max, bin_width)		

		t_bin_boundary.append(t_out[0])
		t_bin_height.append(t_out[1])
		t_bin_width.append(t_out[2])

		ch_bin_boundary.append(ch_out[0])
		ch_bin_height.append(ch_out[1])
		ch_bin_width.append(ch_out[2])

	# t_bin_item is a list of np arrays
	t_bin_boundary, t_bin_height, t_bin_width, \
	ch_bin_boundary, ch_bin_height, ch_bin_width

	print t_bin_width, ch_bin_width

	fwhm_index = []
	fwhm = []
	peak_indexes = []
	hold = []

	#for height_set in range(len(t_bin_height)):

	#	hold.append(get_indexes(t_bin_height[height_set]))
	#	fwhm_index.append(hold[height_set][1]-hold[height_set][0])
	#	fwhm.append(fwhm_index[height_set] * t_bin_width[height_set])
	#	peak_indexes.append(hold[height_set][2])

	#fwhm, peak_indexes
	
	#print "fwhm is ", fwhm
	
	fig1, ax1 = plt.subplots()
	fig2, ax2 = plt.subplots()
	labels = []

	for set in range(len(t_bin_boundary)):
		
		fig1.suptitle("Delta t")
		#ax1.bar(t_bin_boundary[set][:-1], t_bin_height[set],\
		#width= t_bin_width[set])
		ax1.step(t_bin_boundary[set][:-1],t_bin_height[set],\
		where = 'mid',label='%s data' % set)
		ax1.set_xlabel("time [ns]")
		ax2.set_ylabel("Relative Intensity")

		fig2.suptitle("Charge Peak")
		#ax2.bar(ch_bin_boundary[set][:-1], ch_bin_height[set],\
		#width = t_bin_width[set])
		ax2.step(ch_bin_boundary[set][:-1],ch_bin_height[set],\
		where = 'mid',label='%s data' % set)
		ax2.set_xlabel("Charge [pC]")
		ax2.set_ylabel("Relative Intensity")

		labels.append(raw_input("Input legend label for dataset: "))

	ax1.legend(labels)
	ax2.legend(labels)
	plt.show()

if __name__=='__main__':
	
	print "Directory: ", sys.argv[1]
	data = read_dir(sys.argv[1])
	plot(data)
	
