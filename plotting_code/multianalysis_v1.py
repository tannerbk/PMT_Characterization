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
		print "opening file..."
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
	print "the bin width is: ", width
	print "the length of b is ", len(b)
	print "the length of binheight is ", len(bin_height)

	return bin_boundary, bin_height, width

#def stack_peaks(x): this will be a future implementation
#	index_min = max(range(len(values)), key=values.__getitem__)

def plot(x):

	t_bin_boundary = []
	t_bin_height = []
	t_bin_width = []

	ch_bin_boundary = []
	ch_bin_height = []
	ch_bin_width = []

	for e in range(len(x)):
		t_out = normalize(x[e][0], 10, 20, .1)
		ch_out = normalize(x[e][1], 0, 6, .1)		

		t_bin_boundary.append(t_out[0])
		t_bin_height.append(t_out[1])
		t_bin_width.append(t_out[2])

		ch_bin_boundary.append(ch_out[0])
		ch_bin_height.append(ch_out[1])
		ch_bin_width.append(ch_out[2])

	t_bin_boundary, t_bin_height, t_bin_width, \
	ch_bin_boundary, ch_bin_height, ch_bin_width

	fig1, ax1 = plt.subplots()
	fig2, ax2 = plt.subplots()

	for set in range(len(t_bin_boundary)):

		fig1.suptitle("Delta t")
		#ax1.bar(t_bin_boundary[set][:-1], t_bin_height[set],\
		#width= t_bin_width[set])
		ax1.step(t_bin_boundary[set][:-1],t_bin_height[set],\
		where = 'mid')
		ax1.set_xlabel("time [ns]")
		ax2.set_ylabel("Relative Intensity")

		fig2.suptitle("Charge Peak")
		#ax2.bar(ch_bin_boundary[set][:-1], ch_bin_height[set],\
		#width = t_bin_width[set])
		ax2.step(ch_bin_boundary[set][:-1],ch_bin_height[set],\
		where = 'mid')
		ax2.set_xlabel("Charge [pC]")
		ax2.set_ylabel("Relative Intensity")

	plt.show()

if __name__=='__main__':
	
	print "Directory: ", sys.argv[1]
	data = read_dir(sys.argv[1])
	plot(data)
	
