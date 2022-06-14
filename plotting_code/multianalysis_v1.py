# this is a test python code developed by mduce to learn how to use pyroot
# the import file to the code simply must be a .root
# for now, i've dumped ZT0076_1650V_Neck_and_Cone_lappd_0_gr0_ch1.root in .
# currently, this is for TWO DATASETS

import ROOT # import the CERN ROOT module
import sys
from matplotlib import pyplot as plt
import numpy as np

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
		
	return deltat_c, charge_c

def open_root(fi1):
	# open ROOT file
	f1 = ROOT.TFile.Open(fi1)
	# open ROOT tree
	tree1 = f1.Get("output")
	# Get the data from the ROOT tree using open_tree function
	ht1 = open_tree(tree1)
	return ht1

def normalize(data):
	# crop dataset so that all datasets can be plotted together
	# limit deltat values to those between 0 and X ns
	crop_min = 10
	crop_max = 20
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
	b = np.arange(float(min(y1)),float(max(y1)) + .2, .2)

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

def plot(y1, y2, y3, y4):
	f1 = str(input("Input legend value for dataset 1:\n"))
	f2 = str(input("Input legend value for dataset 2:\n"))
	f3 = str(input("Input legend value for dataset 3:\n"))
	f4 = str(input("Input legend value for dataset 4:\n"))	

	bb1, bh1, bw1 = normalize(y1)
	bb2, bh2, bw2 = normalize(y2)
	bb3, bh3, bw3 = normalize(y3)
	bb4, bh4, bw4 = normalize(y4)
	
	plt.figure(1)
	plt.step(bb1[:-1],bh1,'k',linestyle='-',linewidth=2,where='mid')
	plt.bar(bb1[:-1],bh1,width=bw1,linewidth=0,alpha=0.3)
	plt.step(bb2[:-1],bh2,'k',linestyle='--',linewidth=2,where='mid')
	plt.bar(bb2[:-1],bh2,width=bw2,linewidth=0,alpha=0.3)
	plt.step(bb3[:-1],bh3,'k',linestyle=':',linewidth=2,where='mid')
	plt.bar(bb3[:-1],bh3,width=bw3,linewidth=0,alpha=0.3)
	plt.step(bb4[:-1],bh4,'k',linestyle='-.',linewidth=2,where='mid')
	plt.bar(bb4[:-1],bh4,width=bw4,linewidth=0,alpha=0.3)

	#plt.bar(bb1[:-1], bh1, width = bw1)
	#plt.bar(bb2[:-1], bh2, width = bw2)
	#plt.bar(bb3[:-1], bh3, width = bw4)
	#plt.bar(bb4[:-1], bh4, width = bw4)
	
	plt.xlabel("Time [ns]")
	plt.ylabel("Relative Intensity")
	plt.legend([f1, f2, f3, f4])
	plt.show()

def plotnormed(y1):
	# this function plots a histogram using the built-in
	# 'normed' functionality which normalizes such that
	# the binwidth*binheight = 1 for each hist. this 
	# differs from the normalization conducted in
	# normalize where sum(all bin heights) = 1
	plt.figure(2)
	plt.hist(y1,bins=1000,normed=True)
	plt.show()

if __name__=='__main__':
	# identify number of datasets to be compared, print
	numarg = len(sys.argv)-1
	#print "Compared Data", sys.argv
	print "number of Datasets", numarg
	
	#for e through the number of datasets to be compared
	# open_root file, tree, data using open_tree
	#for i in range(numarg):
	dt1, ch1 = (open_root(sys.argv[1]))
	dt2, ch2 = (open_root(sys.argv[2]))
	dt3, ch3 = (open_root(sys.argv[3]))
	dt4, ch4 = (open_root(sys.argv[4]))
	plot(dt1, dt2, dt3, dt4)
	#plotnormed(dt1)
