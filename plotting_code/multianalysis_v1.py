# this is python code was developed by mduce so multiple pmt datasets can be
# compared after they have been run through analysis.py
# input: four .root files outputted from analysis.py
# output: layered histogram of four datasets for desired parameter
# 
# NOTE: currently, this is for 4 datasets and is designed to look at deltat

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

def plot(t1, t2, t3, t4, c1, c2, c3, c4):
	f1 = str(input("Input HV for dataset 1:[numerical only] \n"))
	f2 = str(input("Input HV for dataset 2 [numerical only]:\n"))
	f3 = str(input("Input HV for dataset 3 [numerical only]:\n"))
	f4 = str(input("Input HV for dataset 4 [numerical only]:\n"))	

	tbb1, tbh1, tbw1 = normalize(t1,10, 20, .1)
	tbb2, tbh2, tbw2 = normalize(t2,10, 20, .1)
	tbb3, tbh3, tbw3 = normalize(t3,10, 20, .1)
	tbb4, tbh4, tbw4 = normalize(t4,10, 20, .1)

	cbb1, cbh1, cbw1 = normalize(c1,0, 6, .1)
	cbb2, cbh2, cbw2 = normalize(c2,0, 6, .1)
	cbb3, cbh3, cbw3 = normalize(c3,0, 6, .1)
	cbb4, cbh4, cbw4 = normalize(c4,0, 6, .1)
	
	plt.figure(1)
	plt.step(tbb1[:-1],tbh1,'k',linestyle='-',linewidth=1,where='mid')
	plt.bar(tbb1[:-1],tbh1,width=tbw1,linewidth=0,alpha=0.3)
	plt.step(tbb2[:-1],tbh2,'k',linestyle='--',linewidth=1,where='mid')
	plt.bar(tbb2[:-1],tbh2,width=tbw2,linewidth=0,alpha=0.3)
	plt.step(tbb3[:-1],tbh3,'k',linestyle=':',linewidth=1,where='mid')
	plt.bar(tbb3[:-1],tbh3,width=tbw3,linewidth=0,alpha=0.3)
	plt.step(tbb4[:-1],tbh4,'k',linestyle='-.',linewidth=1,where='mid')
	plt.bar(tbb4[:-1],tbh4,width=tbw4,linewidth=0,alpha=0.3)

	plt.xlabel("Time [ns]")
	plt.ylabel("Relative Intensity")
	plt.legend([f1, f2, f3, f4])

	plt.figure(2)
	plt.step(cbb1[:-1],cbh1,'k',linestyle='-',linewidth=1,where='mid')
	plt.bar(cbb1[:-1],cbh1,width=cbw1,linewidth=0,alpha=0.3)
	plt.step(cbb2[:-1],cbh2,'k',linestyle='--',linewidth=1,where='mid')
	plt.bar(cbb2[:-1],cbh2,width=cbw2,linewidth=0,alpha=0.3)
	plt.step(cbb3[:-1],cbh3,'k',linestyle=':',linewidth=1,where='mid')
	plt.bar(cbb3[:-1],cbh3,width=cbw3,linewidth=0,alpha=0.3)
	plt.step(cbb4[:-1],cbh4,'k',linestyle='-.',linewidth=1,where='mid')
	plt.bar(cbb4[:-1],cbh4,width=cbw4,linewidth=0,alpha=0.3)

	plt.xlabel("Charge [pC]")
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

	plot(dt1, dt2, dt3, dt4, ch1, ch2, ch3, ch4)

	#prompt = raw_input("dt or charge?:").lower()

	#if prompt == "dt":
	#	print "executing dt analysis"
	#	plot_t(dt1, dt2, dt3, dt4)
	#elif prompt == "charge":
	#	print "executing charge analysis"
	#	#plot(ch1, ch2, ch3, ch4)
	#else:
	#	print "improper input, exiting..."
	#	sys.exit()
