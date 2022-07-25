# this is python code was developed by mduce so multiple pmt datasets can be
# compared after they have been run through analysis.py
# input: a directory containing all .root files desired for plotting
#	userinput: identifier for each dataset for plot legend
# output: layered histogram of four datasets
# 
# NOTE: code requires input of a directory with the to-be-compared .root files

import ROOT # import the CERN ROOT module
import sys
import argparse
from matplotlib import pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
import glob
from prettytable import PrettyTable as pretty
from scipy.stats import norm
from scipy.optimize import curve_fit

def open_tree(tree):
	
	# print total number of recorded events
	print "Number of entries:", tree.GetEntries() 
	
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
#	print "Extracting data from root file..."
	
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
def read_dir(dir_name):
	value_list = []
	files = sorted(glob.glob(dir_name + "/*.root"))
	print "Opening directory ", dir_name
	for file in files:
		dt, ch = open_root(file)
		value_list.append((dt,ch))
		print "Data extracted from ", file
	return value_list, files
	# list[1][1] = dt1, list[n] = dtn, chn

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

def gethistvals(x):
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

	t_max_index = []
	ch_max_index = []
	
	# user plug in crop values if non-standard
	#t_min = raw_input("Lower value of time: ")
	#t_max = raw_input("Upper value of time: ")
	#ch_min = raw_input("Lower value of charge: ")
	#ch_max = raw_input("Upper value of charge: ")

	# standard lower and upper lims of histograms
	t_min = 10
	t_max = 20
	ch_min = 0
	ch_max = 7
	
	bin_width = .1

	# crop and normalize datasets, allocate tuple values
	for e in range(len(x)):
		# normalize and crop data
		t_out = normalize(x[e][0], t_min, t_max, bin_width)
		ch_out = normalize(x[e][1], ch_min, ch_max, bin_width)		

		# allocate normalization outputs
		t_bin_boundary.append(t_out[0])
		t_bin_height.append(t_out[1])
		t_bin_width.append(t_out[2])

		t_max_index.append(get_peak_index(t_out[1]))
		ch_max_index.append(get_peak_index(ch_out[1]))

		ch_bin_boundary.append(ch_out[0])
		ch_bin_height.append(ch_out[1])
		ch_bin_width.append(ch_out[2])

	# t_bin_item is a list of np arrays
	t_bin_boundary, t_bin_height, t_bin_width, \
	ch_bin_boundary, ch_bin_height, ch_bin_width, \
	t_max_index, ch_max_index

	return t_bin_boundary, t_bin_height, t_bin_width,  \
	ch_bin_boundary, ch_bin_height, ch_bin_width, \
	t_max_index, ch_max_index

def getfithistvals(x,t_range,c_range):
	# x = tuple of lists of np arrays
	# x[0] = t_bin_boundary, x[1] = ...height, x[2] = ...width
	# x[3] = ch_bin_boundary, x[4] = ...height, x[5] = ...width
	# x[6] = t_max_index, x[7] = ch_max_index
	
	cropped_t_bin_boundary = []
	cropped_t_bin_height = []
	cropped_ch_bin_boundary = []
	cropped_ch_bin_height = []

	for sets in range(len(x[1])):
		pk_indx = int(x[6][sets]), int(x[7][sets])

		t_adjust = int(t_range*pk_indx[0])
		ch_adjust = int(c_range*pk_indx[1])

		cropped_t_bin_boundary.append(crop_data(x[0][sets],pk_indx[0],t_adjust))
		cropped_t_bin_height.append(crop_data(x[1][sets],pk_indx[0],t_adjust))

		cropped_ch_bin_boundary.append(crop_data(x[3][sets],pk_indx[1],ch_adjust))
		cropped_ch_bin_height.append(crop_data(x[4][sets],pk_indx[1],ch_adjust))

	return cropped_t_bin_boundary, cropped_t_bin_height, x[2], \
	cropped_ch_bin_boundary, cropped_ch_bin_height, x[5]

def get_peak_index(x):
	pk_val_indxs = np.where(x == 1)
	# pk is a tuple of arrays containing the indices at which x == 1
	pk_indx = pk_val_indxs[0]
	return pk_indx

def crop_data(x,pk,adjust):
	cropped_x = x[pk-adjust:pk+adjust]
	return cropped_x

def Gauss(c, a, x, sigma):
	return a * np.exp(-(c-x)**2 / (2 * sigma**2))

def center_boundaries(x):
	
	t_half_bin_width = x[2][0]/2
	ch_half_bin_width = x[5][0]/2

	t_centered_boundaries = []
	ch_centered_boundaries = []

	for sets in range(len(x[0])):
		t_centered_boundaries.append([number - t_half_bin_width for number in x[0][sets]])
		ch_centered_boundaries.append([number - ch_half_bin_width for number in x[3][sets]])
	t_centered_boundaries, ch_centered_boundaries

	return t_centered_boundaries, x[1], x[2], ch_centered_boundaries, x[4], x[5]

def apply_gauss_fit(x):
	t_mean = []
	t_std = []
	t_popt = []
	t_pcov = []
	
	c_mean = []
	c_std = []
	c_popt = []
	c_pcov = []

	for set in range(len(x[0])):
		t_mean_out, t_std_out = norm.fit(x[0][set])
		c_mean_out, c_std_out = norm.fit(x[3][set])
		
		t_mean.append(t_mean_out)
		t_std.append(t_std_out)

		c_mean.append(c_mean_out)
		c_std.append(c_std_out)

		t_popt_out, t_pcov_out = curve_fit(Gauss, x[0][set], x[1][set], \
		 p0 = [1, t_mean_out, t_std_out])
		c_popt_out, c_pcov_out = curve_fit(Gauss, x[3][set], x[4][set], \
		 p0 = [1, c_mean_out, c_std_out])

		t_popt.append(t_popt_out)
		t_pcov.append(t_pcov_out)

		c_pcov.append(c_pcov_out)
		c_popt.append(c_popt_out)

	t_mean, t_std, t_popt, t_pcov, c_mean, c_std, c_popt, c_pcov
	
	t_fwhm = [sigma*2.35 for sigma in t_std]
	c_fwhm = [sigma*2.35 for sigma in c_std]

	return t_mean, t_std, t_popt, t_pcov,\
	 c_mean, c_std, c_popt, c_pcov,\
	 t_fwhm, c_fwhm

def plot(x,labels,direc):
	# input: output from normalization function above where
	# x[0] = t_bin_boundary, x[1] = ...height, x[2] = ...width
	# x[3] = ch_bin_boundary, x[4] = ...height, x[5] = ...width
	# x[6] = t_max_index, x[7] = ch_max_index

	# output: 3 plots with all datasets plotted
	# plot 1 Intensity vs TTS
	# plot 2 Intensity vs Charge Peak

	# initialize plots and labels
	fig1, ax1 = plt.subplots()
	fig2, ax2 = plt.subplots()

	# establish colorgradient
	color = iter(cm.rainbow(np.linspace(0, 1, len(x[0]))))

	# for number of datasets generate plots
	for set in range(len(x[0])):

		# increment plot color from gradient
		c = next(color)
		
		fig1.suptitle("Delta t")
		ax1.step(x[0][set][:-1],x[1][set],\
		where = 'mid',label='%s data' % set, c = c)
		ax1.set_xlabel("time [ns]")
		ax2.set_ylabel("Relative Intensity")

		fig2.suptitle("Charge Peak")
		ax2.step(x[3][set][:-1],x[4][set],\
		where = 'mid',label='%s data' % set, c = c)
		ax2.set_xlabel("Charge [pC]")
		ax2.set_ylabel("Relative Intensity")

	ax1.legend(labels)
	ax2.legend(labels)

	fig1.savefig(direc + "/Delta_t.png")
	fig2.savefig(direc + "/Charge_Peak.png")

def plot_shifted_t(x,labels,direc):

	fig3, ax3 = plt.subplots()

	color = iter(cm.rainbow(np.linspace(0, 1, len(x[0]))))

	for set in range(len(x[0])):
		c = next(color)

		ax3.step(x[0][set][:-1]-\
		x[0][set][x[6][set]],\
		x[1][set], where = 'mid',\
		label='%s data' % set, c = c),

	ax3.legend(labels)
	ax3.set_xlabel("time [ns]")
	ax3.set_ylabel("Relative Intnsity")
	fig3.suptitle("Shifted Delta t")

	fig3.savefig(direc + "/Shifted_t.png")

def plot_fits(x, x_full, labels, direc):

	fig4, ax4 = plt.subplots()
	fig5, ax5 = plt.subplots()

	stats = apply_gauss_fit(x)
	t_popt = stats[2]
	c_popt = stats[6]

	color = iter(cm.rainbow(np.linspace(0, 1, len(x[0]))))

	for set in range(len(x[0])):
		c = next(color)

		# full tts plot
		ax4.step(x_full[0][set][:-1],x_full[1][set], where='mid', \
		label='%s data' % set, c = c)
		# plot cropped tts plot (fit area)
		ax4.step(x[0][set],x[1][set], c = c, where = 'mid', \
		linewidth=2, label='_nolegend_')
		# plot gauss fit
		ax4.plot(x[0][set],Gauss(x[0][set], *t_popt[set]),'r-', \
		label='_nolegend_')

		# full Q plot
		ax5.step(x_full[3][set][:-1],x_full[4][set], where='mid', \
		label='%s data' % set, c = c)
		# plot cropped Q plot (fit area)
		ax5.step(x[3][set],x[4][set], c = c, where='mid', \
		linewidth=2, label='_nolegend')
		# plot gauss fit 
		ax5.plot(x[3][set],Gauss(x[3][set],*c_popt[set]),'r-', \
		label = '_nolegend_')

	fig4.suptitle("Zoomed TTS with Gaussian Fit")
	ax4.set_xlabel("time [ns]")
	ax4.set_ylabel("Relative Intnsity")
	ax4.legend(labels)
	fig5.suptitle("Zoomed Charge Peak with Gaussian Fit")
	ax5.set_xlabel("Charge [pC]")
	ax5.set_ylabel("Relative Intensity")
	ax5.legend(labels)

	fig4.savefig(direc + "/TTS_Fit.png")
	fig5.savefig(direc + "/Q_fit.png")

def printstats(x, labels, files, direc):
	
	t = pretty(['Dataset','Label','TTS Mean','TTS FWHM','TTS STD'])
	c = pretty(['Dataset','Label','Q Mean','Q FWHM','Q STD'])
	t.title = 'TTS'
	t.padding_width = 1
	filenames = [name.rsplit('/',1)[-1] for name in files]

	for sets in range(len(labels)):
		t.add_row([filenames[sets], labels[sets], x[0][sets], \
		x[8][sets], x[1][sets]])
	
		c.add_row([filenames[sets], labels[sets], x[4][sets], \
		x[9][sets], x[5][sets]])
	t, c
	print(t)
	print(c)

	ttable = t.get_string()
	ctable = c.get_string()

	with open(direc + '/fit_stats.txt','wb') as f:
		f.write(ttable)
		f.write(ctable)

if __name__=='__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('-d', '--directory', type=str, default="")
	parser.add_argument('-s', '--shift', type=str, default="")
	parser.add_argument('-f', '--fits', type=str, default="")
	parser.add_argument('-p', '--plots',type=str, default="")
	args = parser.parse_args()

	# pull .root files from directory
	if args.directory != "":
		direc = args.directory
		data, files = read_dir(direc)
	
		labels = []

		for sets in range(len(files)):
			labels.append(raw_input("Input legend label for dataset: "))
		labels
		
		# build histogram vals for each set
		hist_vals = gethistvals(data)

	labels
	
	if args.plots == "y":
		plot(hist_vals,labels, direc)

	if args.fits == "y":
		fit_hist_vals = getfithistvals(hist_vals,.20,.40)
		centered_vals = center_boundaries(fit_hist_vals)
		stats = apply_gauss_fit(centered_vals)
		plot_fits(centered_vals,center_boundaries(hist_vals),labels,direc)

		printstats(stats, labels, files, direc)

	if args.shift == "y":
		plot_shifted_t(hist_vals,labels,direc)

	plt.show()







