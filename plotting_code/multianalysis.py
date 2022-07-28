# this is python code was developed by mduce so multiple pmt datasets can be
# compared after they have been run through analysis.py
# contact mduce@gatech.edu if needed
# inputs:
#	-d directory containing all .root files desired for plotting
#	-p for histograms
#	-f for histograms with fits and fit stats
#	-s for deltat histogram where plot peaks are centered on 0
#	-x to save everything that is plotted/printed
#	userinput: identifier for each dataset for plot legend
# output: selected -p, -f, or -s plots (or all)
# 
# NOTE: code requires input of a directory with the to-be-compared .root files
# bring a fork

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
	return value_list, files
	# list[1][1] = dt1, list[n] = dtn, chn

def normalize(data, crop_min, crop_max, bw):
	# here by tygers
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

	# spaghet
	for sets in range(len(x[1])):
		pk_indx = int(x[6][sets]), int(x[7][sets])

		t_adjust = int(round(t_range*pk_indx[0]))
		ch_adjust = int(round(c_range*pk_indx[1]))

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
	# x = tuple of arrays
	# x[0] = t_bin_boundaries, x[1] = t_bin_heights, x[2] = t_bin_width
	# x[3] = c_bin_boundaries, x[4] = c_bin_heights, x[4] = c_bin_width

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
	# x = tuple of arrays
	# x[0] = t_cent_boundaries, x[1] = t_bin_height, x[2] = t_bin_width
	# x[3] = c_cent_boundaries, x[4] = c_bin_height, x[5] = c_bin_width

	t_popt, t_pcov  = [], []
	c_popt, c_pcov = [], []
	t_std, c_std = [], []
	t_mean, c_mean = [], []
	dt_pcov, dc_pcov = [], []
	t_std_err, c_std_err = [], []

	for set in range(len(x[0])):
		# generate approximate mean, std for future fit
		t_mean_out, t_std_out = norm.fit(x[0][set])
		c_mean_out, c_std_out = norm.fit(x[3][set])
		
		# fit data to Gauss func (custom) using scipy curvefit func
		# popt = predicted optimatal parameters (max, mean, std)
		# pcov = covariance array of parameters
		t_popt_out, t_pcov_out = curve_fit(Gauss, x[0][set], x[1][set], \
		 p0 = [1, t_mean_out, t_std_out]) #p0 = best guesses from norm.fit
		c_popt_out, c_pcov_out = curve_fit(Gauss, x[3][set], x[4][set], \
		 p0 = [1, c_mean_out, c_std_out])

		t_popt.append(t_popt_out)
		t_pcov.append(t_pcov_out)

		c_pcov.append(c_pcov_out)
		c_popt.append(c_popt_out)

		t_std.append(t_popt_out[2])
		c_std.append(c_popt_out[2])

		t_mean.append(t_popt_out[1])
		c_mean.append(c_popt_out[1])

		dt_pcov.append(np.diagonal(t_pcov_out))
		dc_pcov.append(np.diagonal(c_pcov_out))

	t_popt, t_pcov, c_popt, c_pcov, t_mean, c_mean, t_std, c_std, dt_pcov, dc_pcov

	for sets in range(len(dt_pcov)):
		t_std_err.append(dt_pcov[sets][2])
		c_std_err.append(dc_pcov[sets][2])
	t_std_err, c_std_err
	
	# caculate fwhm from sigma = fwhm*2.355
	t_fwhm = [sigma*2.355 for sigma in t_std]
	c_fwhm = [sigma*2.355 for sigma in c_std]

	dt_fwhm = [sigma_err*2.355 for sigma_err in t_std_err]
	dc_fwhm = [sigma_err*2.355 for sigma_err in c_std_err]

	return t_popt, t_pcov, c_popt, c_pcov,\
	 t_fwhm, c_fwhm, dt_pcov, dc_pcov,\
	t_mean, c_mean, t_std, c_std,\
	dt_fwhm, dc_fwhm

def plot(x,labels):
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
	color = iter(cm.nipy_spectral(np.linspace(0, 1, 2*len(x[0]))))

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

	return fig1, fig2

def plot_shifted_t(x,labels):
	# x = cropped_t_bin_boundary, cropped_t_bin_height, t_bin_width
	#	cropped_c_bin_boundary, cropped_c_bin_height, c_bin_width
	# 	t_max_index, c_max_index

	fig3, ax3 = plt.subplots()
	color = iter(cm.nipy_spectral(np.linspace(0, 1, 2*len(x[0]))))

	for set in range(len(x[0])):
		c = next(color)

		# subtract dataset from t_max
		ax3.step(x[0][set][:-1] - x[0][set][x[6][set]],\
		x[1][set], where = 'mid',\
		label='%s data' % set, c = c),

	ax3.legend(labels)
	ax3.set_xlabel("time [ns]")
	ax3.set_ylabel("Relative Intnsity")
	fig3.suptitle("Shifted Delta t")

	return fig3

def plot_fits(x, x_full, stats, labels):
	# x = centered vals from cropped dataset
	# x_full = centered vals from full dataset
	# stats = output from apply_gauss_fit
	# labels = user defined labels

	# initalize as many subplots as there are datasets
	fig4, ax4 = plt.subplots(len(x[0]))
	fig5, ax5 = plt.subplots(len(x[0]))

	t_popt = stats[0]
	c_popt = stats[2]

	# colormap to be used
	color = iter(cm.nipy_spectral(np.linspace(0, 1, 2*len(x[0]))))
	
	for set in range(len(x[0])):
		c = next(color)
		
		# full tts plot
		ax4[set].step(x_full[0][set][:-1],x_full[1][set], where='mid', \
		label='%s data' % set, c = c)
		# plot cropped tts plot (fit area)
		ax4[set].step(x[0][set],x[1][set], c = c, where = 'mid', \
		linewidth=2, label='_nolegend_')
		# plot gauss fit
		ax4[set].plot(x[0][set],Gauss(x[0][set], *t_popt[set]),'r-', \
		label='_nolegend_')

		ax4[set].set_xlabel("Time [ns]")
		ax4[set].set_ylabel("Relative Intnsity")
		ax4[set].set_title(labels[set],loc='left',fontsize=10)

		# ----------------------------------

		# full Q plot
		ax5[set].step(x_full[3][set][:-1],x_full[4][set], where='mid', \
		label='%s data' % set, c = c)
		# plot cropped Q plot (fit area)
		ax5[set].step(x[3][set],x[4][set], c = c, where='mid', \
		linewidth=2, label='_nolegend_')
		# plot gauss fit 
		ax5[set].plot(x[3][set],Gauss(x[3][set],*c_popt[set]),'r-', \
		label = '_nolegend_')

		ax5[set].set_xlabel("Charge [pC]")
		ax5[set].set_ylabel("Relative Intensity")
		ax5[set].set_title(labels[set],loc='left',fontsize=10)

	fig4.suptitle("Zoomed TTS, Fit")
	fig5.suptitle("Zoomed Charge Peak, Fit")

	return fig4, fig5

def printstats(stats, labels, files):
	# stats = t_popt, t_pcov, c_popt, c_pcov, t_fwhm, c_fwhm, dt_pcov,\
	#	dc_pcov, t_mean, c_mean, t_std, c_std, dt_fwhm, dc_fwhm
	# labels = user defined labels
	# files = datasets that were analyzed

	# use PrettyTable package
	w = pretty(['Dataset','Label'])
	t = pretty(['Label','TTS Mean','TTS FWHM','TTS STD'])
	c = pretty(['Label','Q Mean','Q FWHM','Q STD'])
	filenames = [name.rsplit('/',1)[-1] for name in files]

	for sets in range(len(labels)):
		# print table with filenames and user-defined labels
		w.add_row([filenames[sets], labels[sets]])

		# print table with TTS fit stats
		t.add_row([labels[sets],stats[8][sets], stats[4][sets], stats[10][sets]])
		t.add_row(["Parameter Error", stats[6][sets][1],\
		 stats[12][sets], stats[6][sets][2]])
			
		# print table with Q fit stats
		c.add_row([labels[sets],stats[9][sets], stats[5][sets], stats[11][sets]])
		c.add_row(["Parameter Error", stats[7][sets][1],\
		stats[13][sets], stats[7][sets][2]])
	t, c, w

	print w
	print t
	print c

	wtable = w.get_string()
	ttable = t.get_string()
	ctable = c.get_string()

	return wtable, ttable, ctable

if __name__=='__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('-d', '--directory', type=str, default="")
	parser.add_argument('-s', '--shift', action='store_true')
	parser.add_argument('-f', '--fits', action='store_true')
	parser.add_argument('-p', '--plots', action='store_true')
	parser.add_argument('-x', '--save', action='store_true')
	args = parser.parse_args()

	# pull .root files from directory
	if args.directory != "":
		direc = args.directory
		data, files = read_dir(direc)
	
		labels = []
		filenames = [name.rsplit('/',1)[-1] for name in files]

		print filenames

		for sets in range(len(files)):
			labels.append(raw_input("Input legend label for dataset: "))
		labels
		

		# build histogram vals for each set
		hist_vals = gethistvals(data)

	labels
	
	if args.plots:
		fig1, fig2 = plot(hist_vals,labels)

		if args.save:
			fig1.savefig(direc + "/Delta_t.png")
			fig2.savefig(direc + "/Charge_Peak.png")

	if args.fits:
		# .20 for tts and .40 charg
		fit_hist_vals = getfithistvals(hist_vals,.18,.30)
		cent_vals = center_boundaries(fit_hist_vals)
		stats = apply_gauss_fit(cent_vals)

		fig4, fig5 = plot_fits(cent_vals,center_boundaries(hist_vals),stats,labels)

		w, t, c = printstats(stats, labels, files)

		if args.save:
			fig4.savefig(direc + "/TTS_Fit.png")
			fig5.savefig(direc + "/Q_fit.png")

			with open(direc + '/fit_stats.txt','wb') as f:
				f.write(w)
				f.write(t)
				f.write(c)

	if args.shift:
		fig3 = plot_shifted_t(hist_vals,labels)

		if args.save:
			fig3.savefig(direc + "/Shifted_t.png")


	plt.show()







