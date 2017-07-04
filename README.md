PMT Characterization
====================

Loops through scope traces and calculates the charge and timing histograms. It writes the histogram to a root file that can be open in TBrowser. The user must specify the length of the pedestal window in time bins. Pulses with dark pulses in the pedestal window can be thrown out using a variance cut. The traces that are identified as 'above noise' for the timing histogram are selected using a hard cut on charge.  
