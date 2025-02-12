import ROOT

def open_tree(tree):
    '''
    Loop over the entries in a ROOT tree and fill a histogram
    with the transit time, which is referred to as "deltat" in
    the data files.
    '''

    # Create the ROOT histogram
    h = ROOT.TH1D("","",1000,-100,100)

    # Loop over the entries in the ROOT tree.
    # There is one entry per triggered event.
    for i in range(tree.GetEntries()):

        tree.GetEntry(i)

        # Select waveforms where the PMT pulse drops
        # below a 5 mV threshold.
        if(tree.peak_voltage > -5.0):
            continue

        # Fill the histogram with the value of "deltat"
        # This is our histogram of the transit time.
        h.Fill(tree.deltat)

    return h


def fit_data(h):
    '''
    Lets fit the ROOT histogram with a Gaussian
    '''
    # First we find where the peak is
    m = h.GetMaximumBin()
    c = h.GetBinCenter(m)

    # Then we fit a Gaussian around the peak
    fit_range = 0.5
    fit = ROOT.TF1("gaus", "gaus", c - fit_range, c + fit_range)
    h.Fit(fit, "Q0", "", c - fit_range, c + fit_range)

    # Return the results of the fit
    return fit


if __name__=='__main__':

    # These variables just point to one data file. Better to make these 
    # a user input (eg, sys.argv[1])
    odir = "/home/www/pmt_characterization_website/app/static/images/"
    directory = odir + "ZC2973_Cherenkov_1700V_None_-40.0mV_10.0pC_0.0Hrs_tts_reproducibility_0/"
    filename = "ZC2973_Cherenkov_1700V_None_-40.0mV_10.0pC_0.0Hrs_tts_reproducibility_0_lappd_0_gr0_ch1.root"

    # Open the ROOT file
    f = ROOT.TFile.Open(directory+filename)
    # Open the ROOT tree
    tree = f.Get("output")

    h = open_tree(tree)
    fit = fit_data(h)

    # We extract the width of the fit
    sigma = fit.GetParameter(2)
    print ("The TTS (sigma of the Gaussian) is %.3f ns" % (sigma))

    # Draw the histogram on a ROOT canvas
    c1 = ROOT.TCanvas("c1","c1",800,600)

    # Lets make the plot a little prettier
    h.GetXaxis().SetTitleFont(132)
    h.GetXaxis().SetLabelFont(132)
    h.GetYaxis().SetTitleFont(132)
    h.GetYaxis().SetLabelFont(132)
    h.GetXaxis().SetTitle("#Deltat (ns)")
    h.GetYaxis().SetTitle("Counts")

    # Manually set the x-axis range 
    h.GetXaxis().SetRangeUser(10, 20)

    # Draw the histogram and the fit
    h.Draw("E") # "E" draws with uncertainties on each bin
    fit.SetLineColor(ROOT.kRed)
    fit.Draw("same") # "same" draws on the same canvas

    c1.Update()

    raw_input()

