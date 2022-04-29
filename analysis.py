import os
import ROOT
import sys
#from subprocess import Popen, PIPE, call
import subprocess
import argparse
import psycopg2

# TO-DO, make inputs
OUTPUT = "/home/www/pmt_characterization_website/app/static/images/"

def connect_to_db():

    conn = psycopg2.connect('host=%s dbname=%s user=%s password=%s' % \
                            ('localhost', 'pmt_testing', 'postgres', 'b33feroni'))

    cursor = conn.cursor()
    return conn, cursor


def write_to_db(source, pmtid, pmt_type, hv,
                tts, lp, ap, pp, dr,
                q_peak, q_width, q_high, q_pv):


    conn, cursor = connect_to_db()

    cursor.execute("INSERT INTO pmt_information "
                   "(source, pmt_id, pmt_type, high_voltage, tts_sigma, late_pulsing_pct, after_pulsing_pct, "
                   "pre_pulsing_pct, dark_rate, charge_peak, charge_width, high_charge_pct, charge_peak_to_valley) "
                   "VALUES ('%s', '%s', '%s', %d, %f, %f, %f, %f, %f, %f, %f, %f, %f)" % \
                   (source, pmtid, pmt_type, hv, tts, lp, ap, pp, dr, q_peak, q_width, q_high, q_pv))

    conn.commit()


def fit_charge(hq):

    c = ROOT.TCanvas("c", "c", 800, 600)

    hq.GetXaxis().SetRangeUser(0.3, 4.0)

    m1 = hq.GetMaximum()
    b1 = hq.GetMaximumBin()
    c1 = hq.GetBinCenter(b1)
    fit_range = 0.5

    fit = ROOT.TF1("gaus", "gaus", c1 - fit_range, c1 + fit_range)
    hq.Fit(fit, "Q0", "", c1 - fit_range, c1 + fit_range)

    hq.Draw("")
    fit.Draw("same")

    qmean = fit.GetParameter(1)
    qsigma = fit.GetParameter(2)

    hq.GetXaxis().SetRangeUser(-0.3, 4.0)
    hq.GetYaxis().SetRangeUser(0.0, m1*1.25)

    c.Print("charge.png")

    b1 = hq.FindBin(qmean - 2*qsigma)
    b2 = hq.FindBin(qmean + 2*qsigma)
    b3 = hq.FindBin(10.0)

    total_int = hq.Integral(b1, b3)
    high_charge_int = hq.Integral(b2, b3)

    high_charge_pct = (high_charge_int/total_int)*100

    print ("Charge mean: %.2f pC" % qmean)
    print ("Charge width: %.2f pC" % qsigma)
    print ("High charge rate: %.2f pct" % high_charge_pct)

    return qmean, qsigma, high_charge_pct


def fit_timing(ht):

    c = ROOT.TCanvas("c", "c", 800, 600)

    m1 = ht.GetMaximumBin()
    c1 = ht.GetBinCenter(m1)
    fit_range = 2.0

    fit = ROOT.TF1("gaus", "gaus", c1 - fit_range, c1 + fit_range)
    ht.Fit(fit, "Q0", "", c1 - fit_range, c1 + fit_range)

    tts = fit.GetParameter(2)

    dark_fit = ROOT.TF1("pol0", "pol0", c1 - 60.0, c1 - 10.0)
    ht.Fit(dark_fit, "Q0", "", c1 - 60.0, c1 - 10.0)

    p0 = dark_fit.GetParameter(0)
    # To-do get bin-width automatically
    dark_rate = p0/(0.1*1e-9*ht.GetEntries())

    c_late_low = c1 + 10.0
    c_late_high = c1 + 60.0
    b_low = ht.FindBin(c_late_low)
    b_high = ht.FindBin(c_late_high)
    late_int = ht.Integral(b_low, b_high)

    ht.Draw("")

    fit.Draw("same")
    dark_fit.Draw("same")

    ht.GetXaxis().SetRangeUser(c1 - 10.0, c1 + 10.0)

    c.Update()

    c.Print("time_zoomed.png")

    ht.GetXaxis().SetRangeUser(0.0, 150.0)

    c.SetLogy()

    c.Update()

    c.Print("time.png")

    print ("TTS: %.2f ns" % tts)
    print ("Dark rate: %.2f Hz" % dark_rate)

    return tts, dark_rate


def open_tree(fname, threshold):

    f = ROOT.TFile.Open(fname)
    t = f.Get("output")

    ht = ROOT.TH1D("time","time",2000,-50,150)
    hq = ROOT.TH1D("charge","charge",600,-0.5,5.5)

    print "Entries:", t.GetEntries()
    for i in range(t.GetEntries()):

        t.GetEntry(i)

        hq.Fill(t.charge - t.charge_empty)

        if(t.peak_voltage > threshold): continue

        ht.Fill(t.deltat)

    return ht, hq


def run_analysis(datafile, output_name, pedestal):

    # TO-DO, make inputs
    command = ("/data/snoplus/home/tannerbk/pmt_characterization/src/run_pmt_characterization %s %s "
               "lappd_0 gr0 ch1 gr0 ch0 gr0 ch2 %d" % (datafile, output_name, pedestal))

    print ("Running:", command)

    commands = command.split()

    subprocess.call(commands)


def create_event_file(directory, ofile):

    event_file = open(ofile, "w")

    for f in sorted(os.listdir(directory)):
        if ".h5" not in f: continue
        if f[-3:] != ".h5": continue
        event_file.write(directory + "/" + f + "\n")

    event_file.close()
    return event_file


def write_root_file(output, ht, hq):

    f = ROOT.TFile.Open(output, "RECREATE")
    ht.Write()
    hq.Write()
    f.Write()


def pretty_plot(h, xname):

    h.SetTitle("")
    h.GetYaxis().SetTitle("Counts")
    h.GetXaxis().SetTitle(xname)
    h.SetLineColor(ROOT.kBlack)
    h.SetMarkerColor(ROOT.kBlack)
    h.GetXaxis().SetLabelFont(132)
    h.GetXaxis().SetTitleFont(132)
    h.GetYaxis().SetLabelFont(132)
    h.GetYaxis().SetTitleFont(132)
    h.SetStats(0)

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--source', type=str, required=True)
    parser.add_argument('-d', '--directory', type=str, required=True)
    parser.add_argument('-v', '--high-voltage', type=int, required=True)
    parser.add_argument('-i', '--pmt-id', type=str, required=True)
    parser.add_argument('-p', '--pmt-type', type=str, required=True)
    parser.add_argument('-w', '--pedestal', type=int, default=200)
    parser.add_argument('-t', '--threshold', type=float, default=-5.0)
    parser.add_argument('-f', '--txt-file', type=str, default="data.txt")
    parser.add_argument('-o', '--root-file', type=str, default="data.root")
    args = parser.parse_args()

    source_options = ["Cherenkov", "LED"]

    if args.source not in source_options:
        print ("Invalid source.")
        print ("Options:", source_options)
        sys.exit(1)

    pmt_options = ['R7081', 'R11780', 'R14688']

    if args.pmt_type not in pmt_options:
        print ("Invalid pmt type.")
        print ("Options:", pmt_options)
        sys.exit(1) 
       

    output_name = args.pmt_id + "_" + str(args.high_voltage) + "V"
    dirname = OUTPUT + output_name

    try:
        os.mkdir(dirname)
    except OSError:
        pass

    os.chdir(dirname)

    event_file = create_event_file(args.directory, args.txt_file)

    datafile = dirname + "/" + args.txt_file

    run_analysis(datafile, output_name, args.pedestal)

    root_file = dirname + "/" + output_name + "_lappd_0_gr0_ch1.root"  

    ht, hq = open_tree(root_file, args.threshold)

    pretty_plot(ht, "Time (ns)")
    pretty_plot(hq, "Charge (pC)")

    print ("PMT ID %s" % args.pmt_id)
    print ("HV: %d V" % args.high_voltage)

    tts, dark_rate = fit_timing(ht)

    q_mean, q_width, high_charge_pct = fit_charge(hq)

    # FIXME..
    write_to_db(args.source, args.pmt_id, args.pmt_type, args.high_voltage, tts, \
                5.0, 1.0, 1.0, dark_rate, q_mean, q_width, high_charge_pct, 2.5)

    write_root_file(args.root_file, ht, hq)

