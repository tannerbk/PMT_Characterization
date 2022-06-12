import os
import ROOT
import sys
import subprocess
import argparse
import psycopg2

OUTPUT = "/home/www/pmt_characterization_website/app/static/images/"

def connect_to_db():
    '''
    Connect to the PMT testing database.
    '''
    conn = psycopg2.connect('host=%s dbname=%s user=%s password=%s' % \
                            ('localhost', 'pmt_testing', 'postgres', 'b33feroni'))

    cursor = conn.cursor()
    return conn, cursor


def write_to_db(source, pmtid, pmt_type, hv,
                tts, lp, ap, pp, dr,
                q_peak, q_width, q_high, q_pv,
                entries, thresh, cr, comp, comment,
                tq_cut, t_thresh, settle):
    '''
    Write to PMT testing information to the database.
    '''
    conn, cursor = connect_to_db()

    cursor.execute("INSERT INTO pmt_information "
                   "(source, pmt_id, pmt_type, high_voltage, tts_sigma, late_pulsing_pct, after_pulsing_pct, "
                   "pre_pulsing_pct, dark_rate, charge_peak, charge_width, high_charge_pct, "
                   "charge_peak_to_valley, entries, threshold, coincidence_rate, magnetic_compensation, " 
                   "comment, trigger_q_cut, trigger_threshold, settling_time)"
                   "VALUES ('%s', '%s', '%s', %d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %d, %f, %f, '%s', '%s', %f, %f, %f)" % \
                   (source, pmtid, pmt_type, hv, tts, lp, ap, pp, dr, \
                    q_peak, q_width, q_high, q_pv, entries, thresh, cr, \
                    comp, comment, tq_cut, t_thresh, settle))

    conn.commit()


def fit_charge(hq):
    '''
    Fit the charge distribution to extract relevant parameters.
    '''
    c = ROOT.TCanvas("c", "c", 800, 600)

    hq.GetXaxis().SetRangeUser(0.3, 4.0)

    mmax = hq.GetMaximum()
    b1 = hq.GetMaximumBin()
    c1 = hq.GetBinCenter(b1)
    fit_range = 0.4

    fit = ROOT.TF1("gaus", "gaus", c1 - fit_range, c1 + fit_range)
    hq.Fit(fit, "Q0", "", c1 - fit_range, c1 + fit_range)

    hq.Draw("")
    fit.Draw("same")

    qmean = fit.GetParameter(1)
    qsigma = fit.GetParameter(2)

    hq.GetXaxis().SetRangeUser(0.2, 0.6)
    mmin = hq.GetMinimum()

    p_to_v = 0
    if mmin > 0:
        p_to_v = mmax/mmin

    hq.GetXaxis().SetRangeUser(-0.3, 7.0)
    hq.GetYaxis().SetRangeUser(0.0, mmax*1.25)

    c.Print("charge.png")

    b1 = hq.FindBin(qmean - 3*qsigma)
    b2 = hq.FindBin(qmean + 3*qsigma)
    b3 = hq.FindBin(10.0)

    total_int = hq.Integral(b1, b3)
    high_charge_int = hq.Integral(b2, b3)

    high_charge_pct = (high_charge_int/total_int)*100

    print ("Charge mean: %.2f pC" % qmean)
    print ("Charge width: %.2f pC" % qsigma)
    print ("High charge rate: %.2f pct" % high_charge_pct)
    print ("Peak to value: %.2f" % p_to_v)

    return qmean, qsigma, high_charge_pct, p_to_v


def fit_timing(ht, entries):

    c = ROOT.TCanvas("c", "c", 800, 600)

    m1 = ht.GetMaximumBin()
    c1 = ht.GetBinCenter(m1)
    fit_range = 0.5

    fit = ROOT.TF1("gaus", "gaus", c1 - fit_range, c1 + fit_range)
    ht.Fit(fit, "Q0", "", c1 - fit_range, c1 + fit_range)

    tts = fit.GetParameter(2)

    df_low = 16
    df_high = 6

    dark_fit = ROOT.TF1("pol0", "pol0", c1 - df_low, c1 - df_high)
    ht.Fit(dark_fit, "LQ0", "", c1 - df_low, c1 - df_high)

    p0 = dark_fit.GetParameter(0)
    # To-do get bin-width automatically
    dark_rate = p0/(0.1*1e-9*entries)

    c_prompt_low = c1 - 4.0
    c_late_low = c1 + 4.0
    c_late_high = c1 + 40.0

    bp_low = ht.FindBin(c_prompt_low)
    b_low = ht.FindBin(c_late_low)
    b_high = ht.FindBin(c_late_high)

    late_int = ht.Integral(b_low, b_high)
    prompt_int = ht.Integral(bp_low, b_low)

    late_int_drc = (c_late_high - c_late_low)*1e-9*dark_rate*entries
    prompt_int_drc = (c_late_low - c_prompt_low)*1e-9*dark_rate*entries

    total_late = (late_int - late_int_drc) 
    total_prompt = (prompt_int - prompt_int_drc)

    fr_late = float(total_late)*100/total_prompt

    ht.Draw("")

    fit.Draw("same")
    dark_fit.Draw("same")

    ht.GetXaxis().SetRangeUser(c1 - 6.0, c1 + 10.0)

    c.Update()

    c.Print("time_zoomed.png")

    ht.GetXaxis().SetRangeUser(-20.0, 100.0)

    c.SetLogy()

    c.Update()

    c.Print("time.png")

    print ("TTS (sigma): %.2f ns" % tts)
    print ("TTS (FWHM): %.2f ns" % (float(tts)*2.355))
    print ("Dark rate: %.1f Hz" % dark_rate)
    print ("Late fraction: %.2f pct" % fr_late)

    return tts, dark_rate, fr_late


def open_tree(fname, threshold, trigger_q_cut):

    f = ROOT.TFile.Open(fname)
    t = f.Get("output")

    ht = ROOT.TH1D("time","time",2000,-50,150)
    hq = ROOT.TH1D("charge","charge",600,-0.5,11.5)

    ht.SetDirectory(0)
    hq.SetDirectory(0)

    print "Entries:", t.GetEntries()
    entries = 0
    coincidence_rate = 0.0
    for i in range(t.GetEntries()):

        t.GetEntry(i)

        # Bad pedestal window
        if(t.stddev > 0.04): continue

        hq.Fill(t.charge - t.charge_empty)

        if(t.trigger_charge < trigger_q_cut): continue

        entries += 1

        if(t.peak_voltage > threshold): continue

        coincidence_rate += 1.0

        ht.Fill(t.deltat)

    coincidence_rate /= float(t.GetEntries())

    return ht, hq, entries, coincidence_rate


def run_analysis(datafile, output_name, pedestal):

    # TO-DO, make inputs
    command = ("/data/snoplus/home/tannerbk/pmt_characterization/src/run_pmt_characterization %s %s "
               "lappd_0 gr0 ch1 gr0 ch0 gr0 ch2 %d" % (datafile, output_name, pedestal))

    print ("Running: %s" % command)

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
    h.GetYaxis().SetTitleOffset(1.2)
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
    parser.add_argument('-c', '--magnetic-compensation', type=str, required=True)
    parser.add_argument('-n', '--note', type=str, default="")
    parser.add_argument('-m', '--settle-time', type=float, default=0.0)
    parser.add_argument('-q', '--trigger-q-cut', type=float, default=10.0)
    parser.add_argument('-w', '--pedestal', type=int, default=200)
    parser.add_argument('-t', '--threshold', type=float, default=-5.0)
    parser.add_argument('-r', '--trigger-threshold', type=float, default=-40.0)
    parser.add_argument('-f', '--txt-file', type=str, default="data.txt")
    parser.add_argument('-o', '--root-file', type=str, default="data.root")
    parser.add_argument('-x', '--save', action="store_true")
    args = parser.parse_args()

    source = args.source.upper()
    pmt_id = args.pmt_id.upper()
    pmt_type = args.pmt_type.upper() 

    source_options = ["CHERENKOV", "LED"]

    if source not in source_options:
        print ("Invalid source.")
        print ("Options:", source_options)
        sys.exit(1)

    pmt_options = ['R7081', 'R11780', 'R14688', 'H11934']

    if pmt_type not in pmt_options:
        print ("Invalid pmt type:", pmt_type)
        print ("Options:", pmt_options)
        sys.exit(1) 
       

    output_name = pmt_id + "_" + str(args.high_voltage) + "V" + "_" + args.magnetic_compensation
    output_name += "_" + str(args.trigger_threshold) + "mV"
    output_name += "_" + str(args.trigger_q_cut) + "pC"
    output_name += "_" + args.note

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

    ht, hq, entries, coinc_rate = open_tree(root_file, args.threshold, args.trigger_q_cut)

    pretty_plot(ht, "Time (ns)")
    pretty_plot(hq, "Charge (pC)")

    print ("PMT ID %s" % pmt_id)
    print ("HV: %d V" % args.high_voltage)

    tts, dark_rate, fr_late = fit_timing(ht, entries)

    q_mean, q_width, high_charge_pct, p_to_v = fit_charge(hq)

    if args.save:
        write_to_db(args.source, pmt_id, pmt_type, args.high_voltage, tts, \
                    fr_late, 1.0, 1.0, dark_rate, q_mean, q_width, high_charge_pct, \
                    p_to_v, entries, args.threshold, coinc_rate, args.magnetic_compensation, \
                    args.note, args.trigger_q_cut, args.trigger_threshold, args.settle_time)

        write_root_file(args.root_file, ht, hq)

