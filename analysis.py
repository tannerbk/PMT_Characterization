import os
import ROOT
import sys
import subprocess
import settings
import psycopg2


def connect_to_db():
    '''
    Connect to the PMT testing database.
    '''
    conn = psycopg2.connect('host=%s dbname=%s user=%s password=%s' % \
                            (settings.dbhost, settings.dbname, settings.dbuser, settings.dbpass))

    cursor = conn.cursor()
    return conn, cursor


def check_db(source, pmtid, pmt_type, hv, comp, comment, settle):
    '''
    Check whether an entry already exists in the database
    '''
    conn, cursor = connect_to_db()

    if not comment:
        cursor.execute("SELECT key FROM pmt_information WHERE source='%s' AND "
                       "pmt_id='%s' AND pmt_type='%s' AND high_voltage=%s AND "
                       "magnetic_compensation='%s' AND settling_time=%s" % \
                       (source, pmtid, pmt_type, hv, comp, settle))

    try:
        return int(cursor.fetchone()[0]) 
    except Exception:
        return None


def update_db(key, source, pmtid, pmt_type, hv,
              tts, lp, ap, pp, dr,
              q_peak, q_width, q_high, q_pv,
              entries, thresh, cr, comp, comment,
              tq_cut, t_thresh, settle, tts_err, dr_err, odir):
    '''
    Write to PMT testing information to the database.
    '''
    conn, cursor = connect_to_db()

    mc = comp.replace("_", " ")

    # Insert the data into the database
    cursor.execute("UPDATE pmt_information SET source='%s', pmt_id='%s', pmt_type='%s', high_voltage=%s, tts_sigma=%s, late_pulsing_pct=%s, after_pulsing_pct=%s, pre_pulsing_pct=%s, dark_rate=%s, charge_peak=%s, charge_width=%s, high_charge_pct=%s, charge_peak_to_valley=%s, entries=%s, threshold=%s, coincidence_rate=%s, magnetic_compensation='%s', comment='%s', trigger_q_cut=%s, trigger_threshold=%s, settling_time=%s, tts_sigma_err=%s, dark_rate_err=%s, directory='%s' WHERE key=%s" % (source, pmtid, pmt_type, hv, tts, lp, ap, pp, dr, q_peak, q_width, q_high, q_pv, entries, thresh, cr, mc, comment, tq_cut, t_thresh, settle, tts_err, dr_err, odir, key))

    conn.commit()


def write_to_db(source, pmtid, pmt_type, hv,
                tts, lp, ap, pp, dr,
                q_peak, q_width, q_high, q_pv,
                entries, thresh, cr, comp, comment,
                tq_cut, t_thresh, settle, tts_err, dr_err, odir):
    '''
    Write to PMT testing information to the database.
    '''
    conn, cursor = connect_to_db()

    mc = comp.replace("_", " ")

    # Insert the data into the database
    cursor.execute("INSERT INTO pmt_information "
                   "(source, pmt_id, pmt_type, high_voltage, tts_sigma, late_pulsing_pct, after_pulsing_pct, "
                   "pre_pulsing_pct, dark_rate, charge_peak, charge_width, high_charge_pct, "
                   "charge_peak_to_valley, entries, threshold, coincidence_rate, magnetic_compensation, " 
                   "comment, trigger_q_cut, trigger_threshold, settling_time, tts_sigma_err, dark_rate_err, directory)"
                   "VALUES ('%s', '%s', '%s', %d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %d, %f, %f, '%s', "
                   "'%s', %f, %f, %f, %f, %f, '%s')" % (source, pmtid, pmt_type, hv, tts, lp, ap, pp, dr, q_peak, \
                   q_width, q_high, q_pv, entries, thresh, cr, mc, comment, tq_cut, t_thresh, settle, tts_err, \
                   dr_err, odir))

    conn.commit()


def fit_charge_led(hq):
    '''
    Fit the charge distribution to extract relevant parameters.
    '''
    c = ROOT.TCanvas("c", "c", 800, 600)

    hq.Rebin(50)

    mmax = hq.GetMaximum()
    b1 = hq.GetMaximumBin()
    c1 = hq.GetBinCenter(b1)
    fit_range = c1*0.2

    fit = ROOT.TF1("gaus", "gaus", c1 - fit_range, c1 + fit_range)
    hq.Fit(fit, "Q0", "", c1 - fit_range, c1 + fit_range)

    hq.Draw("")
    fit.Draw("same")

    qmean = fit.GetParameter(1)
    qsigma = fit.GetParameter(2)

    c.Print("charge.png")

    return qmean, qsigma, 0, 0


def fit_charge(hq, hq_cut, interactive):
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
    hq_cut.SetLineColor(ROOT.kBlue)
    hq_cut.Draw("same")
    fit.Draw("same")

    qmean = fit.GetParameter(1)
    qsigma = fit.GetParameter(2)

    #hq.GetXaxis().SetRangeUser(0.2, 0.6)
    mb = hq.GetXaxis().FindBin(0.5)
    mmin = hq.GetBinContent(mb)

    p_to_v = 0
    if mmin > 0:
        p_to_v = mmax/mmin

    hq.GetXaxis().SetRangeUser(-0.3, 7.0)
    hq.GetYaxis().SetRangeUser(0.0, mmax*1.25)

    b1 = hq.FindBin(qmean - 3*qsigma)
    b2 = hq.FindBin(qmean + 3*qsigma)
    b3 = hq.FindBin(10.0)

    total_int = hq.Integral(b1, b3)
    high_charge_int = hq.Integral(b2, b3)

    high_charge_pct = (high_charge_int/total_int)*100

    print ("Charge mean: %.2f pC" % qmean)
    print ("Charge width: %.2f pC" % qsigma)
    print ("High charge rate: %.2f pct" % high_charge_pct)
    print ("Peak to valley: %.2f" % p_to_v)

    c.Update()

    if interactive:
        print ("Hit enter to continue.")
        raw_input()
    else:
        c.Print("charge.png")

    return qmean, qsigma, high_charge_pct, p_to_v


def fit_timing(ht, entries, interactive):
    '''
    Fit the PMT timing histogram and extract relevant paramters.
    '''
    c = ROOT.TCanvas("c", "c", 800, 600)

    m1 = ht.GetMaximumBin()
    c1 = ht.GetBinCenter(m1)
    fit_range = 0.5

    fit = ROOT.TF1("gaus", "gaus", c1 - fit_range, c1 + fit_range)
    ht.Fit(fit, "Q0", "", c1 - fit_range, c1 + fit_range)

    tts = fit.GetParameter(2)
    tts_unc = fit.GetParError(2)

    df_low = 10
    df_high = 60

    dark_fit = ROOT.TF1("pol0", "pol0", df_low, df_high)
    ht.Fit(dark_fit, "LQ0", "", df_low, df_high)

    p0 = dark_fit.GetParameter(0)
    p0_err = dark_fit.GetParError(0)
    # To-do get bin-width automatically
    dark_rate = p0/(0.1*1e-9*entries) # Bins are 0.1ns wide
    dark_rate_error = p0_err/(0.1*1e-9*entries)

    ht.GetXaxis().SetRangeUser(df_low, df_high)

    ht.GetXaxis().SetRangeUser(-50, 150)

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

    if interactive:
        print ("Hit enter to continue.")
        raw_input()
    else:
        c.Print("time_zoomed.png")

    m1 = ht.GetMaximum()
    ht.GetXaxis().SetRangeUser(0.0, 140.0)
    ht.GetYaxis().SetRangeUser(1e-4, 1.2*m1)

    c.SetLogy()

    c.Update()

    if interactive:
        print ("Hit enter to continue.")
        raw_input()
    else:
        c.Print("time.png")

    print ("TTS (sigma): %.2f ns +/- %.2f") % (tts, tts_unc)
    print ("TTS (FWHM): %.2f ns" % (float(tts)*2.355))
    print ("Dark rate: %.1f Hz" % dark_rate)
    print ("Late fraction: %.2f pct" % fr_late)

    return tts, dark_rate, fr_late, tts_unc, dark_rate_error


def open_tree(fname, threshold, trigger_threshold, trigger_q_cut, source):
    '''
    Open the processed .root files and extract the timing/charge histograms
    '''
    f = ROOT.TFile.Open(fname)
    t = f.Get("output")

    ht = ROOT.TH1D("time","time",2000,-50,150)
    hq = ROOT.TH1D("charge","charge",6000,-0.5,119.5)
    hq_cut = ROOT.TH1D("charge_cut","charge_cut",6000,-0.5,119.5)

    ht.SetDirectory(0)
    hq.SetDirectory(0)
    hq_cut.SetDirectory(0)

    print "Processing", t.GetEntries(), "events"
    entries = 0
    coincidence_rate = 0.0
    for i in range(t.GetEntries()):

        t.GetEntry(i)

        # Bad pedestal window
        if(t.stddev > 0.04): continue

        q = t.charge - t.charge_empty

        hq.Fill(q)

        if(source == "Cherenkov" and t.trigger_charge < trigger_q_cut): continue

        entries += 1

        #if(t.peak_voltage > threshold): continue
        #if(t.peak_voltage_trigger > trigger_threshold): continue
        if(q > 3.0 or q < 0.3): continue
        #if(t.samples_above_threshold < 3): continue

        coincidence_rate += 1.0

        hq_cut.Fill(t.charge - t.charge_empty)
        ht.Fill(t.deltat)

    coincidence_rate /= float(t.GetEntries())
    print ("Coincidence rate after analysis cuts: %.2f pct" % (coincidence_rate*100))

    return ht, hq, hq_cut, entries, coincidence_rate


def run_analysis(wd, datafile, output_name, pedestal, source, channel):
    '''
    Run the analysis code over the data files
    '''
    if source == "Cherenkov": led = 0
    if source == "LED": led = 1

    # The command to run the C++ analysis code
    command = ("%s/src/run_pmt_characterization %s %s %s gr0 %s gr0 ch0 gr0 ch2 %d %d" % \
              (wd, datafile, output_name, settings.digit_name, channel, pedestal, led))

    print ("Running: %s" % command)

    commands = command.split()

    # Run the analysis code
    subprocess.call(commands)


def create_event_file(directory, ofile, max_files):
    '''
    Write the list of hdf5 data files into a file to process.
    '''
    event_file = open(ofile, "w")

    count = 0
    for f in sorted(os.listdir(directory)):
        if ".h5" not in f: continue
        if f[-3:] != ".h5": continue
        if max_files != 0 and count > max_files: continue
        event_file.write(directory + "/" + f + "\n")
        count+=1

    event_file.close()
    return event_file


def write_root_file(output, ht, hq, hq_cut):
    '''
    Write the output root file
    '''
    f = ROOT.TFile.Open(output, "RECREATE")
    ht.Write()
    hq.Write()
    hq_cut.Write()
    f.Write()


def pretty_plot(h, xname):
    '''
    Beautify the root histograms
    '''
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

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--source', type=str, required=True)
    parser.add_argument('-d', '--directory', type=str, required=True)
    parser.add_argument('-v', '--high-voltage', type=int, required=True)
    parser.add_argument('-i', '--pmt-id', type=str, required=True)
    parser.add_argument('-p', '--pmt-type', type=str, required=True)
    parser.add_argument('-c', '--magnetic-compensation', type=str, required=True)
    parser.add_argument('-k', '--channel', type=str, required=True)
    parser.add_argument('-n', '--note', type=str, default="")
    parser.add_argument('-m', '--settle-time', type=float, default=0.0)
    parser.add_argument('-q', '--trigger-q-cut', type=float, default=10.0)
    parser.add_argument('-w', '--pedestal', type=int, default=200)
    parser.add_argument('-t', '--threshold', type=float, default=-5.0)
    parser.add_argument('-r', '--trigger-threshold', type=float, default=-40.0)
    parser.add_argument('-f', '--txt-file', type=str, default="data.txt")
    parser.add_argument('-o', '--root-file', type=str, default="data.root")
    parser.add_argument('-x', '--save', action="store_true")
    parser.add_argument('-y', '--interactive', action="store_true")
    parser.add_argument('-z', '--max_files', type=int, default=0)
    parser.add_argument('-a', '--force-overwrite', action="store_true") 
    args = parser.parse_args()

    # Prepare the string formatting
    if args.source != "LED":
        source = args.source.capitalize()
    else:
        source = args.source
    pmt_id = args.pmt_id.upper()
    pmt_type = args.pmt_type.upper() 
    magnetic_compensation = args.magnetic_compensation.capitalize()

    # Check the source and PMT selected are valid
    if source not in settings.source_options:
        print ("Invalid source.")
        print ("Options:", settings.source_options)
        sys.exit(1)

    if pmt_type not in settings.pmt_options:
        print ("Invalid pmt type:", pmt_type)
        print ("Options:", settings.pmt_options)
        sys.exit(1) 

    # Setup the output filename
    output_name = pmt_id + "_" + str(args.source) + "_" 
    output_name += str(args.high_voltage) + "V" + "_" + magnetic_compensation
    output_name += "_" + str(args.trigger_threshold) + "mV"
    output_name += "_" + str(args.trigger_q_cut) + "pC"
    output_name += "_" + str(args.settle_time) + "Hrs"
    output_name += "_" + args.note

    wd = os.getcwd()
    update_database = False

    # Create output directory
    dirname = settings.output_dir + output_name
    try:
        os.makedirs(dirname,0777)
    except OSError:
        key = check_db(source, pmt_id, pmt_type, args.high_voltage, magnetic_compensation, args.note, args.settle_time)
        if not args.force_overwrite:
            print ("Error directory already exists. Use -a to force overwrite the data.")
            sys.exit(1)
        elif key and args.save:
            print ("You will also override database entry:", key)
            update_database = True

    print ("PMT ID %s" % pmt_id)
    print ("HV: %d V" % args.high_voltage)

    os.chdir(dirname)

    # Create list of hdf5 datafiles
    event_file = create_event_file(args.directory, args.txt_file, args.max_files)
    datafile = dirname + "/" + args.txt_file

    # Run the analysis code over the datafiles
    run_analysis(wd, datafile, output_name, args.pedestal, source, args.channel)

    # Processed .root file output name
    root_file = dirname + "/" + output_name + "_" + settings.digit_name + \
                "_gr0_" + args.channel + ".root"  

    # Open processed .root file and extract the timing and charge histograms
    ht, hq, hq_cut, entries, coinc_rate = open_tree(root_file, args.threshold, args.trigger_threshold, args.trigger_q_cut, source)

    # Beautify the plots
    pretty_plot(ht, "Time (ns)")
    pretty_plot(hq, "Charge (pC)")
    pretty_plot(hq_cut, "Charge (pC)")

    # Now fit the timing and charge figures, which wedo diggerently for the different sources
    if source != "LED":
        tts, dark_rate, fr_late, tts_err, dark_rate_err = fit_timing(ht, entries, args.interactive)
        q_mean, q_width, high_charge_pct, p_to_v = fit_charge(hq, hq_cut, args.interactive)
    else:
        tts, dark_rate, fr_late, tts_err, dark_rate_err = 0, 0, 0, 0, 0
        q_mean, q_width, high_charge_pct, p_to_v = fit_charge_led(hq)

    # Save information to the database
    if args.save and not update_database:
        print "Inserting into database."
        write_to_db(args.source, pmt_id, pmt_type, args.high_voltage, tts, \
                    fr_late, 0.0, 0.0, dark_rate, q_mean, q_width, high_charge_pct, \
                    p_to_v, entries, args.threshold, coinc_rate, magnetic_compensation, \
                    args.note, args.trigger_q_cut, args.trigger_threshold, args.settle_time, \
                    tts_err, dark_rate_err, args.directory)
    elif args.save and update_database:
        print "Updating database."
        update_db(key, args.source, pmt_id, pmt_type, args.high_voltage, tts, \
                    fr_late, 0.0, 0.0, dark_rate, q_mean, q_width, high_charge_pct, \
                    p_to_v, entries, args.threshold, coinc_rate, magnetic_compensation, \
                    args.note, args.trigger_q_cut, args.trigger_threshold, args.settle_time, \
                    tts_err, dark_rate_err, args.directory)

    write_root_file(args.root_file, ht, hq, hq_cut)

    try:
        os.chmod(dirname + "/" + args.root_file, 0777)
        os.chmod(dirname + "/" + args.txt_file, 0777)
        os.chmod(root_file, 0777) 
        os.chmod(dirname + "/" + "charge.png", 0777)
        if source != "LED":
            os.chmod(dirname + "/" + "time_zoomed.png", 0777)
            os.chmod(dirname + "/" + "time.png", 0777)
            os.chmod(dirname + "/" + "time.png", 0777)
    except OSError:
        pass
