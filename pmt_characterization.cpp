#include "pmt_characterization.h"
#include "tools.h"

int main (int argc, char* argv[])
{
    // Datafile
    datafile = argv[1];

    // Set the running mode
    // mode = 0: normal
    // mode = 1: debugging output
    // mode = 2: separate pulse types in timing
    // mode = 3: write individual waveforms with triple lognormal fits to file
    if(argc > 4){
       mode = atoi(argv[4]);
       cout << "Setting run mode: " << mode << endl;
    }
    else{
       mode = 0;
    }

    // Fitter for waveform fitting to lognormal pulses
    TVirtualFitter::SetDefaultFitter("Minuit");
    TVirtualFitter::SetMaxIterations(50000);
    TFile *ff = 0;
    if(mode == 3){
      ff = TFile::Open("waveforms.root","new");
    }

    double sampling_time = GetSampling(datafile)*1e9;
    size_t length_trace = GetLength(datafile);


    // Create ROOT histograms
    TH1F *average_waveform = new TH1F("Average_Waveform","",
        length_trace, 0, length_trace*sampling_time);

    TH1F *average_waveform_trigger = new TH1F("Average_Waveform_Trigger","",
        length_trace, 0, length_trace*sampling_time);
    TH1F *Timing = new TH1F("Timing","",(length_trace)*sampling_time*2, 0,
        (length_trace)*sampling_time);
    
    TH1F *prompt_hits = new TH1F("prompt","",(length_trace)*sampling_time*2, 0,
        (length_trace)*sampling_time);
    TH1F *late_hits = new TH1F("late","",(length_trace)*sampling_time, 0,
        (length_trace)*sampling_time);
    TH1F *double_hits_prompt = new TH1F("double_prompt","",(length_trace)*sampling_time, 0,
        (length_trace)*sampling_time);
    TH1F *double_hits_late = new TH1F("double_late","",(length_trace)*sampling_time, 0,
        (length_trace)*sampling_time);
    TH1F *dark_hits = new TH1F("dark","",(length_trace)*sampling_time, 0,
            (length_trace)*sampling_time);

    TH1F *fit_waveform = new TH1F("Fit Waveform","", length_trace, 0, length_trace);

    TH1F * parameter0 = new TH1F("p0","",20,-100,100);
    TH1F * parameter3 = new TH1F("p3","",20,-100,100);
    TH1F * parameter6 = new TH1F("p6","",20,-100,100);

    TH1F * parameter1 = new TH1F("p1","",200,-0.2,0.2);
    TH1F * parameter4 = new TH1F("p4","",200,-0.2,0.2);
    TH1F * parameter7 = new TH1F("p7","",200,-0.2,0.2);

    TH1F * parameter2 = new TH1F("p2","",200,-0.2,5);
    TH1F * parameter5 = new TH1F("p5","",200,-0.2,5);
    TH1F * parameter8 = new TH1F("p8","",200,-0.2,5);

    waveform_voltage.resize(length_trace, 0);
    waveform_voltage_trigger.resize(length_trace,0);

    // Find a starting point for the window using the average waveforms
    int count = GetWindowStart(datafile);

    // Build the average waveform for both PMTs
    for(size_t i = 0; i < length_trace; i++){
      waveform_voltage[i] = waveform_voltage[i]/count;
      average_waveform->SetBinContent(i+1, waveform_voltage[i]);
    }

    for(size_t i = 0; i < length_trace; i++){
      waveform_voltage_trigger[i] = waveform_voltage_trigger[i]/count;
      average_waveform_trigger->SetBinContent(i+1, waveform_voltage_trigger[i]);
    }

    // Calculate good starting points
    unsigned long pre_window = 12 / (dx*1e9);
    unsigned long pre_trigger_window = 20 / (dx*1e9);
    unsigned long window_start = average_waveform->GetMinimumBin() - pre_window;
    unsigned long window_start_trigger = average_waveform_trigger->GetMinimumBin() - pre_trigger_window;
    cout << " ----- Window Start Times ----- " << endl;
    cout << "Integration start time " << window_start*dx*1e9 << " ns." << endl;
    cout << "Trigger start time " << window_start_trigger*dx*1e9 << " ns." << endl;
    size_t number_windows = length_trace/(2*size_signal_window);

    // Size of the pedestal window, in samples
    // Either user set or taken from the average waveform
    if(argc > 3){
      window_width = atoi(argv[3]);
    }
    else{
      window_width = min(window_start, window_start_trigger) - 10/(dx*1e9);
    }

    int file_count = 0;
    int window_count = 0;
    ifstream ifs (argv[1] , ifstream::in); // Open File List Stream
    ifs >> filename;

    while (ifs.good()){ // While Stream Is Open, Analyze New Files.

      file.openFile(filename, H5F_ACC_RDONLY); //Open HDF5 File
      ifs >> filename;

      // Open dataset for characterization PMT
      dataset = file.openDataSet(pmt_channel);
      horiz_interval = dataset.openAttribute("horiz_interval");
      vertical_gain = dataset.openAttribute("vertical_gain");
      horiz_interval.read(PredType::NATIVE_DOUBLE, &dx);
      vertical_gain.read(PredType::NATIVE_DOUBLE, &dy);

      // Open dataset for trigger PMT. Sampling time is same for both channels.
      dataset_trigger = file.openDataSet(trigger_channel);
      vertical_gain_trigger = dataset_trigger.openAttribute("vertical_gain");
      vertical_gain_trigger.read(PredType::NATIVE_DOUBLE, &dy2);

      DataCluster * datacluster = Init_Data(&dataset);
      DataCluster * datacluster_trigger = Init_Data(&dataset_trigger);

      size_t window_length = datacluster->trace_length;
      if(file_count == 0){
          cout << " ----- PMT Characterization Parameters ----- " << endl;
          cout << "The time bin width is " << dx*1e9 << " ns" << endl;
          cout << "The trace length is " << window_length*dx*1e9 << " ns" << endl;
          cout << "The pedestal window is 0 - " << window_width*dx*1e9 << " ns" << endl;
          cout << "The trigger window is " << window_start_trigger*dx*1e9 << " - " << window_start_trigger*dx*1e9 + 40 << " ns" << endl;
          cout << "The pmt window is " << window_start*dx*1e9 << " - " << window_start*dx*1e9 + size_signal_window*dx*1e9 << " ns" << endl;
          cout << "The vertical resolution is " << dy*1000 << " mV" << endl;
          cout << "The vertical resoution of the trigger is " << dy2*1000 << " mV" << endl;
          cout << "Analyzing Files: " << endl;
      }
      cout << "File " << file_count << ": " << filename << endl;

      float variance;
      float max_time_entry = 0.0;
      size_t time_bin_counter = 0;
      size_t time_bin_counter2 = 0;
      size_t time_bin_counter_trigger = 0;
      size_t time_bin_counter2_trigger = 0;
      float prompt_peak = 0;

      for(size_t j = 0; j < datacluster->n_traces; j++){

        Read_Trace(datacluster,j);
        Read_Trace(datacluster_trigger,j);

        float pedestal = TMath::Mean (window_width, datacluster->data_out)*dy;
        float trigger_pedestal = TMath::Mean (window_width, datacluster_trigger->data_out)*dy2;

        float window_counter = 0;
        double ncharge = 0.0;
        variance = 0.0;
        int k = 0;
        int nwindows = 0;
        int found_late_pulse = 0;

        for(size_t i = 0; i < window_width; i++){
          float voltage = ((float)datacluster->data_out[i]*dy-pedestal);
          variance += voltage * voltage;
        }

        // Only cut on variance of test PMT
        // traces that don't pass this cut might have dark pulses in the ped window
        if(variance < variance_cut) {

          for(size_t i = window_start; i < window_start + size_signal_window; i++){
    	    float voltage = ((float)datacluster->data_out[i]*dy-pedestal);
    	    ncharge += (voltage*((-1000.0*dx*1e9)/termination_ohms));
          }

          for(size_t window = 0; window < number_windows; window++){

            // need to initialze these all for every window
            double window_charge = 0.0;
            max_time_entry = 0.0;
            time_bin_counter = 0.0;
            time_bin_counter2 = 0.0;

    	    for(size_t i = window_counter + window_width; i < window_counter + window_width + size_signal_window; i++){
    	      float voltage = ((float)datacluster->data_out[i]*dy-pedestal);
    	      window_charge += (voltage*((-1000.0*dx*1e9)/termination_ohms));
    	    }

    	    if(window_charge > charge_cut_low && window_charge < charge_cut_high){

    	      for(size_t i = window_counter + window_width; i < window_counter + window_width + size_signal_window; i++){
    	        float voltage = ((float)datacluster->data_out[i]*dy-pedestal);
    	        if(voltage < max_time_entry){
    	          max_time_entry = voltage;
    	          time_bin_counter = i;
    	        }
    	      }

    	      prompt_peak = ((float)datacluster->data_out[time_bin_counter]*dy-pedestal);
              peaks->Fill(prompt_peak);

              for(size_t i = time_bin_counter - 20/(dx*1e9); i < time_bin_counter; i++){
    	        float voltage = ((float)datacluster->data_out[i]*dy-pedestal);
    	        if(voltage > prompt_peak*prompt_peak_fraction){
    	          time_bin_counter2 = i;
    	        }
    	      }

              max_time_entry = 0;
              time_bin_counter_trigger = 0;
              time_bin_counter2_trigger = 0;

              for(size_t i = window_start_trigger; i < window_start_trigger + 40/(dx*1e9); i++){
                float voltage = ((float)datacluster_trigger->data_out[i]*dy2-trigger_pedestal);
                if(voltage < max_time_entry){
    	            max_time_entry = voltage;
    	            time_bin_counter_trigger = i;
                }
              }

              prompt_peak = ((float)datacluster_trigger->data_out[time_bin_counter_trigger]*dy2-trigger_pedestal);

              // Look 20ns before peak for 20% crossing
              for(size_t i = time_bin_counter_trigger - 20/(dx*1e9); i < time_bin_counter_trigger; i++){
                float voltage = ((float)datacluster_trigger->data_out[i]*dy2-trigger_pedestal);
                if(voltage > prompt_peak*prompt_peak_fraction){
    	            time_bin_counter2_trigger = i;
                }
              }

    	      pmt_time.push_back(time_bin_counter2*dx*1e9);
              trigger_time.push_back(time_bin_counter2_trigger*dx*1e9);
              double dt = (double(time_bin_counter2) - double(time_bin_counter2_trigger))*dx*1e9;
              kTime.push_back(dt);
              kCharge.push_back(window_charge);
              

              if(mode == 2){
                if(nwindows == 0){
                  prompt_pulse.push_back(dt);
                  prompt_q.push_back(window_charge);
                  nwindows++;
    	          k++;
    	          window_counter = window_counter + size_signal_window;
                  continue;
                }
                if(nwindows == 1){
                  late_pulse.push_back(dt);
                  late_q.push_back(window_charge);
                  found_late_pulse = 1;
                  nwindows++;
                }
              }
            }

            if(mode == 2){
              // Got to the end of the window and didn't find a double pulse
              if(nwindows == 1 && found_late_pulse == 0 && window == number_windows - 1){
                prompt_pulse.pop_back();
                prompt_q.pop_back();
              }
            }
 
    	    k++;
    	    window_counter = window_counter + size_signal_window;
          }

          cout << "Analyzing trace " << j+1 << " out of "
               <<datacluster->n_traces<< " total traces " <<'\r';
          cout.flush();


          if(mode == 3 && ncharge > 0.5){
            //TCanvas *cn = new TCanvas("","",800,600);
            float std = TMath::RMS (window_width, datacluster->data_out);
            for(size_t i = 0; i < window_length; i++){
    	      float voltage = ((float)datacluster->data_out[i]*dy-pedestal);
              fit_waveform->SetBinContent(i, voltage);
              fit_waveform->SetBinError(i, std*dy);
            }
            double xmin = window_width;
            double xmax = (window_width + 30/dx);
            double min_bin = fit_waveform->GetMinimumBin();
            double max_time = fit_waveform->GetBinCenter(min_bin);
            double max_voltage = fit_waveform->GetBinContent(min_bin);
            TF1* waveform3 = new TF1("waveform3",TripleLognormal,xmin,xmax,9);
            double bottom_fit = 0.0;
            double top_fit = 0.0;
            for(size_t i = max_time - 150; i < max_time; i++){
              float voltage = ((float)datacluster->data_out[i]-pedestal);
              if(voltage > max_voltage*0.05){
                bottom_fit = i;
              }
            }
            for(size_t i = max_time + 150; i > max_time; i--){
               float voltage = ((float)datacluster->data_out[i]-pedestal);
               if(voltage > max_voltage*0.15){
                 top_fit = i;
               }
            }
            double bf = bottom_fit - 300;
            double tf = top_fit + 300;
            waveform3->SetParameter(0,max_time);
            waveform3->SetParameter(1,0.02);
            waveform3->SetParameter(2,0.1);
            waveform3->SetParameter(3,max_time);
            waveform3->SetParameter(4,0.05);
            waveform3->SetParameter(5,0.1);
            waveform3->SetParameter(6,max_time);
            waveform3->SetParameter(7,0.05);
            waveform3->SetParameter(8,0.1);
            waveform3->SetLineColor(kRed);
            fit_waveform->Fit("waveform3","Q","",bf,tf);
            double chi2_3 = waveform3->GetChisquare();
            double NDF_3 = waveform3->GetNDF();
            const char* conv3 = gMinuit->fCstatu.Data();
            if(chi2_3/NDF_3 < 2.5 && strcmp(conv3,converged)==0){
              //fit_waveform->Write();
              //cn->Update();
              if(chi2_3/NDF_3 < 2.5 && strcmp(conv3,converged)==0){
                vector<double> parameters(9,0);
                for(int p = 0; p < 9; p++){
                  parameters.at(p) = waveform3->GetParameter(p);
                }
                if(parameters.end() == std::find(parameters.begin(), parameters.end(), true)){
                   parameter0->Fill(parameters[0] - max_time);
                   parameter1->Fill(parameters[1]);
                   parameter2->Fill(parameters[2]);
                   parameter3->Fill(parameters[3] - max_time);
                   parameter4->Fill(parameters[4]);
                   parameter5->Fill(parameters[5]);
                   parameter6->Fill(parameters[6] - max_time);
                   parameter7->Fill(parameters[7]);
                   parameter8->Fill(parameters[8]);
                }
              }
            }
          }


          window_count++; // number of traces
          charges_signal->Fill(ncharge);
          pedestals->Fill(pedestal);
          trigger_pedestals->Fill(trigger_pedestal);
        }
        //store histograms
        variances->Fill(variance);
      }
      file_count++;
      file.close();
    }
    ifs.close();

    if(mode == 1){
      cout  << "The number of Traces accepted is " << window_count << endl;
      cout  << "The number of files is " << file_count  << endl;
    }

    // If this code wasn't bad enough, its about to get worse
    // continue reading at your own peril

    TAxis *axis = charges_signal->GetXaxis();

    // spe fit
    float begin_spe_peak = 0.4;
    float end_spe_peak = 10.0;
    int first_bin = axis->FindBin(begin_spe_peak);
    int last_bin = axis->FindBin(end_spe_peak);
    double charge_finder[last_bin - first_bin];
    int f = 0.0;
    int b = 0.0;
    float max_charge = 0.0;
    float max_charge_entry = 0.0;
    float first_bin2 = first_bin;

    while(first_bin < last_bin){
      charge_finder[f] = charges_signal->GetBinContent(first_bin);
      if(charge_finder[f] > max_charge_entry){
        max_charge_entry = charge_finder[f];
        b = f;
      }
      f++;
      first_bin++;
    }

    max_charge = axis->GetBinCenter(first_bin2 + b);

    Double_t bottom_charge_fit = 2.0/3.0 * max_charge;
    Double_t top_charge_fit = 1.50 * max_charge;

    // looking at the gaussian noise
    Int_t zero_bin = axis->FindBin(0.0);
    Double_t entries_zero_bin = charges_signal->GetBinContent(zero_bin);
    Double_t entrie_checker = entries_zero_bin;
    while(entrie_checker >= entries_zero_bin/2.0){
      zero_bin++;
      entrie_checker = charges_signal->GetBinContent(zero_bin);
    }
    Double_t charge_half_noise_height = charges_signal->GetBinCenter(zero_bin);

    // Fit the charge peak and the charge valley
    TCanvas *c1 = new TCanvas("c1","Charge",200,10,500,500);
    charges_signal->Fit("gaus","q 0","",-charge_half_noise_height, charge_half_noise_height);
    TF1 *myfunc2=charges_signal->GetFunction("gaus");
    Double_t p22=myfunc2->GetParameter(2);
    charges_signal->Fit("pol2","q 0","", 6*p22, bottom_charge_fit);
    TF1 *myfunc3=charges_signal->GetFunction("pol2");
    Double_t p000=myfunc3->GetParameter(0);
    Double_t p111=myfunc3->GetParameter(1);
    Double_t p222=myfunc3->GetParameter(2);

    charges_signal->Fit("gaus","q","",bottom_charge_fit, top_charge_fit);
    TF1 *myfunc = charges_signal->GetFunction("gaus");

    Double_t spe_peak= myfunc->GetParameter(0);
    Double_t spe_mean=myfunc->GetParameter(1);
    Double_t spe_sigma=myfunc->GetParameter(2);
    charges_signal->SetLineColor(kBlack);
    charges_signal->SetStats(0);
    charges_signal->Draw();
    c1->Update();

    // Extremely hacky way of finding the charge "valley"
    Float_t x = 0.0;
    Double_t min_function = 0.0;
    Double_t d_function[20000] = {0.0};
    Int_t s = 0.0;

    while(x < 2.0){
      d_function[s] = p111 + 2*p222*x;
      if(d_function[s-1] < 0.0 && d_function[s] > 0.0){
        min_function = p000 + p111*x + p222*x*x;
      }
      x += 0.001;
      s++;
    }

    // use to count number of entries above 3 times the electronic noise width
    Double_t low_charge_bin = axis->FindBin(3*p22);
    // use to count number of entries 3 sigma above the charge peak
    Double_t high_charge_bin = axis->FindBin(spe_mean + 3*spe_sigma);
    Double_t low_charge_entry = 0.0;
    Double_t high_charge_entry = 0.0;

    for(int i = low_charge_bin; i < charges_signal->GetNbinsX(); i++){
      low_charge_entry += charges_signal->GetBinContent(i);
      if(i > high_charge_bin){
        high_charge_entry += charges_signal->GetBinContent(i);
      }
    }

    Double_t FWHM =  spe_sigma*2*sqrt(2*log(2));
    Double_t HCT = high_charge_entry * 100 / low_charge_entry;
    Double_t elec_width = p22;
    Double_t Peak_to_valley = spe_peak/min_function;
    TPaveStats *ptstats = new TPaveStats(0.6,0.6,0.98,0.98,"brNDC");
    ptstats->SetName("stats");
    ptstats->SetBorderSize(2);
    ptstats->SetTextAlign(12);
    ptstats->SetTextFont(42);
    ptstats->SetTextSize(0.035);
    ptstats->SetShadowColor(0);

    TText *text_t = ptstats->AddText(box_title);
    text_t->SetTextSize(0.05);
    TString hv_c = Form("Operating Voltage = %d V",  high_voltage);
    text_t = ptstats->AddText(hv_c);
    TString ent = Form("Entries = %d", window_count );
    text_t = ptstats->AddText(ent);
    TString noise_t = Form("Electronics Noise Width = %.3f pC", elec_width);
    text_t = ptstats->AddText(noise_t);
    TString peak_t = Form("Charge Peak = %.3f pC", spe_mean);
    text_t = ptstats->AddText(peak_t);
    TString FWHM_t = Form("Charge FWHM = %.3f pC", FWHM);
    text_t = ptstats->AddText(FWHM_t);
    TString ptv = Form("Peak-to-valley = %.3f", Peak_to_valley);
    text_t = ptstats->AddText(ptv);
    TString HCT_t = Form("High Charge Tail = %.3f%%", HCT);
    text_t = ptstats->AddText(HCT_t);

    ptstats->SetOptStat(0);
    ptstats->SetOptFit(0);
    ptstats->SetFillColor(0);

    // This finds a good range for the plot
    int bmin_charge = charges_signal->GetXaxis()->FindBin(0.5);
    int bmax_charge = charges_signal->GetXaxis()->FindBin(3.5);
    charges_signal->GetXaxis()->SetRange(bmin_charge, bmax_charge);
    int max_bin_charge = charges_signal->GetMaximumBin();
    int max_charge_range = 1.20*charges_signal->GetBinContent(max_bin_charge);
    charges_signal->GetXaxis()->SetRange();

    charges_signal->GetListOfFunctions()->Add(ptstats);

    charges_signal->GetXaxis()->SetTitle("Charge (pC)");
    charges_signal->GetXaxis()->SetTitleFont(132);
    charges_signal->GetXaxis()->SetTitleColor(1);
    charges_signal->GetXaxis()->SetLabelFont(132);
    charges_signal->GetXaxis()->SetRangeUser(-1.0,6.0);
    charges_signal->GetYaxis()->SetRangeUser(0.0,max_charge_range);
    charges_signal->GetYaxis()->SetTitle("Events");
    charges_signal->GetYaxis()->SetTitleOffset(1.2);
    charges_signal->GetYaxis()->SetLabelFont(132);
    charges_signal->GetYaxis()->SetTitleFont(132);
    c1->Modified();

    cout << " ----- PMT Charge Parameters ----- " << endl;
    if(mode == 1){
      cout << "Electronic Noise Width is " << p22 << "pC" << endl;
      cout << "Peak " << spe_mean << ", Valley " << min_function << endl;
      cout << "The Charge FWHM is " << spe_sigma*2*sqrt(2*log(2)) << "pC" << endl;
      cout << "High Charge Tail " << high_charge_entry * 100 / low_charge_entry << "%" << endl;
    }
    cout << "The charge peak is " << spe_mean << "pC" << endl;
    cout << "The charge peak-to-valley is " << Peak_to_valley << endl;

    // get time difference
    for(size_t j = 0; j < trigger_time.size(); j++){
       Timing->Fill(pmt_time[j] - trigger_time[j]);
    }

    float time_bin_max = Timing->GetMaximumBin();
    float time_bin_max_content = Timing->GetBinContent(time_bin_max);
    TAxis *time_axis = Timing->GetXaxis();

    // Find the half-heights to integrate around
    float low_time_bin_fit = 0;
    for(int i = time_bin_max; i > time_bin_max - 10; i--){
       if(Timing->GetBinContent(i) < time_bin_max_content/10){
          low_time_bin_fit = time_axis->GetBinCenter(i);
          break;
       }
    }

    float high_time_bin_fit = 0;
    for(int i = time_bin_max; i < time_bin_max + 10; i++){
       if(Timing->GetBinContent(i) < time_bin_max_content/10){
          high_time_bin_fit = time_axis->GetBinCenter(i);
          break;
       }
    }

    // draw the fits
    TCanvas *c2 = new TCanvas("c2","Timing",200,10,700,500);
    Timing->Fit("gaus","q","",low_time_bin_fit, high_time_bin_fit);
    TF1 *fitfunc=Timing->GetFunction("gaus");
    Double_t f1 = fitfunc->GetParameter(1);
    Double_t f2 = fitfunc->GetParameter(2);
    gPad->SetLogy();
    Timing->SetStats(0);
    Timing->SetLineColor(kBlack);
    Timing->Draw();
    c2->SetLogy();
    c2->Update();

    // Draw the charge for prompt hits
    for(size_t i = 0; i < kTime.size(); i++){
       if(kTime[i] < f1 + 1 && kTime[i] > f1 - 1){
          prompt_charge->Fill(kCharge[i]);
       }
    }

    Double_t Prompt_sigma = f2;
    Double_t Prompt_FWHM = 2*f2*sqrt(2*log(2));

    double bmin = charges_signal->GetXaxis()->FindBin(charge_cut_integral);
    double bmax = charges_signal->GetXaxis()->FindBin(100);
    double coincidence_integral = charges_signal->Integral(bmin,bmax);
    double coincidence_pct = (coincidence_integral)/(window_count)*100;

    double time_bin_prompt_low = Timing->FindBin(f1-3*f2);
    double time_bin_prompt_high = Timing->FindBin(f1+3*f2);
    double coincidence_timing_integral = Timing->Integral(time_bin_prompt_low, time_bin_prompt_high);
    double coincidence_timing_pct = (coincidence_timing_integral)/(window_count)*100;

    double time_bin_late_low = Timing->FindBin(f1+10);
    double time_bin_late_high = Timing->FindBin(f1+length_late);
    double late_timing_integral = Timing->Integral(time_bin_late_low, time_bin_late_high);
    double late_pct = (late_timing_integral)/(coincidence_timing_integral)*100;

    double time_bin_dark_low = Timing->FindBin(f1-45); //45-15ns before prompt peak
    double time_bin_dark_high = Timing->FindBin(f1-15);
    double dark_timing_integral = Timing->Integral(time_bin_dark_low, time_bin_dark_high);
    double dark_rate = (dark_timing_integral)/(30*1.0e-9*window_count);

    cout << " ----- PMT Rate Parameters ----- " << endl;
    cout << "Coincidence percent of all hits " << coincidence_pct << endl;
    cout << "Coincidence perecent of prompt hits " << coincidence_timing_pct << endl;
    cout << "Late Hit Percent " << late_pct << endl;
    cout << "Dark Rate " << dark_rate << endl;

    double dark_rate_correction = (coincidence_integral - dark_timing_integral)/(window_count)*100;
    cout << dark_timing_integral << endl;
    cout << coincidence_integral << endl;
    cout << window_count << endl;
    cout << dark_rate_correction << endl;


    if(mode == 2){
      for(size_t j = 0; j < trigger_time.size(); j++){
        double dt = (pmt_time[j] - trigger_time[j]);
        double q = kCharge[j];
        if(dt > f1 - 8 && dt < f1 + 8){
          prompt_hits->Fill(dt);
          prompt_hits_charge->Fill(q);
        }
        else if(dt > f1 + 8 && dt < f1 + length_late){
          late_hits->Fill(dt);
          late_hits_charge->Fill(q);
        }
        else{
          dark_hits->Fill(dt);
          dark_hits_charge->Fill(q);
        }
      }
      for(size_t j = 0; j < prompt_pulse.size(); j++){
        // Only select pulses with correct late && prompt pulsing time
        if(late_pulse[j] > f1 + 8 && late_pulse[j] < f1 + length_late &&  
          prompt_pulse[j] > f1 - 8 && prompt_pulse[j] < f1 + 8){
          double_hits_prompt->Fill(prompt_pulse[j]);
          double_hits_prompt_charge->Fill(prompt_q[j]);
          prompt_hits->Fill(prompt_pulse[j], -1);
          prompt_hits_charge->Fill(prompt_q[j], -1);
          double_hits_late->Fill(late_pulse[j]);
          double_hits_late_charge->Fill(late_q[j]);
          late_hits_charge->Fill(late_q[j], -1);
          late_hits->Fill(late_pulse[j], -1);
        }
      }
    }

    // make the statistics legend
    TPaveStats *ttstats = new TPaveStats(0.6,0.6,0.98,0.98,"brNDC");
    ttstats->SetName("stats");
    ttstats->SetBorderSize(2);
    ttstats->SetTextAlign(12);
    ttstats->SetTextFont(132);
    ttstats->SetTextSize(0.035);
    ttstats->SetShadowColor(0);

    // statistics for timing histogram
    TText *text_tt = ttstats->AddText(box_title);
    text_tt->SetTextSize(0.05);
    TString hv_t = Form("Operating Voltage = %d V",  high_voltage);
    text_tt = ttstats->AddText(hv_t);
    TString t_ent = Form("Hits above noise = %.0f", Timing->GetEntries() );
    text_tt = ttstats->AddText(t_ent);
    TString pulse_sigma_t = Form("Prompt Sigma = %.3f ns", Prompt_sigma);
    text_tt = ttstats->AddText(pulse_sigma_t);
    TString pulse_fwhm_t = Form("Prompt FWHM = %.3f ns", Prompt_FWHM);
    text_tt = ttstats->AddText(pulse_fwhm_t);
    TString coincidence_t = Form("Prompt Coincidence Rate = %.3f%%", coincidence_timing_pct);
    text_tt = ttstats->AddText(coincidence_t);
    TString darkrate_t = Form("Dark Rate = %.0f Hits/s", dark_rate);
    text_tt = ttstats->AddText(darkrate_t);
    TString late_ratio_t = Form("Late Ratio =%.3f%%", late_pct);
    text_tt = ttstats->AddText(late_ratio_t);

    ttstats->SetOptStat(0);
    ttstats->SetOptFit(0);
    ttstats->SetFillColor(0);

    Timing->GetListOfFunctions()->Add(ttstats);
    Timing->GetXaxis()->SetTitle("#Deltat (ns)");
    Timing->GetXaxis()->SetTitleFont(132);
    Timing->GetXaxis()->SetTitleColor(1);
    Timing->GetXaxis()->SetLabelFont(132);
    Timing->GetYaxis()->SetTitle("Events");
    Timing->GetYaxis()->SetLabelFont(132);
    Timing->GetYaxis()->SetTitleFont(132);
    Timing->GetXaxis()->SetRangeUser(40.0, 105.0);
    c2->Modified();

    cout << " ----- PMT Rate Parameters ----- " << endl;
    cout << "Prompt Sigma " << f2 << endl;
    if(mode == 1){
      cout << "Prompt Mean " << f1 << endl;
    }

    // draw the fits
    TCanvas *c3 = new TCanvas("c3","c3",800,600);
    Timing->Fit("gaus","q","",low_time_bin_fit, high_time_bin_fit);
    gPad->SetLogy();
    prompt_hits->SetStats(0);
    prompt_hits->SetLineColor(kBlack);
    prompt_hits->Draw();
    late_hits->SetLineColor(kRed);
    late_hits->SetLineStyle(2);
    late_hits->Draw("same");
    double_hits_prompt->SetLineColor(kBlue);
    double_hits_prompt->SetLineStyle(3);
    double_hits_prompt->Draw("same");
    double_hits_late->SetLineColor(kCyan);
    double_hits_late->SetLineStyle(5);
    double_hits_late->SetMarkerStyle(5);
    double_hits_late->Draw("same");
    dark_hits->SetLineColor(kGreen);
    dark_hits->SetLineStyle(4);
    dark_hits->Draw("same");
    TLegend *t1 = new TLegend(0.6, 0.6, 0.8, 0.8);
    t1->SetBorderSize(0);;
    t1->SetTextFont(132);
    t1->SetFillColor(0);
    t1->AddEntry(prompt_hits, "Prompt");
    t1->AddEntry(late_hits, "Late");
    t1->AddEntry(double_hits_prompt, "Double (1)");
    t1->AddEntry(double_hits_late, "Double (2)");
    t1->AddEntry(dark_hits, "Dark");
    t1->Draw("same");
    c3->SetLogy();
    c3->Update();

    // Output Histograms to File
    if (argc > 2){
      TFile f(argv[2],"new");
      cout << "Writing data to: " << argv[2] << endl;

      pedestals->Write();
      variances->Write();
      charges_signal->Write();
      peaks->Write();
      average_waveform->Write();

      trigger_pedestals->Write();
      average_waveform_trigger->Write();
      Timing->Write();
      if(mode == 2){
        prompt_hits->Write();
        late_hits->Write();
        dark_hits->Write();
        double_hits_prompt->Write();
        double_hits_late->Write();
        prompt_hits_charge->Write();
        late_hits_charge->Write();
        dark_hits_charge->Write();
        double_hits_prompt_charge->Write();
        double_hits_late_charge->Write();
        prompt_charge->Write();
      }
      if(mode == 3){
        parameter0->Write();
        parameter1->Write();
        parameter2->Write();
        parameter3->Write();
        parameter4->Write();
        parameter5->Write();
        parameter6->Write();
        parameter7->Write();
        parameter8->Write();
      }
    }

    // clean-up
    delete pedestals;
    delete variances;
    delete charges_signal;
    delete peaks;
    delete average_waveform;

    delete parameter0;
    delete parameter1;
    delete parameter2;
    delete parameter3;
    delete parameter4;
    delete parameter5;
    delete parameter6;
    delete parameter7;
    delete parameter8;

    delete trigger_pedestals;
    delete average_waveform_trigger;
    delete Timing;
    delete prompt_hits;
    delete late_hits;
    delete dark_hits;
    delete double_hits_prompt;
    delete double_hits_late;
    delete prompt_hits_charge;
    delete late_hits_charge;
    delete dark_hits_charge;
    delete double_hits_prompt_charge;
    delete double_hits_late_charge;
    delete prompt_charge;
    delete fit_waveform;

    return 0;
}

double GetSampling(char* datafile){

    ifstream ifs(datafile , ifstream::in); // Open File List Stream
    ifs >> filename;

    while (ifs.good()){ // While Stream Is Open, Analyze New Files.

      file.openFile(filename, H5F_ACC_RDONLY); //Open HDF5 File
      ifs >> filename;

      // Open dataset for characterization PMT
      dataset = file.openDataSet(pmt_channel);
      horiz_interval = dataset.openAttribute("horiz_interval");
      horiz_interval.read(PredType::NATIVE_DOUBLE, &dx);
      break;
    }

    return dx;
} 

size_t GetLength(char* datafile){

    ifstream ifs(datafile , ifstream::in); // Open File List Stream
    ifs >> filename;
    size_t window_length = 0;
    while (ifs.good()){ // While Stream Is Open, Analyze New Files.

      file.openFile(filename, H5F_ACC_RDONLY); //Open HDF5 File
      ifs >> filename;

      dataset = file.openDataSet(pmt_channel);

      DataCluster * datacluster = Init_Data(&dataset);
      window_length = datacluster->trace_length;
      break;
    }
    return window_length;
}

int GetWindowStart(char* datafile){

    ifstream ifs(datafile , ifstream::in); // Open File List Stream
    ifs >> filename;

    int window_count = 0;

    while (ifs.good()){ // While Stream Is Open, Analyze New Files.

      file.openFile(filename, H5F_ACC_RDONLY); //Open HDF5 File
      ifs >> filename;

      // Open dataset for characterization PMT
      dataset = file.openDataSet(pmt_channel);
      horiz_interval = dataset.openAttribute("horiz_interval");
      vertical_gain = dataset.openAttribute("vertical_gain");
      horiz_interval.read(PredType::NATIVE_DOUBLE, &dx);
      vertical_gain.read(PredType::NATIVE_DOUBLE, &dy);

      // Open dataset for trigger PMT. Sampling time is same for both channels.
      dataset_trigger = file.openDataSet(trigger_channel);
      vertical_gain_trigger = dataset_trigger.openAttribute("vertical_gain");
      vertical_gain_trigger.read(PredType::NATIVE_DOUBLE, &dy2);

      DataCluster * datacluster = Init_Data(&dataset);
      DataCluster * datacluster_trigger = Init_Data(&dataset_trigger);

      unsigned long window_length = datacluster->trace_length;

      cout << "Finding peak of average waveforms." << endl;

      for(size_t j = 0; j < datacluster->n_traces; j++){

        Read_Trace(datacluster,j);
        Read_Trace(datacluster_trigger,j);
        window_count++;

        for(size_t i = 0; i < window_length; i++){
          float voltage = ((float)datacluster->data_out[i]*dy);
          waveform_voltage[i] += voltage;
          float voltage_trigger = ((float)datacluster_trigger->data_out[i]*dy2);
          waveform_voltage_trigger[i] += voltage_trigger;
        }
      }
      file.close();
      break; // one file
    }
    ifs.close();

    return window_count;
}
