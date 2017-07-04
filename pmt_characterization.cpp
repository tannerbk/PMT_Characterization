#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <numeric>

// Root Libraries
#include <TFile.h>
#include <TMath.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>
#include <TROOT.h>
#include <TPaveStats.h>
#include <TStyle.h>
#include "TLatex.h"
// HDF5 Library
#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

using namespace std;

const int RANK_OUT = 2;

/* Scope Settings */
const int length_trace = 10002;
const unsigned long number_entries = 100000;
const unsigned long number_files = 20;
const float sampling_time = 0.05;
const string trigger_channel = "channel2";
const string pmt_channel = "channel3";

/* Plot and Output Settings */
const int debug = 1; // Set to 0 for debugging output
const int high_voltage = 1840;
TString box_title = Form("R5912-MOD");

/* Analysis settings */
const int length_late = 60;
const float variance_cut = 0.005;
// Charge cuts to fill timing histogram
const float charge_cut_low = 1.0;
const float charge_cut_high = 3.0;
// Charge cuts to integrate above to count coincidence
const float charge_cut_integral = 0.3;
const int size_signal_window = 400;
const int number_of_windows = (length_trace)/(size_signal_window);

static double time_max_trig[number_files][number_entries] = {{ 0 }};
static double time_max_trig2[number_files][number_entries][number_of_windows] = {{{ 0 }}};
static double delta_T[number_files][number_entries][number_of_windows] = {{{ 0 }}};

struct DataCluster{
  DataSet *dataset; // Dataset pointer
  DataSpace dataspace; // DataSet's DataSpace
  DataSpace memspace; // MemSpace Object for Data Extraction
  hsize_t offset[RANK_OUT]; // Data Extraction Parameters...
  hsize_t count[RANK_OUT];
  hsize_t offset_out[RANK_OUT];
  hsize_t count_out[RANK_OUT];
  unsigned long trace_length; // Length of a Scope Trace
  unsigned long n_traces; // Number of traces in DataSet
  char * data_out; // Pointer to Data Buffer
};

typedef struct DataCluster DataCluster;

/* DataCluster Methods */
DataCluster * Init_Data(DataSet *dataset);
int Read_Trace(DataCluster *datacluster, unsigned long trace_index);

int main (int argc, char* argv[])
{

  unsigned long window_width = atoi(argv[3]);
  const int termination_ohms = 50;

  try{
      cout << "Characterization PMT: " << box_title << endl;

      // Attribute Variables
      Attribute horiz_interval;
      Attribute vertical_gain;
      Attribute vertical_gain_trigger;

      // Attribute Value Variables
      double dx,dy,dy2;

      TH1F *pedestals = new TH1F("Pedestal","",100000,-1,1);
      TH1F *variances = new TH1F("Variance","",500000,-0.1,1.0);
      TH1F *charges_signal = new TH1F("Charge","",800,-4.0,12.0);
      TH1F *average_waveform = new TH1F("Average_Waveform","", length_trace, 0, length_trace*sampling_time);
      TH1F *trigger_pedestals = new TH1F("Pedestal_Trigger","",100000,-1.0,1.0);
      TH1F *average_waveform_trigger = new TH1F("Average_Waveform_Trigger","", length_trace, 0, length_trace*sampling_time);
      TH1F *Timing = new TH1F("Timing","",(length_trace/2)*sampling_time*2, 0, (length_trace/2)*sampling_time);

      int window = 0; // window count
      int r = 0; // file count
      int dark_count_window = 0;

      size_t i;
      size_t k;
      size_t window_count = 0.0;
      string filename;
      H5File file;
      DataSet dataset;
      DataSet dataset_trigger;

      // average waveform
      float waveform_voltage[length_trace] = {0};
      float waveform_voltage_trigger[length_trace] = {0};

      ifstream ifs ( argv[1] , ifstream::in ); // Open File List Stream
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

        unsigned long window_length = datacluster->trace_length;

        if(float(dx*1e9) != float(sampling_time)){
           cout << "WARNING: Sampling time mismatch: " << dx*1e9 << " " << sampling_time << endl;
        }
        if(window_length != length_trace){
           cout << "WARNING: Window length mismatch: " << window_length << " " << length_trace << endl;
        }

        cout << "File " << r  << ", finding peak of average waveform." << endl;

        for(size_t j = 0; j < datacluster->n_traces; j++){

          Read_Trace(datacluster,j);
          Read_Trace(datacluster_trigger,j);
          float pedestal = TMath::Mean (window_width, datacluster->data_out)*dy;
          float trigger_pedestal = TMath::Mean (window_width, datacluster_trigger->data_out)*dy2;
          window_count++;

          for(i=0; i < window_length; i++){
            float voltage = ((float)datacluster->data_out[i]*dy-pedestal);
            waveform_voltage[i] += voltage;
            float voltage_trigger = ((float)datacluster_trigger->data_out[i]*dy2-trigger_pedestal);
            waveform_voltage_trigger[i] += voltage_trigger;
          }
        }
        file.close();
        break; // one file
      }
      ifs.close();

      Float_t avg_waveform_voltage[length_trace];
      for(i = 0; i < length_trace; i++){
	avg_waveform_voltage[i] = waveform_voltage[i]/window_count;
	average_waveform->SetBinContent(i+1, avg_waveform_voltage[i]);
      }

      Float_t avg_waveform_voltage_trigger[length_trace];
      for(i = 0; i < length_trace; i++){
	avg_waveform_voltage_trigger[i] = waveform_voltage_trigger[i]/window_count;
	average_waveform_trigger->SetBinContent(i+1, avg_waveform_voltage_trigger[i]);
      }

      window_count = 0;
      unsigned long pre_window = 12 / (dx*1e9);
      unsigned long pre_trigger_window = 20 / (dx*1e9);
      unsigned long window_start = average_waveform->GetMinimumBin() - pre_window;
      unsigned long window_start_trigger = average_waveform_trigger->GetMinimumBin() - pre_trigger_window;
      cout << " ----- Window Start Times ----- " << endl;
      cout << "Integration start time " << window_start*dx*1e9 << " ns." << endl;
      cout << "Trigger start time " << window_start_trigger*dx*1e9 << " ns." << endl;
      const int number_windows = (length_trace)/(2*size_signal_window);

      ifstream ifs2 ( argv[1] , ifstream::in ); // Open File List Stream
      ifs2 >> filename;

      while (ifs2.good()){ // While Stream Is Open, Analyze New Files.

        file.openFile(filename, H5F_ACC_RDONLY); //Open HDF5 File
        ifs2 >> filename;

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
        if(r == 0){
           cout << " ----- PMT Characterization Parameters ----- " << endl;
           cout << "File " << r << ": " << filename << endl;
           cout << "The time bin width is " << dx*1e9 << " ns" << endl;
           cout << "The trace length is " << window_length*dx*1e9 << " ns" << endl;
           cout << "The pedestal window is 0 - " << window_width*dx*1e9 << " ns" << endl;
           cout << "The trigger window is " << window_start_trigger*dx*1e9 << " - " << window_start_trigger*dx*1e9 + 40 << " ns" << endl;
           cout << "The pmt window is " << window_start*dx*1e9 << " - " << window_start*dx*1e9 + 30 << " ns" << endl;
           cout << "The vertical resolution is " << dy*1000 << " mV." << endl;
           cout << "The vertical resoution of the trigger is " << dy2*1000 << " mV." << endl;
           cout << "Analyzing Files: " << endl;
        }
        else{
           cout << "File " << r ": " << filename << endl;
        }

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
          double dcharge = 0.0;
          variance = 0.0;
          k = 0;

          for(i=0; i < window_width; i++){
            float voltage = ((float)datacluster->data_out[i]*dy-pedestal);
            variance += voltage * voltage;
          }

          // Only cut on variance of test PMT
          // traces that don't pass this cut might have dark pulses in the ped window
          if(variance < variance_cut) {

            for(i = window_start; i < window_start + size_signal_window; i++){
      	      float voltage = ((float)datacluster->data_out[i]*dy-pedestal);
      	      ncharge += (voltage*((-1000.0*dx*1e9)/termination_ohms));
            }

            for(i = 0; i < size_signal_window; i++){
      	      float voltage = ((float)datacluster->data_out[i]*dy-pedestal);
      	      dcharge += (voltage*((-1000.0*dx*1e9)/termination_ohms));
            }

            if(dcharge > charge_cut_integral){
              dark_count_window++;
            }

            for(window = 0; window < number_windows; window++){

              // need to initialze these all for every window
              double window_charge = 0.0;
              max_time_entry = 0.0;
              time_bin_counter = 0;
              time_bin_counter2 = 0;

      	      for(i = window_counter + window_width; i < window_counter + window_width + size_signal_window; i++){
      	        float voltage = ((float)datacluster->data_out[i]*dy-pedestal);
      	        window_charge += (voltage*((-1000.0*dx*1e9)/termination_ohms));
      	      }

      	      if(window_charge > charge_cut_low && window_charge < charge_cut_high){
      	        for(i = window_counter + window_width; i < window_counter + window_width + size_signal_window; i++){
      	          float voltage = ((float)datacluster->data_out[i]*dy-pedestal);
      	          if(voltage < max_time_entry){
      	            max_time_entry = voltage;
      	            time_bin_counter = i;
      	          }
      	        }
      	      }

      	      prompt_peak = ((float)datacluster->data_out[time_bin_counter]*dy-pedestal);

      	      if(window_charge > charge_cut_low && window_charge < charge_cut_high){
                for(i = time_bin_counter - 20/(dx*1e9); i < time_bin_counter; i++){
      	          float voltage = ((float)datacluster->data_out[i]*dy-pedestal);
      	          if(voltage > prompt_peak*0.2){
      	            time_bin_counter2 = i;
      	          }
      	        }
      	      }

      	      // convert bin to time
      	      time_max_trig2[r][j][k] =  time_bin_counter2*dx*1e9;
      	      k++;
      	      window_counter = window_counter + size_signal_window;
            }

            max_time_entry = 0;
            time_bin_counter_trigger = 0;
            time_bin_counter2_trigger = 0;

            for(i = window_start_trigger; i < window_start_trigger + 40/(dx*1e9); i++){
              float voltage = ((float)datacluster_trigger->data_out[i]*dy2-trigger_pedestal);
              if(voltage < max_time_entry){
      	        max_time_entry = voltage;
      	        time_bin_counter_trigger = i;
              }
            }

            prompt_peak = ((float)datacluster_trigger->data_out[time_bin_counter_trigger]*dy2-trigger_pedestal);

            // Look 20ns before peak for 20% crossing
            for(i = time_bin_counter_trigger - 20/(dx*1e9); i < time_bin_counter_trigger; i++){
              float voltage = ((float)datacluster_trigger->data_out[i]*dy2-trigger_pedestal);
              if(voltage > prompt_peak*0.2){
      	        time_bin_counter2_trigger = i;
              }
            }

            // convert bin to time
            time_max_trig[r][j] = time_bin_counter2_trigger*dx*1e9;

            window_count++; // number of traces
            charges_signal->Fill(ncharge);
            pedestals->Fill(pedestal);
            trigger_pedestals->Fill(trigger_pedestal);
          }
          //store histograms
          variances->Fill(variance);
        }
        r++;
        file.close();
      }
      ifs2.close();

      if(debug == 0){
        cout  << "The number of Traces accepted is " << window_count << endl;
        cout  << "The number of files is " << r  << endl;
      }

      // If this code wasn't bad enough, its about to get worse
      // continue reading at your own peril

      TAxis *axis = charges_signal->GetXaxis();

      // spe fit
      Int_t begin_spe_peak = 1.0;
      Int_t end_spe_peak = 10.0;
      Int_t first_bin = axis->FindBin(begin_spe_peak);
      Int_t last_bin = axis->FindBin(end_spe_peak);
      Double_t charge_finder[last_bin - first_bin];
      Int_t f = 0.0;
      Int_t b = 0.0;
      Double_t max_charge = 0.0;
      Double_t max_charge_entry = 0.0;
      Double_t first_bin2 = first_bin;

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

      for(i = low_charge_bin; i < size_t(charges_signal->GetNbinsX()); i++){
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
      TString ent = Form("Entries = %zu", window_count );
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
      if(debug == 0){
        cout << "Electronic Noise Width is " << p22 << "pC" << endl;
        cout << "Peak " << spe_mean << ", Valley " << min_function << endl;
        cout << "The Charge FWHM is " << spe_sigma*2*sqrt(2*log(2)) << "pC" << endl;
        cout << "High Charge Tail " << high_charge_entry * 100 / low_charge_entry << "%" << endl;
      }
      cout << "The charge peak is " << spe_mean << "pC" << endl;
      cout << "The charge peak-to-valley is " << Peak_to_valley << endl;

      // get time difference
      for(size_t j = 0; j < number_files; j++){
	for(i=0; i < number_entries; i++){
	  for(k=0; k < number_of_windows; k++){
            if(time_max_trig2[j][i][k] != 0){
              delta_T[j][i][k] = time_max_trig2[j][i][k] - time_max_trig[j][i];
	      Timing->Fill(delta_T[j][i][k]);
	    }
	  }
	}
      }

      float time_bin_max = Timing->GetMaximumBin();
      float time_bin_max_content = Timing->GetBinContent(time_bin_max);
      TAxis *time_axis = Timing->GetXaxis();

      // Find the half-heights to integrate around
      float low_time_bin_fit = 0;
      for(i = time_bin_max; i > time_bin_max - 10; i--){
         if(Timing->GetBinContent(i) < time_bin_max_content/10){
            low_time_bin_fit = time_axis->GetBinCenter(i);
            break;
         }
      }

      float high_time_bin_fit = 0;
      for(i = time_bin_max; i < time_bin_max + 10; i++){
         if(Timing->GetBinContent(i) < time_bin_max_content/10){
            high_time_bin_fit = time_axis->GetBinCenter(i);
            break;
         }
      }

      // draw the fits
      TCanvas *c2 = new TCanvas("c2","Timing",200,10,700,500);
      Timing->Fit("gaus","q","",low_time_bin_fit, high_time_bin_fit);
      TF1 *fitfunc=Timing->GetFunction("gaus");
      Double_t f1=fitfunc->GetParameter(1);
      Double_t f2=fitfunc->GetParameter(2);
      gPad->SetLogy();
      Timing->Draw();
      c2->SetLogy();
      c2->Update();

      // get rates
      Double_t coincidence_count = 0.0;
      Double_t dark_count = 0.0;
      Double_t late_pulse_count = 0.0;
      Double_t pre_pulse_count = 0.0;
      Double_t post_pulse_count = 0.0;

      /* I do the coincidence and dark count two different ways. This 
         way is less robust to the details of the analysis. The other way
         just uses integrals of the charge distributions to count, so 
         it should be better for the absolute numbers */
      for(size_t j = 0; j < number_files; j++){
	for(i = 0; i < number_entries; i++){
	  for(k = 0; k < number_of_windows; k++){
	    if(time_max_trig2[j][i][k] != 0){
              // if delta t is within 3 sigma of prompt peak
	      if(delta_T[j][i][k] < f1 + 3*f2 && delta_T[j][i][k] > f1 - 3*f2){
		coincidence_count++;
	      }
              // 5 sigma after prompt peak to length late
	      else if(delta_T[j][i][k] > f1 + 5*f2 && delta_T[j][i][k] < f1 + length_late){
		late_pulse_count++;
	      }
	      else if(delta_T[j][i][k] < f1 + 5*f2 && delta_T[j][i][k] > f1 + 3*f2){
		post_pulse_count++;
	      }
	      else if(delta_T[j][i][k] < f1 - 3*f2 && delta_T[j][i][j] > f1 - 10){
		pre_pulse_count++;
	      }
	      else{
		dark_count++;
	      }
	    }
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

      Double_t delta_t_dark = (f1 - 10) + ((length_trace/2)*0.1 - f1 - length_late);
      Double_t delta_t_coincidence = 6*f2;
      Double_t delta_t_late = length_late - 5*f2;
      Double_t delta_t_pre = 10 - 3*f2;
      Double_t t_not_used = 2*f2;

      Double_t ns = 1.0e-9;
      Double_t Prompt_sigma = f2;
      Double_t Prompt_FWHM = 2*f2*sqrt(2*log(2));
      Double_t Dark_Rate = dark_count/(delta_t_dark*ns*window_count);
      Double_t dark_pulse_count1 = Dark_Rate*delta_t_coincidence*ns*window_count;
      Double_t dark_pulse_count2 = Dark_Rate*delta_t_late*ns*window_count;
      Double_t Coincidence_Percent = (coincidence_count - dark_pulse_count1)/(window_count)*100;
      Double_t Late_Hit_Percent = (late_pulse_count - dark_pulse_count2)/(Timing->GetEntries())*100;

      if(debug == 0){
        cout << "Coincidence Hits " << coincidence_count << endl;
        cout << "Late Pulse Hits " << late_pulse_count << endl;
        cout << "Dark Pulse Hits " << dark_count << endl;
        cout << "Pre Pulse Hits " << pre_pulse_count << endl;
        cout << "Post Pulse Hits " << post_pulse_count << endl;
        cout << "Total hits = " << coincidence_count + late_pulse_count + dark_count + pre_pulse_count + post_pulse_count << endl;
        cout << "T coinc " << delta_t_coincidence  << endl;
        cout << "T late " << delta_t_late  << endl;
        cout << "T dark " << delta_t_dark << endl;
        cout << "Total time = " << delta_t_dark + delta_t_coincidence + delta_t_late + delta_t_pre + t_not_used << endl;
      }

      // A Better way to do the coincidence count is integrate the charge distribution above
      // some value
      double coincidence_pct = 0;
      double bmin = charges_signal->GetXaxis()->FindBin(charge_cut_integral);
      double bmax = charges_signal->GetXaxis()->FindBin(100);
      double coincidence_integral = charges_signal->Integral(bmin,bmax); 
      coincidence_pct = (coincidence_integral - dark_count_window)/(window_count)*100;

      cout << " ----- PMT Rate Parameters ----- " << endl;      
      cout << "Coincidence pulses " << coincidence_integral << " Dark pulses " << dark_count_window << endl;
      cout << "Coincidence percent based on charge integral " << coincidence_pct << endl;

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
      TString coincidence_t = Form("Prompt Coincidence Rate = %.3f%%", Coincidence_Percent);
      text_tt = ttstats->AddText(coincidence_t);
      TString darkrate_t = Form("Dark Rate = %.0f Hits/s", Dark_Rate);
      text_tt = ttstats->AddText(darkrate_t);
      TString late_ratio_t = Form("Late Ratio =%.3f%%", Late_Hit_Percent);
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
      c2->Modified();

      cout << " ----- PMT Rate Parameters ----- " << endl;      
      cout << "Prompt Mean " << f1 << endl;
      cout << "Prompt Sigma " << f2 << endl;
      if(debug == 0){
        cout << "Coincidence Percent " << Coincidence_Percent << endl;
        cout << "Late Hit Percent " << Late_Hit_Percent << endl;
        cout << "Dark rate correction to coincidence number " <<  dark_pulse_count1 << ", number of coincidence events " << coincidence_count << endl;
        cout << "Dark rate correction to late number " <<  dark_pulse_count2 << ", number of late events " << late_pulse_count << endl;
      }

      // Output Histograms to File
      if (argc > 2){
	TFile f(argv[2],"new");

	pedestals->Write();
	variances->Write();
	charges_signal->Write();
	average_waveform->Write();

	trigger_pedestals->Write();
	average_waveform_trigger->Write();
	Timing->Write();
      }

      // clean-up
      delete pedestals;
      delete variances;
      delete charges_signal;
      delete average_waveform;

      delete trigger_pedestals;
      delete average_waveform_trigger;
      delete Timing;
    } // end of try block

      // catch failure caused by the H5File operations
  catch( FileIException error )
    {
      error.printError();
      return -1;
    }

  // catch failure caused by the DataSet operations
  catch( DataSetIException error )
    {
      error.printError();
      return -1;
    }

  // catch failure caused by the DataSpace operations
  catch( DataSpaceIException error )
    {
      error.printError();
      return -1;
    }

  // catch failure caused by the DataSpace operations
  catch( DataTypeIException error )
    {
      error.printError();
      return -1;
    }

  return 0; // successfully terminated
}

DataCluster * Init_Data(DataSet *dataset){

  DataCluster * datacluster = new DataCluster[1];

  /*
   * Get dataspace of the dataset.
   */
  datacluster->dataset = dataset;
  datacluster->dataspace = datacluster->dataset->getSpace();

  /*
   * Get the dimension size of each dimension in the dataspace and
   * display them.
   */

  hsize_t dims_out[2];
  datacluster->dataspace.getSimpleExtentDims( dims_out, NULL);
  datacluster->trace_length = (unsigned long)(dims_out[1]);
  datacluster->n_traces = (unsigned long)(dims_out[0]);

  // cout << "Reading " << datacluster->n_traces << " traces of length " << datacluster->trace_length << "..." << endl;

  // Data Buffer
  datacluster->data_out = new char[datacluster->trace_length]; // Scope data is size char
  for (unsigned long i = 0; i < datacluster->trace_length; i++) datacluster->data_out[i]= 0;

  /*
   * Define hyperslab in the dataset.
   */

  datacluster->offset[0] = 0;
  datacluster->offset[1] = 0;
  datacluster->count[0] = 1;
  datacluster->count[1] = datacluster->trace_length;
  datacluster->dataspace.selectHyperslab( H5S_SELECT_SET, datacluster->count, datacluster->offset );

  /*
   * Define the memory dataspace.
   */

  hsize_t dimsm[2]; /* memory space dimensions */
  dimsm[0] = dims_out[0];
  dimsm[1] = dims_out[1];
  datacluster->memspace = DataSpace( RANK_OUT, dimsm );

  /*
   * Define memory hyperslab.
   */

  datacluster->offset_out[0] = 0;
  datacluster->offset_out[1] = 0;
  datacluster->count_out[0] = 1;
  datacluster->count_out[1] = datacluster->trace_length;
  datacluster->memspace.selectHyperslab( H5S_SELECT_SET, datacluster->count_out, datacluster->offset_out );

  return datacluster;
}

/* Method: Read_Trace(DataCluster *datacluster, unsigned long trace_index)
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Updates a DataCluster datacluster so that its buffer contains trace number trace_index
 */

int Read_Trace(DataCluster *datacluster, unsigned long trace_index){
  datacluster->offset[0]= (hsize_t)trace_index;
  datacluster->dataspace.selectHyperslab( H5S_SELECT_SET, datacluster->count, datacluster->offset );
  datacluster->memspace.selectHyperslab( H5S_SELECT_SET, datacluster->count_out, datacluster->offset_out );
  datacluster->dataset->read( datacluster->data_out, PredType::NATIVE_CHAR, datacluster->memspace, datacluster->dataspace );
  return 0; // No protection...
}

