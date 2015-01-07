#include <iostream>
#include <string>
#include <fstream>
#include <math.h>

// Root Libraries
#include <TFile.h>
#include <TMath.h>
#include <TRandom2.h>
#include <TMathBase.h>
#include <TTree.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TROOT.h>
#include <TPaveStats.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TLine.h>
#include <TGaxis.h>

// HDF5 Library
#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

using namespace std;

const int RANK_OUT = 2; // Data Rank

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

  try
    {

      // Attribute Variables
      Attribute horiz_interval;
      Attribute vertical_gain;
      Attribute horiz_offset;
      Attribute vertical_offset;
      Attribute max_value;
      Attribute min_value;

      // Attribute Value Variables
      double dx,dy,xoffset,yoffset;
      double Vmin, Vmax;

      //Define ROOT histograms here 
      TH1F *voltages = new TH1F("Voltages over Pedestal Window for ETL PMT","",5000,-1.0,1.0);
      TH1F *noise_voltages = new TH1F("Voltages over Signal Window for ETL PMT","",5000,-1.0,1.0);
      TH1F *variances = new TH1F("Variance over Pedestal Window for ETL PMT","",500000,-0.1,1.0);
      TH1F *noise_variances = new TH1F("Variance over Signal Window for ETL PMT","",500000,-0.1,1.0);
      TH1F *noise_pedestals = new TH1F("Pedestal for ETL PMT","",100000,-1,1); 
      TH1F *corrected_noise_pedestals = new TH1F("Corrected Pedestal for ETL PMT","",100000,-1,1);
      TH1F *charges_noise = new TH1F("Charge ETL PMT","",20000,-150,200);   
      TH1F *acceptable_waveform = new TH1F("Acceptable Waveform ETL PMT","", 2502, 0, 502); 
      TH1F *pulse_waveform = new TH1F("Unacceptable Waveform ETL PMT","", 2502, 0, 502); 
      TH1F *average_waveform = new TH1F("Average Waveform ETL PMT","", 2502, 0, 502);  
      
      TH1F *big_noise_pedestals = new TH1F("Pedestal for Trigger PMT","",5000,-1.0,1.0); 
      TH1F *big_noise_voltages =  new TH1F("Voltages over Signal Window for Trigger PMT","",5000,-1.0,1.0);
      TH1F *big_noise_charges =  new TH1F("Signal Charge for Trigger PMT","",5000,-150,200); 
      TH1F *big_average_waveform = new TH1F("Average Waveform for Trigger PMT","", 2502, 0, 502); 

      TH1F *Timing = new TH1F("Timing histogram","",300, 0, 250); 

      unsigned long i; 
      unsigned long j; 
      unsigned int window_count = 0.0;
      unsigned int new_window_count = 0.0;  
      string filename; 
      H5File file; 
      DataSet dataset;

      // average waveform
      const unsigned long length_trace = 2502;  // hard coded in, shitty 
      const unsigned long number_entries = 200000; // hard coded in, shitty 
      float waveform_voltage[length_trace] = {0};
      float big_waveform_voltage[length_trace] = {0};
      float timing_voltage[length_trace] = {0};
      float timing_voltage2[length_trace] = {0};
      double time_max_trig[number_entries] = {0}; 
      double time_max_trig2[number_entries] = {0}; 
      double delta_T[number_entries] = {0}; 

      // seperates file input 
      int g = 0; 
      int p = 0;
      int a = 0; 

      ofstream chargefile;
      chargefile.open("info.text"); // used mostly for debugging 
	  
      for(g=0; g<1; g++){

	ifstream ifs ( argv[1] , ifstream::in ); // Open File List Stream
	ifs >> filename;

	while (ifs.good()){ // While Stream Is Open, Analyze New Files.

	  cout<<filename<<endl;
	  file.openFile(filename, H5F_ACC_RDONLY); //Open HDF5 File 
	  ifs >> filename;
	  dataset = file.openDataSet("channel3"); // Open HDF5 DataSet

	  // Reading in Attributes
	  horiz_interval = dataset.openAttribute("horiz_interval");
	  vertical_gain = dataset.openAttribute("vertical_gain");
	  horiz_offset = dataset.openAttribute("horiz_offset");
	  vertical_offset = dataset.openAttribute("vertical_offset");
	  max_value = dataset.openAttribute("max_value");
	  min_value = dataset.openAttribute("min_value");

	  horiz_interval.read(PredType::NATIVE_DOUBLE, &dx);
	  vertical_gain.read(PredType::NATIVE_DOUBLE, &dy);
	  horiz_offset.read(PredType::NATIVE_DOUBLE, &xoffset);
	  vertical_offset.read(PredType::NATIVE_DOUBLE, &yoffset);
	  max_value.read(PredType::NATIVE_DOUBLE, &Vmax);
	  min_value.read(PredType::NATIVE_DOUBLE, &Vmin);
	  
	  DataCluster * datacluster = Init_Data(&dataset); //initialize datacluser
	
	  unsigned long window_length = datacluster->trace_length; 
	  unsigned long integrate_length = window_length - window_width; 
	  cout << " the dx is " << dx << ", the integrate length is " << integrate_length << ", the trace length is " << window_length << ", the dy is " << dy << endl;
		  
	  // for PMT1
	  double ncharge = 0.0; 
	  float voltage = 0.0; 
	  float noise_voltage = 0.0; 
	  float window_voltage = 0.0;
	  float variance = 0.0; 
	  float noise_variance = 0.0;  
	  float noise_pedestal = 0.0; 
	  float corrected_noise_pedestal = 0.0; 

	  for (j = 0; j < datacluster->n_traces; j++){
	    cout << "Analyzing trace " << j+1 << " out of " 
		 <<datacluster->n_traces<< " total traces " << "for channel 3 ETL"<<'\r';
	    cout.flush(); 

	    Read_Trace(datacluster,j); 	    	    
	    noise_pedestal = TMath::Mean (window_width, datacluster->data_out)*dy;
	  
	    // initialize 
	    ncharge = 0.0;  
	    variance = 0.0; 
	    noise_variance = 0.0; 
	    noise_voltage = 0.0; 	   
	    
	    for(i=0; i < window_width; i++){
	      voltage = ((float)datacluster->data_out[i]*dy-noise_pedestal);
	      voltages->Fill(voltage);
	      variance = variance + voltage * voltage; 
	    }
	    
	    if(variance < 0.00056) { // traces that don't pass this cut have dark pulses! 
	      corrected_noise_pedestal = TMath::Mean (window_width, datacluster->data_out)*dy; 
	      for(i = window_width; i < window_width + integrate_length ; i++){ // signal window 
		noise_voltage = ((float)datacluster->data_out[i]*dy-corrected_noise_pedestal);
		noise_voltages->Fill(noise_voltage); 
		noise_variance = noise_variance + noise_voltage * noise_voltage; 
		ncharge = ncharge+(noise_voltage*((-1000.0*dx*1e9)/termination_ohms)); // charge 
	      }  
	      
	      for(i=0; i < window_length; i++){ // entire window 
		window_voltage = ((float)datacluster->data_out[i]*dy-corrected_noise_pedestal);  
		waveform_voltage[i] = waveform_voltage[i] + window_voltage; 
	      }	   
	      
	      for(i=0; i < window_length; i++){ // draws an acceptable trace 
		acceptable_waveform->SetBinContent(i+1,(float)datacluster->data_out[i]*dy-corrected_noise_pedestal); 
	      }

	      window_count++; // number of traces 
	  
	      charges_noise->Fill(ncharge); 
	      noise_variances->Fill(noise_variance); 
	      corrected_noise_pedestals->Fill(corrected_noise_pedestal); 
	    }

	    if(variance > 0.00056){ // draw a waveform that was cut, should have a dark pulse 
	      for(i=0; i<window_length; i++){
		pulse_waveform->SetBinContent(i+1,(float)datacluster->data_out[i]*dy-noise_pedestal);
	      }
	    }
	    
	    //store histograms 
	    variances->Fill(variance); 
	    noise_pedestals->Fill(noise_pedestal);
	  }
	  
	  file.close(); 
	}
	ifs.close(); 
      }

 
      Float_t avg_waveform_voltage[length_trace];  
      for(i=0; i<length_trace; i++){ 
	avg_waveform_voltage[i] = waveform_voltage[i]/window_count;
	average_waveform->SetBinContent(i+1, avg_waveform_voltage[i]);  
      }

      cout  << "\n The number of Traces accepted is " << window_count << " for the ETL PMT." << endl;
      
      for(p=0;p<1;p++){

	ifstream ifs ( argv[1] , ifstream::in ); // Open File List Stream
	ifs >> filename;
	
	while (ifs.good()){ // While Stream Is Open, Analyze New Files.

	  cout<<filename<< " trig pmt" << endl;
	  file.openFile(filename, H5F_ACC_RDONLY); //Open HDF5 File 
	  ifs >> filename;

	  dataset = file.openDataSet("channel1"); // open channel 1
	  // re-set for new channel 
	  horiz_interval = dataset.openAttribute("horiz_interval");
	  vertical_gain = dataset.openAttribute("vertical_gain");
	  horiz_offset = dataset.openAttribute("horiz_offset");
	  vertical_offset = dataset.openAttribute("vertical_offset");
	  max_value = dataset.openAttribute("max_value");
	  min_value = dataset.openAttribute("min_value");

	  horiz_interval.read(PredType::NATIVE_DOUBLE, &dx);
	  vertical_gain.read(PredType::NATIVE_DOUBLE, &dy);
	  horiz_offset.read(PredType::NATIVE_DOUBLE, &xoffset);
	  vertical_offset.read(PredType::NATIVE_DOUBLE, &yoffset);
	  max_value.read(PredType::NATIVE_DOUBLE, &Vmax);
	  min_value.read(PredType::NATIVE_DOUBLE, &Vmin);

	  DataCluster * datacluster = Init_Data(&dataset); //initialize datacluser
 
	  unsigned long window_length = datacluster->trace_length;
	  unsigned long integrate_length = window_length - window_width; 

	  cout << " the dx is " << dx << ", the integrate length is " << integrate_length << ", the trace length is " << window_length << ", the dy is " << dy << endl;
	  
	  // for trig PMT 
	  double big_noise_charge = 0.0;
	  float big_noise_voltage = 0.0; 
	  float big_noise_pedestal = 0.0; 
	  float big_window_voltage = 0.0; 
	  unsigned long time_bin_counter = 0.0; 
	  float max_time_entry = 0.0;
	
	  for (j = 0; j < datacluster->n_traces; j++){
	    cout << "Analyzing trace " << j+1 << " out of " 
	    	 <<datacluster->n_traces<< " total traces " << " for channel 1 trigger PMT."<<'\r';
	    cout.flush(); 

	    Read_Trace(datacluster,j); 
	    big_noise_pedestal = TMath::Mean (500, datacluster->data_out)*dy; // pedestal correction
	    
	    big_noise_charge = 0.0; 
	    max_time_entry = 0.0; 
	    time_bin_counter = 0.0; 

	    for(i = 500; i < 500 + integrate_length; i++){ // signal window for trigger pmt, arbitrary 
	      big_noise_voltage = ((float)datacluster->data_out[i]*dy-big_noise_pedestal);
	      big_noise_voltages->Fill(big_noise_voltage); 
	      big_noise_charge = big_noise_charge+(big_noise_voltage*((-1000.0*dx*1e9)/termination_ohms)); 
	    }

	    for(i = 0; i < window_length ; i++){
	      big_window_voltage = ((float)datacluster->data_out[i]*dy-big_noise_pedestal);
	      big_waveform_voltage[i] = big_waveform_voltage[i] + big_window_voltage; 
	    } 

	    for(i= 0; i < window_length ; i++){
	      timing_voltage[i] = ((float)datacluster->data_out[i]*dy-big_noise_pedestal);
	      if(timing_voltage[i] < max_time_entry){
		max_time_entry = timing_voltage[i]; 
		time_bin_counter = i; // bin corresponding to maximum voltage 
	      }
	    }
	     
	    // convert bin to time 
	    time_max_trig[j] = time_bin_counter*0.2; // time-bins are 0.2ns long 

	    new_window_count++; 
	    big_noise_pedestals->Fill(big_noise_pedestal); 
	    big_noise_charges->Fill(big_noise_charge); 
	  }

	  file.close(); 

	}
	ifs.close(); 
      }

      Float_t big_avg_waveform_voltage[length_trace];  
      for(i=0; i<length_trace; i++){ 
	big_avg_waveform_voltage[i] = big_waveform_voltage[i]/new_window_count; 
	big_average_waveform->SetBinContent(i+1, big_avg_waveform_voltage[i]);  
      }

      cout  << "\n The number of Traces accepted is " << new_window_count << " for the trigger PMT." << endl;


      TAxis *axis = charges_noise->GetXaxis();

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
      chargefile << "first bin " << first_bin << " last bin " << last_bin << endl;      
      while(first_bin < last_bin){
	charge_finder[f] = charges_noise->GetBinContent(first_bin);
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
      Double_t entries_zero_bin = charges_noise->GetBinContent(zero_bin); 
      Double_t entrie_checker = entries_zero_bin;
      while(entrie_checker >= entries_zero_bin/2.0){
	zero_bin++;
	entrie_checker = charges_noise->GetBinContent(zero_bin); 
      }      
      Double_t charge_half_noise_height = charges_noise->GetBinCenter(zero_bin); 

      // do the fits, don't draw 
      TCanvas *c1 = new TCanvas("c1","Charge",200,10,700,500);
      charges_noise->Fit("gaus","0","",bottom_charge_fit, top_charge_fit);
      TF1 *myfunc=charges_noise->GetFunction("gaus");
      Double_t p0=myfunc->GetParameter(0); 
      Double_t p1=myfunc->GetParameter(1); 
      Double_t p2=myfunc->GetParameter(2); 
      charges_noise->Fit("gaus","0","",-charge_half_noise_height, charge_half_noise_height);
      TF1 *myfunc2=charges_noise->GetFunction("gaus");
      //Double_t p00=myfunc2->GetParameter(0);
      //Double_t p11=myfunc2->GetParameter(1); 
      Double_t p22=myfunc2->GetParameter(2);
      charges_noise->Fit("pol2","0","",5*p22, bottom_charge_fit);
      TF1 *myfunc3=charges_noise->GetFunction("pol2");
      Double_t p000=myfunc3->GetParameter(0); 
      Double_t p111=myfunc3->GetParameter(1); 
      Double_t p222=myfunc3->GetParameter(2);
      charges_noise->Draw(); 
      c1->Update(); 

      axis->SetRange(0.3, 1.0); 
      Double_t valley = charges_noise->GetMinimumBin(); 
      Double_t valley_entries = charges_noise->GetBinContent(valley);
      cout << "peak " << p0 << " valley " << valley_entries << endl; 

      Double_t Low_charge_counter = axis->FindBin(3*p22); // use to count number of entries above 3 times the electronic noise width
      Double_t High_charge_counter = axis->FindBin(p1 + 3*p2); // use to count number of entries 3 sigma above the charge peak 
      Double_t low_charge_entry = 0.0; 
      Double_t high_charge_entry = 0.0; 
      Double_t Low_bin_counter = Low_charge_counter; 
      Double_t High_bin_counter = High_charge_counter; 
           
      for(i = Low_bin_counter; i < Low_bin_counter + 1000; i++)
	{
	  low_charge_entry = low_charge_entry + charges_noise->GetBinContent(Low_charge_counter); 
	  Low_charge_counter++;
	}

      for(i = High_bin_counter; i < High_bin_counter + 1000; i++)
	{
	  high_charge_entry = high_charge_entry + charges_noise->GetBinContent(High_charge_counter);
	  High_charge_counter++;
	}
      
      // make the statistics legend
      TPaveStats *ptstats = new TPaveStats(0.6,0.6,0.98,0.98,"brNDC");
      ptstats->SetName("stats");
      ptstats->SetBorderSize(2);
      ptstats->SetTextAlign(12);
      ptstats->SetTextFont(42);
      ptstats->SetTextSize(0.035);
      ptstats->SetShadowColor(0);
      
      Double_t FWHM =  p2*2*sqrt(2*log(2));
      Double_t HCT = high_charge_entry * 100 / low_charge_entry;
      Double_t elec_width = p22; 
      Double_t Peak_to_valley = p0/valley_entries; 
      
      TString box_title = "Charge"; 
      TText *text_t = ptstats->AddText(box_title);
      text_t->SetTextSize(0.05); 
      TString ent = Form("Entries = %.0f", charges_noise->GetEntries() );
      text_t = ptstats->AddText(ent);
      TString noise_t = Form("Electronics Noise Width = %.3f pC", elec_width);
      text_t = ptstats->AddText(noise_t); 
      TString peak_t = Form("Charge Peak = %.3f pC", p1);
      text_t = ptstats->AddText(peak_t);
      TString FWHM_t = Form("Charge FWHM = %.3f pC", FWHM); 
      text_t = ptstats->AddText(FWHM_t);
      TString HCT_t = Form("High Charge Tail = %.3f%%", HCT); 
      text_t = ptstats->AddText(HCT_t); 

      ptstats->SetOptStat(0);
      ptstats->SetOptFit(0);
      ptstats->SetFillColor(0);

      charges_noise->GetListOfFunctions()->Add(ptstats);    
      charges_noise->GetXaxis()->SetTitle("Charge (pC)");
      charges_noise->GetXaxis()->SetTitleFont(42);
      charges_noise->GetXaxis()->SetTitleColor(1);
      charges_noise->GetXaxis()->SetLabelFont(42);
      charges_noise->GetXaxis()->SetRangeUser(-1.0,5.0);
      charges_noise->GetYaxis()->SetTitle("Events");
      charges_noise->GetYaxis()->SetTitleOffset(1.2);
      charges_noise->GetYaxis()->SetLabelFont(42);
      charges_noise->GetYaxis()->SetTitleFont(42);
      c1->Modified();
	  
      cout << "Electronic Noise Width is " << p22 << endl; 
      cout << "The Charge Peak is " << p1 << endl; 
      cout << "The Charge FWHM is " << p2*2*sqrt(2*log(2)) << endl; 
      cout << "High Charge Tail " << high_charge_entry * 100 / low_charge_entry << "%" << endl; 
      

      // loop for timing 
      for(a=0; a<1; a++){

	ifstream ifs ( argv[1] , ifstream::in ); // Open File List Stream
	ifs >> filename;

	while (ifs.good()){ // While Stream Is Open, Analyze New Files.

	  cout<<filename<<endl;
	  file.openFile(filename, H5F_ACC_RDONLY); //Open HDF5 File 
	  ifs >> filename;
	  dataset = file.openDataSet("channel3"); // Open HDF5 DataSet

	  // Reading in Attributes
	  horiz_interval = dataset.openAttribute("horiz_interval");
	  vertical_gain = dataset.openAttribute("vertical_gain");
	  horiz_offset = dataset.openAttribute("horiz_offset");
	  vertical_offset = dataset.openAttribute("vertical_offset");
	  max_value = dataset.openAttribute("max_value");
	  min_value = dataset.openAttribute("min_value");

	  horiz_interval.read(PredType::NATIVE_DOUBLE, &dx);
	  vertical_gain.read(PredType::NATIVE_DOUBLE, &dy);
	  horiz_offset.read(PredType::NATIVE_DOUBLE, &xoffset);
	  vertical_offset.read(PredType::NATIVE_DOUBLE, &yoffset);
	  max_value.read(PredType::NATIVE_DOUBLE, &Vmax);
	  min_value.read(PredType::NATIVE_DOUBLE, &Vmin);
	  
	  DataCluster * datacluster = Init_Data(&dataset); //initialize datacluser
	
	  unsigned long window_length = datacluster->trace_length; 
	  unsigned long integrate_length = window_length - window_width; 
	  cout << " the dx is " << dx << ", the integrate length is " << integrate_length << ", the trace length is " << window_length << ", the dy is " << dy << endl;

	  
	  double time_charge = 0.0; 
	  float time_voltage = 0.0; 
	  float time_pedestal = 0.0;
	  float max_time_entry2 = 0.0; 
	  unsigned long time_bin_counter2 = 0.0;
	  //float peak = 0.0;

	  for (j = 0; j < datacluster->n_traces; j++){
	    cout << "Analyzing trace " << j+1 << " out of " 
	    	 <<datacluster->n_traces<< " total traces " << " for channel 3 ETL PMT"<<'\r';
	    cout.flush(); 

	    time_pedestal = TMath::Mean (window_width, datacluster->data_out)*dy;		    
	    Read_Trace(datacluster,j);
	    
	    time_charge = 0.0; 
	    max_time_entry2 = 0.0; 
	    time_bin_counter2 = 0.0; 
	    
	    for(i = window_width; i < window_width + integrate_length ; i++){ // noise window
	      time_voltage = ((float)datacluster->data_out[i]*dy-time_pedestal);
	      time_charge = time_charge+(time_voltage*((-1000.0*dx*1e9)/termination_ohms));		  
	    }
	    	    
	    if(time_charge > 3*elec_width){
	      for(i = window_width; i < window_width + integrate_length ; i++){ // noise window
		timing_voltage2[i] = ((float)datacluster->data_out[i]*dy-time_pedestal);
		if(timing_voltage2[i] < max_time_entry2){
		  max_time_entry2 = timing_voltage2[i];
		  time_bin_counter2 = i; 
		}
	      }
	    }
	    
	    //peak = ((float)datacluster->data_out[time_bin_counter2]*dy-time_pedestal);

	    // convert bin to time 
	    time_max_trig2[j] =  time_bin_counter2*0.2; 

	  }
	  file.close(); 
	}
	ifs.close(); 
      }
      
      // get time difference
      for(i=0; i < number_entries; i++){
	if(time_max_trig2[i] != 0){
	  delta_T[i] = time_max_trig2[i] - time_max_trig[i];
	  chargefile << "time max trig1 " << time_max_trig[i] << " time max trig2" << time_max_trig2[i] << " time difference " << delta_T[i] << " for bin " << i << endl; 
	  Timing->Fill(delta_T[i]); 
	
	}
      }

      float time_bin_max = Timing->GetMaximumBin(); 
      TAxis *time_axis = Timing->GetXaxis();
      Double_t time_max = time_axis->GetBinCenter(time_bin_max);

      // draw the fits 
      TCanvas *c2 = new TCanvas("c2","Charge",200,10,700,500);
      Timing->Fit("gaus","0","",time_max - 3.0, time_max + 3.0);
      TF1 *fitfunc=Timing->GetFunction("gaus");
      Double_t f1=fitfunc->GetParameter(1); cout << " Prompt Mean " << f1 << endl; 
      Double_t f2=fitfunc->GetParameter(2); cout << " Prompt Sigma " << f2 << endl; 
      Timing->Draw(); 
      c2->Update(); 

      // get coincidence rate 
      Double_t coincidence_count = 0.0;
      for(i=0; i < number_entries; i++){
	if(time_max_trig2[i] != 0){
	  if(delta_T[i] < f1 + 3*f2 && delta_T[i] > f1 - 3*f2){ // if delta t is within 4 sigma of prompt peak 
	    coincidence_count++;
	  }
	}
      }
    
      // make the statistics legend
      TPaveStats *ttstats = new TPaveStats(0.6,0.6,0.98,0.98,"brNDC");
      ttstats->SetName("stats");
      ttstats->SetBorderSize(2);
      ttstats->SetTextAlign(12);
      ttstats->SetTextFont(42);
      ttstats->SetTextSize(0.035);
      ttstats->SetShadowColor(0);

      Double_t Prompt_sigma = f2; 
      Double_t Prompt_FWHM = 2*f2*sqrt(2*log(2));
      Double_t Coincidence_Percent = coincidence_count/number_entries*100; // depends on the hard coded number of entries, shitty 
      
      // statistics for timing histogram 
      TString box_title2 = "Timing"; 
      TText *text_tt = ttstats->AddText(box_title2);
      text_tt->SetTextSize(0.05); 
      TString t_ent = Form("Hits above noise = %.0f", Timing->GetEntries() );
      text_tt = ttstats->AddText(t_ent);
      TString pulse_sigma_t = Form("Prompt Sigma = %.3f ns", Prompt_sigma);
      text_tt = ttstats->AddText(pulse_sigma_t); 
      TString pulse_fwhm_t = Form("Prompt FWHM = %.3f ns", Prompt_FWHM);
      text_tt = ttstats->AddText(pulse_fwhm_t); 
      TString coincidence_t = Form("Prompt Coincidence Rate = %.3f%%", Coincidence_Percent);
      text_tt = ttstats->AddText(coincidence_t); 

      ttstats->SetOptStat(0);
      ttstats->SetOptFit(0);
      ttstats->SetFillColor(0);

      Timing->GetListOfFunctions()->Add(ttstats);  
      Timing->GetXaxis()->SetTitle("#Deltat (ns)");
      Timing->GetXaxis()->SetTitleFont(42);
      Timing->GetXaxis()->SetTitleColor(1);
      Timing->GetXaxis()->SetLabelFont(42);
      Timing->GetYaxis()->SetTitle("Events");
      Timing->GetYaxis()->SetLabelFont(42);
      Timing->GetYaxis()->SetTitleFont(42);
      c2->Modified();
      	  
      chargefile.close(); 

      // Output Histograms to File
      if (argc > 2){
	TFile f(argv[2],"new");
	//ETL PMT
	noise_pedestals->Write();
	corrected_noise_pedestals->Write();
	voltages->Write(); 
	noise_voltages->Write(); 
	variances->Write();
	noise_variances->Write();
	charges_noise->Write();
	acceptable_waveform->Write();
	pulse_waveform->Write();
	average_waveform->Write(); 
	  
	//Trig PMT 
	big_noise_pedestals->Write(); 
	big_noise_voltages->Write();
	big_noise_charges->Write(); 
	big_average_waveform->Write(); 

	Timing->Write();
	
      }
  
      // clean-up 
      delete noise_pedestals;
      delete corrected_noise_pedestals; 
      delete voltages;
      delete noise_voltages; 
      delete variances; 
      delete noise_variances; 
      delete charges_noise;   
      delete acceptable_waveform; 
      delete pulse_waveform;
      delete average_waveform; 

      delete big_noise_pedestals; 
      delete big_noise_voltages; 
      delete big_noise_charges; 
      delete big_average_waveform; 

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
