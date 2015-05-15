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
#include <TMinuit.h>
#include <TFitter.h>

// HDF5 Library
#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

using namespace std;

const int RANK_OUT = 2; // Data Rank

const unsigned long length_trace = 10002;  
const unsigned long number_entries = 200000; 
const unsigned long number_files = 19; 
const unsigned long number_windows = 4; 
const int high_voltage = 2085; 
//const unsigned long window_start = 499;
const int length_late = 50; 
const float variance_cut = 0.0006; //0.0002; 
const float charge_cut = 0.8; 
//**** FOR EFFICIENCY ***
// ETL105 = window_start = 699, variance_cut = 0.00019, length_late = 95, HV 1340 
// ETL142 = window_start = 649 
// ZN0116 - window_start = 499 
// ZT0070 - window_start = 399 , length late 80 
// SNO tube - window_Start 1699, length late 50, variance cut 0.0004, HV 2024

//TString box_title = Form("ZN0116 12-inch HQE"); 
//TString box_title = Form("ZT0072 10-inch HQE");
//TString box_title = Form("11 Inch ET PMT: D784 KFLB, 134"); 
TString box_title = Form("8 inch SNO PMT");					

static double time_max_trig[number_files][number_entries] = {{ 0 }}; 
static double time_max_trig2[number_files][number_entries][number_windows] = {{{ 0 }}}; 
static double delta_T[number_files][number_entries][number_windows] = {{{ 0 }}}; 

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
//static double spe_CrystalBall(double *x,double *par);


int main (int argc, char* argv[])
{
  unsigned long window_width = atoi(argv[3]);
  //const unsigned long high_voltage = atoi(argv[4]);
  //TString box_title = argv[5]; 
  const int termination_ohms = 50; 

  try
    {

      cout << box_title << endl; 
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

      //Define ROOT histograms here, the variable names are confusing right now, look at the descriptions 
      TH1F *noise_pedestals = new TH1F("Pedestal for ETL PMT","",100000,-1,1); 
      TH1F *corrected_noise_pedestals = new TH1F("Corrected Pedestal for ETL PMT","",100000,-1,1);
      TH1F *voltages = new TH1F("Voltages over Pedestal Window for ETL PMT","",5000,-1.0,1.0);
      TH1F *noise_voltages = new TH1F("Voltages over Signal Window for ETL PMT","",5000,-1.0,1.0);
      TH1F *variances = new TH1F("Variance over Pedestal Window for ETL PMT","",500000,-0.1,1.0);
      TH1F *charges_signal = new TH1F("Charge ETL PMT","",20000,-150,200);   
      TH1F *acceptable_waveform = new TH1F("Acceptable Waveform ETL PMT","", length_trace, 0, length_trace*0.1); 
      TH1F *pulse_waveform = new TH1F("Unacceptable Waveform ETL PMT","", length_trace, 0, length_trace*0.1); 
      TH1F *average_waveform = new TH1F("Average Waveform ETL PMT","", length_trace, 0, length_trace*0.1);  
      
      TH1F *big_noise_pedestals = new TH1F("Pedestal for Trigger PMT","",100000,-1.0,1.0); 
      TH1F *big_noise_voltages =  new TH1F("Voltages over Signal Window for Trigger PMT","",5000,-1.0,1.0);
      TH1F *big_noise_charges =  new TH1F("Signal Charge for Trigger PMT","",5000,-150,200); 
      TH1F *big_average_waveform = new TH1F("Average Waveform for Trigger PMT","", length_trace, 0, length_trace*0.1); 

      TH1F *Timing = new TH1F("Timing histogram","",(length_trace/2)*0.1, 0, (length_trace/2)*0.1);  
      TH1F *check_average_waveform = new TH1F("Used to find min bin","", length_trace, 0, length_trace*0.1);  

      
      // seperates file input 
      int g = 0; 
      int p = 0;
      int z = 0; 
         
      unsigned long i; 
      unsigned long j; 
      unsigned long k;
      unsigned int window_count = 0.0;
      unsigned int new_window_count = 0.0;  
      unsigned int check_window_count = 0.0;
      string filename; 
      H5File file; 
      DataSet dataset;

      // get rates
      Double_t coincidence_count = 0.0;
      Double_t dark_count = 0.0; 
      Double_t late_pulse_count = 0.0; 
      Double_t pre_pulse_count = 0.0;
      Double_t post_pulse_count = 0.0;

      int r = 0.0; // file count

      // average waveform

      float waveform_voltage[length_trace] = {0};
      float waveform_voltage2[length_trace] = {0};
      float window_voltage[length_trace] = {0};
      float window_voltage2[length_trace] = {0};
      float big_waveform_voltage[length_trace] = {0};
      float timing_voltage[length_trace] = {0};
      float timing_voltage2[length_trace] = {0};
      float timing_voltage4[length_trace] = {0};
      float timing_voltage6[length_trace] = {0};
      float timing_voltage8[length_trace] = {0};

	
      for(z=0; z<1; z++){
	
	ifstream ifs ( argv[1] , ifstream::in ); // Open File List Stream
	ifs >> filename;

	while (ifs.good()){ // While Stream Is Open, Analyze New Files.

	  cout<<filename<<endl;
	  file.openFile(filename, H5F_ACC_RDONLY); //Open HDF5 File 
	  ifs >> filename;
	  dataset = file.openDataSet("channel2"); // Open HDF5 DataSet

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
	  cout << "The time bin width is " << dx << ", the window width is " << window_width << ", the integrate length is " << integrate_length << ", the trace length is " << window_length << ", the vertical resolution is " << dy << "." << endl;

	  float check_pedestal = 0.0;
	  float fix_check_pedestal = 0.0;
	  float c_variance = 0.0;
	  float c_voltage = 0.0;

	  for (j = 0; j < datacluster->n_traces; j++){
	    cout << "Analyzing trace " << j+1 << " out of " 
		 <<datacluster->n_traces<< " total traces " << "for channel 3 ETL"<<'\r';
	    cout.flush(); 

	    Read_Trace(datacluster,j);
	    check_pedestal = TMath::Mean (window_width, datacluster->data_out)*dy;
	    c_variance = 0.0;

	    for(i=0; i < window_width; i++){
	      c_voltage = ((float)datacluster->data_out[i]*dy-check_pedestal);
	      c_variance = c_variance + c_voltage * c_voltage; 
	    }
	    
	    if(c_variance < variance_cut) { // traces that don't pass this cut have dark pulses! 
	      fix_check_pedestal = TMath::Mean (window_width, datacluster->data_out)*dy;
	      check_window_count++;
	    }

	    for(i=0; i < window_length; i++){ // entire window 
	      window_voltage2[i] = ((float)datacluster->data_out[i]*dy-fix_check_pedestal);  
	      waveform_voltage2[i] = waveform_voltage2[i] + window_voltage2[i]; 
	    }

	  }
	  file.close(); 
	}
	ifs.close(); 
      }
      
      cout  << "\n The number of Traces accepted is " << check_window_count << " for the test PMT." << endl;
      
      Float_t check_avg_waveform_voltage[length_trace]; 
      for(i=0; i<length_trace; i++){ 
	check_avg_waveform_voltage[i] = waveform_voltage2[i]/check_window_count;
	check_average_waveform->SetBinContent(i+1, check_avg_waveform_voltage[i]);  
      }

      unsigned long window_start = check_average_waveform->GetMinimumBin() - length_trace/2 - 120;

      cout << "Window start " << window_start << endl;

      ofstream chargefile;
      chargefile.open("info.text"); // used mostly for debugging 
	  
      for(g=0; g<1; g++){

	ifstream ifs ( argv[1] , ifstream::in ); // Open File List Stream
	ifs >> filename;

	while (ifs.good()){ // While Stream Is Open, Analyze New Files.

	  cout<<filename<<endl;
	  file.openFile(filename, H5F_ACC_RDONLY); //Open HDF5 File 
	  ifs >> filename;
	  dataset = file.openDataSet("channel2"); // Open HDF5 DataSet

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
	  cout << "The time bin width is " << dx << ", the window width is " << window_width << ", the integrate length is " << integrate_length << ", the trace length is " << window_length << ", the vertical resolution is " << dy << "." << endl;
	  
	  unsigned long prompt_window_start = window_width + window_start;
		  
	  // for PMT1
	  double ncharge = 0.0; 
	  double tcharge = 0.0;
	  double tcharge2 = 0.0;
	  double tcharge3 = 0.0;
	  float voltage = 0.0; 
	  float noise_voltage = 0.0; 
	  float test_voltage = 0.0;
	  float test_voltage2 = 0.0;
	  float test_voltage3 = 0.0;
	  float variance = 0.0; 
	  float noise_pedestal = 0.0; 
	  float corrected_noise_pedestal = 0.0; 

	  float timing_voltage3 = 0.0; 
	  float timing_voltage5 = 0.0;
	  float timing_voltage7 = 0.0;
	  float timing_voltage9 = 0.0;

	  float max_time_entry2 = 0.0; 
	  float max_time_entry3 = 0.0;
	  float max_time_entry4 = 0.0;
	  float max_time_entry5 = 0.0;

	  unsigned long time_bin_counter2 = 0.0;
	  unsigned long time_bin_counter3 = 0.0;
	  unsigned long time_bin_counter4 = 0.0;
	  unsigned long time_bin_counter5 = 0.0;
	  unsigned long time_bin_counter6 = 0.0;
	  unsigned long time_bin_counter7 = 0.0;
	  unsigned long time_bin_counter8 = 0.0;
	  unsigned long time_bin_counter9 = 0.0;

	  float prompt_peak = 0.0; 
	  float prompt_peak2 = 0.0;
	  float prompt_peak3 = 0.0;
	  float prompt_peak4 = 0.0;

	  for (j = 0; j < datacluster->n_traces; j++){
	    cout << "Analyzing trace " << j+1 << " out of " 
		 <<datacluster->n_traces<< " total traces " << "for channel 3 ETL"<<'\r';
	    cout.flush(); 

	    Read_Trace(datacluster,j); 	    	    
	    noise_pedestal = TMath::Mean (window_width, datacluster->data_out)*dy;
	  
	    // initialize 
	    ncharge = 0.0;
	    tcharge = 0.0;
	    tcharge2 = 0.0; 
	    tcharge3 = 0.0;
	    variance = 0.0; 

	    max_time_entry2 = 0.0; 
	    max_time_entry3 = 0.0;
	    max_time_entry4 = 0.0;
	    max_time_entry5 = 0.0; 

	    time_bin_counter2 = 0.0; 
	    time_bin_counter3 = 0.0;

	    time_bin_counter4 = 0.0; 
	    time_bin_counter5 = 0.0;

	    time_bin_counter6 = 0.0; 
	    time_bin_counter7 = 0.0;	    

	    time_bin_counter8 = 0.0;	    
	    time_bin_counter9 = 0.0;	    

	    prompt_peak = 0.0; 
	    prompt_peak2 = 0.0;
	    prompt_peak3 = 0.0; 
	    prompt_peak4 = 0.0;

	    k = 0; 
	    
	    for(i=0; i < window_width; i++){
	      voltage = ((float)datacluster->data_out[i]*dy-noise_pedestal);
	      voltages->Fill(voltage);
	      variance = variance + voltage * voltage; 
	    }
	    
	    if(variance < variance_cut) { // traces that don't pass this cut have dark pulses! 
	      

	      corrected_noise_pedestal = TMath::Mean (window_width, datacluster->data_out)*dy; 

	      // PROMPT WINDOW 

	      for(i = prompt_window_start; i < prompt_window_start + 300; i++){ // signal window
		noise_voltage = ((float)datacluster->data_out[i]*dy-corrected_noise_pedestal);
		noise_voltages->Fill(noise_voltage); 
		ncharge = ncharge+(noise_voltage*((-1000.0*dx*1e9)/termination_ohms)); // charge 
	      }

	      if(ncharge > charge_cut){ // hard cut on what is accepted as a prompt pulse 
		for(i = prompt_window_start; i < prompt_window_start + 300; i++){ // signal window
		  timing_voltage2[i] = ((float)datacluster->data_out[i]*dy-corrected_noise_pedestal);
		  if(timing_voltage2[i] < max_time_entry2){
		    max_time_entry2 = timing_voltage2[i];
		    time_bin_counter2 = i; // gets the bin of max peak
		  }
		}
	      }
	    
	      prompt_peak = ((float)datacluster->data_out[time_bin_counter2]*dy-corrected_noise_pedestal); // gets voltage of max peak
	      
	      if(ncharge > charge_cut){
		for(i = time_bin_counter2 - 150; i < time_bin_counter2; i++){ // window 15ns before the max peak
		  timing_voltage3 = ((float)datacluster->data_out[i]*dy-corrected_noise_pedestal); // get voltage
		  if(timing_voltage3 > prompt_peak*0.2){ // if voltage is less negative than 20% of promp peak
		    time_bin_counter3 = i; 
		  }
		}
	      }
		
	      // convert bin to time 
	      time_max_trig2[r][j][k] =  time_bin_counter3*0.1;
	      k++;
	    
	      // DARK WINDOW
	      
	      for(i = window_width; i < prompt_window_start; i++){ // dark window
		test_voltage3 = ((float)datacluster->data_out[i]*dy-corrected_noise_pedestal); 
		tcharge3 = tcharge3+(test_voltage3*((-1000.0*dx*1e9)/termination_ohms)); // charge 
	      }
	      
	      if(tcharge3 > charge_cut){
		for(i = window_width; i < prompt_window_start; i++){ // dark window
		  timing_voltage4[i] = ((float)datacluster->data_out[i]*dy-corrected_noise_pedestal);
		  if(timing_voltage4[i] < max_time_entry3){
		    max_time_entry3 = timing_voltage4[i];
		    time_bin_counter4 = i; // gets the bin of max peak
		  }
		}
	      }

	      prompt_peak2 = ((float)datacluster->data_out[time_bin_counter4]*dy-corrected_noise_pedestal); // gets voltage of max peak

	      if(tcharge3 > charge_cut){
		for(i = time_bin_counter4 - 150; i < time_bin_counter4; i++){ // window 15ns before the max peak
		  timing_voltage5 = ((float)datacluster->data_out[i]*dy-corrected_noise_pedestal); // get voltage
		  if(timing_voltage5 > prompt_peak2*0.2){ // if voltage is less negative than 20% of promp peak
		    time_bin_counter5 = i; 
		  }
		}
	      }
	      
	      // convert bin to time 
	      time_max_trig2[r][j][k] =  time_bin_counter5*0.1;
	      k++;

	      // LATE WINDOW 

	      for(i = prompt_window_start + 300; i < prompt_window_start + 1100; i++){
		test_voltage = ((float)datacluster->data_out[i]*dy-corrected_noise_pedestal); 
		tcharge = tcharge+(test_voltage*((-1000.0*dx*1e9)/termination_ohms)); // charge 
	      }

	      if(tcharge > charge_cut){
		for(i = prompt_window_start + 300; i < prompt_window_start + 1100; i++){ // dark window
		  timing_voltage6[i] = ((float)datacluster->data_out[i]*dy-corrected_noise_pedestal);
		  if(timing_voltage6[i] < max_time_entry4){
		    max_time_entry4 = timing_voltage6[i];
		    time_bin_counter6 = i; // gets the bin of max peak
		  }
		}
	      }

	      prompt_peak3 = ((float)datacluster->data_out[time_bin_counter6]*dy-corrected_noise_pedestal); // gets voltage of max peak

	      if(tcharge > charge_cut){
		for(i = time_bin_counter6 - 150; i < time_bin_counter6; i++){ // window 15ns before the max peak
		  timing_voltage7 = ((float)datacluster->data_out[i]*dy-corrected_noise_pedestal); // get voltage
		  if(timing_voltage7 > prompt_peak3*0.2){ // if voltage is less negative than 20% of promp peak
		    time_bin_counter7 = i; 
		  }
		}
	      }
	      
	      // convert bin to time 
	      time_max_trig2[r][j][k] =  time_bin_counter7*0.1;
	      k++;

	      // LAST DARK WINDOW

	      for(i = prompt_window_start + 1100; i < window_length; i++){
		test_voltage2 = ((float)datacluster->data_out[i]*dy-corrected_noise_pedestal); 
		tcharge2 = tcharge2 +(test_voltage2*((-1000.0*dx*1e9)/termination_ohms)); // charge 
	      }

	      if(tcharge2 > charge_cut){
		for(i = prompt_window_start + 1100; i < window_length; i++){ // dark window
		  timing_voltage8[i] = ((float)datacluster->data_out[i]*dy-corrected_noise_pedestal);
		  if(timing_voltage8[i] < max_time_entry5){
		    max_time_entry5 = timing_voltage8[i];
		    time_bin_counter8 = i; // gets the bin of max peak
		  }
		}
	      }

	      prompt_peak4 = ((float)datacluster->data_out[time_bin_counter8]*dy-corrected_noise_pedestal); // gets voltage of max peak

	      if(tcharge2 > charge_cut){
		for(i = time_bin_counter8 - 150; i < time_bin_counter8; i++){ // window 15ns before the max peak
		  timing_voltage9 = ((float)datacluster->data_out[i]*dy-corrected_noise_pedestal); // get voltage
		  if(timing_voltage9 > prompt_peak4*0.2){ // if voltage is less negative than 20% of promp peak
		    time_bin_counter9 = i; 
		  }
		}
	      }
	      
	      // convert bin to time 
	      time_max_trig2[r][j][k] =  time_bin_counter9*0.1;
	      
	      // END OF WINDOWS

	      
	      for(i=0; i < window_length; i++){ // entire window 
		window_voltage[i] = ((float)datacluster->data_out[i]*dy-corrected_noise_pedestal);  
		waveform_voltage[i] = waveform_voltage[i] + window_voltage[i]; 
	      }   

	      for(i=0; i < window_length; i++){ // draws an acceptable trace 
		acceptable_waveform->SetBinContent(i+1,(float)datacluster->data_out[i]*dy-corrected_noise_pedestal); 
	      }
 

	      window_count++; // number of traces 
	  
	      charges_signal->Fill(ncharge); 
	      corrected_noise_pedestals->Fill(corrected_noise_pedestal); 
	    }

	    if(variance > 0.001){ // draw a waveform that was cut, should have a dark pulse 
	      for(i=0; i<window_length; i++){
		pulse_waveform->SetBinContent(i+1,(float)datacluster->data_out[i]*dy-noise_pedestal);
	      }
	    }
	    
	    //store histograms 
	    variances->Fill(variance); 
	    noise_pedestals->Fill(noise_pedestal);
	  }
	  r++;
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
      
      r = 0.0;
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

	  cout << "The file number is " << r << ", the time bin width is " << dx << ", the integrate length is " << integrate_length << ", the trace length is " << window_length << ", the vertical resolution is " << dy << "." << endl;

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

	    for(i = window_length/2 - 60; i < window_length/2 + 60; i++){ // around the trigger pulse 
	      timing_voltage[i] = ((float)datacluster->data_out[i]*dy-big_noise_pedestal);
	      if(timing_voltage[i] < max_time_entry){
		max_time_entry = timing_voltage[i]; 
		time_bin_counter = i; // bin corresponding to maximum voltage 
	      }
	    }
	     
	    // convert bin to time 
	    time_max_trig[r][j] = time_bin_counter*0.1; // time-bins are 0.2ns long , which is the maximum horizontal resolution for this scope 

	    new_window_count++; 
	    big_noise_pedestals->Fill(big_noise_pedestal); 
	    big_noise_charges->Fill(big_noise_charge); 
	  }
	  r++; // count files 
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
      cout  << " The number of files is " << r  << endl; 

      
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
	cout << zero_bin << endl; 
	entrie_checker = charges_signal->GetBinContent(zero_bin); 
      }      
      Double_t charge_half_noise_height = charges_signal->GetBinCenter(zero_bin); 

      // do the fits, don't draw 
      TCanvas *c1 = new TCanvas("c1","Charge",200,10,700,500);
      charges_signal->Fit("gaus","0","",-charge_half_noise_height, charge_half_noise_height);
      TF1 *myfunc2=charges_signal->GetFunction("gaus");
      //Double_t p00=myfunc2->GetParameter(0);
      //Double_t p11=myfunc2->GetParameter(1); 
      Double_t p22=myfunc2->GetParameter(2);
      charges_signal->Fit("pol2","0","", 6*p22, bottom_charge_fit);
      TF1 *myfunc3=charges_signal->GetFunction("pol2");
      Double_t p000=myfunc3->GetParameter(0); 
      Double_t p111=myfunc3->GetParameter(1); 
      Double_t p222=myfunc3->GetParameter(2); 
      /*
	TF1 *spefunction=new TF1("spefunction",spe_CrystalBall,1.0,3.0,5);
      
	spefunction->SetParameters(1.6, -0.8, -2.0, 1.8, 100.0); 
	spefunction->SetParLimits(0, 1.5, 1.7);
	spefunction->SetParLimits(1, -2, 2);
	spefunction->SetParLimits(2, -10, 10);
	spefunction->SetParLimits(3, 0, 10);
	spefunction->SetParLimits(4, 0, 200);*/
      
      /*Double_t par[9]; 
      TF1 *g1 = new TF1("g1","gaus",bottom_charge_fit,top_charge_fit);
      TF1 *g2 = new TF1("g2","expo",0.5,bottom_charge_fit);
      TF1 *g3 = new TF1("g3","pol4",top_charge_fit,6);
      TF1 *total = new TF1("total","gaus(0)+expo(3)+pol4(5)",0.5,6);
      total->SetLineColor(2);
      charges_signal->Fit(g1,"R");
      charges_signal->Fit(g2,"R+");
      charges_signal->Fit(g3,"R+");
      g1->GetParameters(&par[0]);
      g2->GetParameters(&par[3]);
      g3->GetParameters(&par[5]);
      total->SetParameters(par);
      charges_signal->Fit(total,"R+");*/
      charges_signal->Fit("gaus","","",bottom_charge_fit, top_charge_fit);
      TF1 *myfunc = charges_signal->GetFunction("gaus"); 
      //TF1 *bottom_myfunc = charges_signal->GetFunction("exp");
      //TF1 *top_myfunc = charges_signal->GetFunction("pol4");
      //Double_t spe_mean = g1->GetParameter(1); cout << spe_mean  << endl;
      //Double_t spe_sigma = g1->GetParameter(2); cout << spe_sigma << endl;
      //Double_t gaus_norm = g2->GetParameter(0); cout << gaus_norm << endl; 
      //charges_signal->Fit("spefunction","","",1.0, 3.0);
      //charges_signal->Fit("exp","","",0.5, bottom_charge_fit);
      //Double_t expo_p0 = g2->GetParameter(0); cout << expo_p0 << endl;
      //Double_t expo_p1 = g2->GetParameter(1); cout << expo_p1 << endl;
      //charges_signal->Fit("pol4","","",top_charge_fit,5.0);
 
      //Double_t poly_p1 = g3->GetParameter(0); cout << poly_p1 << endl;
      //Double_t poly_p2 = g3->GetParameter(1); cout << poly_p2 << endl;
      //Double_t poly_p3 = g3->GetParameter(2); cout << poly_p3 << endl;
      //Double_t poly_p4 = g3->GetParameter(3); cout << poly_p4 << endl;
      //TF1 *myfunc=charges_signal->GetFunction("spefunction"); 


      /*
    
	TFitter *minuit= new TFitter(5);
	{
	double p1=-1;
	minuit->ExecuteCommand("SET PRINTOUT",&p1,1);
	minuit->ExecuteCommand("Set NOWARNINGS",&p1,0);
	}
      minuit->SetFCN(spe_CrystalBall);
      minuit->SetParameter(0,"spe",1.6,1,0,0);
      minuit->SetParameter(1,"sigma",1,1,0,0);
      minuit->SetParameter(2,"alpha",2,1,0,0);
      minuit->SetParameter(3,"n",1.8,1,0,0);
      minuit->SetParameter(4,"Norm",100,1,0,0);
      minuit->ExecuteCommand("SCAN",0,0);
      minuit->ExecuteCommand("MINIMIZE",0,0);*/
     
      Double_t spe_peak= myfunc->GetParameter(0); cout << spe_peak << endl;  
      Double_t spe_mean=myfunc->GetParameter(1); cout << spe_mean << endl; 
      Double_t spe_sigma=myfunc->GetParameter(2); cout << spe_sigma << endl; 
      /*Double_t alpha=myfunc->GetParameter(2); cout << alpha << endl;
      Double_t n=myfunc->GetParameter(3); cout << n << endl;
      Double_t N=myfunc->GetParameter(4); cout << N << endl;  */
      /*
      Double_t spe_mean = minuit->GetParameter(0); cout << spe_mean << endl; 
      Double_t spe_sigma=minuit->GetParameter(1); cout << spe_sigma << endl; 
      Double_t alpha=minuit->GetParameter(2); cout << alpha << endl;
      Double_t n=minuit->GetParameter(3); cout << n << endl;
      Double_t N=minuit->GetParameter(4); cout << N << endl;*/

      charges_signal->Draw(); 
      c1->Update(); 

      Float_t x = 0.0; 
      Float_t charge_min = 0.0;
      Double_t min_function = 0.0; // = p000 + p111*x + p222*x**2; 
      Double_t d_function[20000] = {0.0}; 
      Int_t s = 0.0; 

      while(x < 2.0){
	d_function[s] = p111 + 2*p222*x; 
	//chargefile << d_function[s] << endl; 
	if(d_function[s-1] < 0.0 && d_function[s] > 0.0){
	  min_function = p000 + p111*x + p222*x*x;
	  charge_min = x; 
	  //cout << "min_function " << min_function << ", charge min " << charge_min << endl; 
	}
	x += 0.001;
	s++; 
      }    

      Double_t Low_charge_counter = axis->FindBin(3*p22); // use to count number of entries above 3 times the electronic noise width
      Double_t High_charge_counter = axis->FindBin(spe_mean + 3*spe_sigma); // use to count number of entries 3 sigma above the charge peak 
      Double_t low_charge_entry = 0.0; 
      Double_t high_charge_entry = 0.0; 
      Double_t Low_bin_counter = Low_charge_counter; 
      Double_t High_bin_counter = High_charge_counter; 
           
      for(i = Low_bin_counter; i < Low_bin_counter + 1000; i++)
	{
	  low_charge_entry = low_charge_entry + charges_signal->GetBinContent(Low_charge_counter); 
	  Low_charge_counter++;
	}

      for(i = High_bin_counter; i < High_bin_counter + 1000; i++)
	{
	  high_charge_entry = high_charge_entry + charges_signal->GetBinContent(High_charge_counter);
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
      
      Double_t FWHM =  spe_sigma*2*sqrt(2*log(2));
      Double_t HCT = high_charge_entry * 100 / low_charge_entry;
      Double_t elec_width = p22; 
      Double_t Peak_to_valley = spe_peak/min_function; 
      
      
      TText *text_t = ptstats->AddText(box_title);
      text_t->SetTextSize(0.05); 
      TString hv_c = Form("Operating Voltage = %d V",  high_voltage);
      text_t = ptstats->AddText(hv_c);
      TString ent = Form("Entries = %.0f", charges_signal->GetEntries() );
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

      charges_signal->GetListOfFunctions()->Add(ptstats);    
      charges_signal->GetXaxis()->SetTitle("Charge (pC)");
      charges_signal->GetXaxis()->SetTitleFont(42);
      charges_signal->GetXaxis()->SetTitleColor(1);
      charges_signal->GetXaxis()->SetLabelFont(42);
      charges_signal->GetXaxis()->SetRangeUser(-1.0,6.0);
      charges_signal->GetYaxis()->SetRangeUser(0.0,1200.0); 
      charges_signal->GetYaxis()->SetTitle("Events");
      charges_signal->GetYaxis()->SetTitleOffset(1.2);
      charges_signal->GetYaxis()->SetLabelFont(42);
      charges_signal->GetYaxis()->SetTitleFont(42);
      c1->Modified();
	  
      cout << "Electronic Noise Width is " << p22 << "pC" << endl; 
      cout << "The Charge Peak is " << spe_mean << "pC" << endl; 
      cout << "The Charge FWHM is " << spe_sigma*2*sqrt(2*log(2)) << "pC" << endl; 
      cout << "Peak " << spe_mean << ", Valley " << min_function << endl; 
      cout << "The Peak-to-valley is " << Peak_to_valley << endl; 
      cout << "High Charge Tail " << high_charge_entry * 100 / low_charge_entry << "%" << endl; 
      
      //TIMING
      
      // get time difference
      for(j=0; j < number_files; j++){
	for(i=0; i < number_entries; i++){
	  for(k=0; k < number_windows; k++){
	    if(time_max_trig2[j][i][k] != 0){
	      delta_T[j][i][k] = time_max_trig2[j][i][k] - time_max_trig[j][i];	   
	      //chargefile << "Time ETL 20 percent pulse " << time_max_trig2[j][i][k] << ", Time Trigger Pulse " << time_max_trig[j][i] << ", Delta T " << delta_T[j][i][k] << endl;
	      Timing->Fill(delta_T[j][i][k]); 
	    }	
	  }
	}
      }
     
      float time_bin_max = Timing->GetMaximumBin(); 
      TAxis *time_axis = Timing->GetXaxis();
      Double_t time_max = time_axis->GetBinCenter(time_bin_max);

      // draw the fits 
      TCanvas *c2 = new TCanvas("c2","Timing",200,10,700,500);
      Timing->Fit("gaus","0","",time_max - 3.0, time_max + 3.0);
      TF1 *fitfunc=Timing->GetFunction("gaus");
      Double_t f1=fitfunc->GetParameter(1);
      Double_t f2=fitfunc->GetParameter(2); 
      Timing->Draw(); 
      c2->SetLogy(); 
      c2->Update(); 

      
      
      for(j=0; j < number_files; j++){
	for(i=0; i < number_entries; i++){
	  for(k=0; k < number_windows; k++){
	    if(time_max_trig2[j][i][k] != 0){
	      if(delta_T[j][i][k] < f1 + 3*f2 && delta_T[j][i][k] > f1 - 3*f2){ // if delta t is within 3 sigma of prompt peak 
		coincidence_count++;
	      }
	      else if(delta_T[j][i][k] > f1 + 5*f2 && delta_T[j][i][k] < f1 + length_late){ // 5 sigma after prompt peak to length late
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
      
      cout << "Coincidence Hits " << coincidence_count << endl;
      cout << "Late Pulse Hits " << late_pulse_count << endl;
      cout << "Dark Pulse Hits " << dark_count << endl;
      cout << "Pre Pulse Hits " << pre_pulse_count << endl;
      cout << "Post Pulse Hits " << post_pulse_count << endl; 
      cout << "Total hits = " << coincidence_count + late_pulse_count + dark_count + pre_pulse_count + post_pulse_count << endl; 
      
      // make the statistics legend
      TPaveStats *ttstats = new TPaveStats(0.6,0.6,0.98,0.98,"brNDC");
      ttstats->SetName("stats");
      ttstats->SetBorderSize(2);
      ttstats->SetTextAlign(12);
      ttstats->SetTextFont(42);
      ttstats->SetTextSize(0.035);
      ttstats->SetShadowColor(0);
        
      Double_t delta_t_dark = (f1 - 10) + ((length_trace/2)*0.1 - f1 - length_late); cout << "T dark " << delta_t_dark << endl; 
      Double_t delta_t_coincidence = 6*f2; cout << "T coinc " << delta_t_coincidence  << endl; 
      Double_t delta_t_late = length_late - 5*f2; cout << "T late " << delta_t_late  << endl; 
      Double_t delta_t_pre = 10 - 3*f2;
      Double_t t_not_used = 2*f2;
      cout << "Total time = " << delta_t_dark + delta_t_coincidence + delta_t_late + delta_t_pre + t_not_used << endl; 
      
      Double_t ns = 1.0e-9; 
      Double_t Prompt_sigma = f2; 
      Double_t Prompt_FWHM = 2*f2*sqrt(2*log(2));
      Double_t Dark_Rate = dark_count/(delta_t_dark*ns*window_count);  
      Double_t dark_pulse_count1 = Dark_Rate*delta_t_coincidence*ns*window_count;  //number_entries*number_files;
      Double_t dark_pulse_count2 = Dark_Rate*delta_t_late*ns*window_count; //number_entries*number_files;
      Double_t Coincidence_Percent = (coincidence_count - dark_pulse_count1)/(window_count)*100; //(number_entries*number_files)*100;
      chargefile << "Corrected Coincidence Hits " << coincidence_count - dark_pulse_count1 << endl;
      Double_t Late_Hit_Percent = (late_pulse_count - dark_pulse_count2)/(Timing->GetEntries())*100; 
      
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
      Timing->GetXaxis()->SetTitleFont(42);
      Timing->GetXaxis()->SetTitleColor(1);
      Timing->GetXaxis()->SetLabelFont(42);
      Timing->GetYaxis()->SetTitle("Events");
      Timing->GetYaxis()->SetLabelFont(42);
      Timing->GetYaxis()->SetTitleFont(42);
      c2->Modified();
      
      cout << "Prompt Mean " << f1 << endl; 
      cout << "Prompt Sigma " << f2 << endl; 
      cout << "Coincidence Percent " << Coincidence_Percent << endl; 
      cout << "Late Hit Percent " << Late_Hit_Percent << endl; 
      cout << "Dark rate correction to coincidence number " <<  dark_pulse_count1 << ", number of coincidence events " << coincidence_count << endl; 
      cout << "Dark rate correction to late number " <<  dark_pulse_count2 << ", number of late events " << late_pulse_count << endl; 
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
	charges_signal->Write();
	acceptable_waveform->Write();
	pulse_waveform->Write();
	average_waveform->Write(); 
	  
	//Trig PMT 
	big_noise_pedestals->Write(); 
	big_noise_voltages->Write();
	big_noise_charges->Write(); 
	big_average_waveform->Write(); 

	Timing->Write();
	check_average_waveform->Write();
	
      }
  
      // clean-up 
      delete noise_pedestals;
      delete corrected_noise_pedestals; 
      delete voltages;
      delete noise_voltages; 
      delete variances; 
      delete charges_signal;   
      delete acceptable_waveform; 
      delete pulse_waveform;
      delete average_waveform; 

      delete big_noise_pedestals; 
      delete big_noise_voltages; 
      delete big_noise_charges; 
      delete big_average_waveform; 

      delete Timing; 
      delete check_average_waveform;

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
/*
// SPE distibution 
static double spe_CrystalBall(double *x, double *par){
  double CrystalBall=0.0;
  double t=(x[0]-par[0])/(par[1]);
  if (par[2] < 0.0) {
    t=-t;
  }
  double absAlpha=abs(par[2]);
  if (t >= -absAlpha) {
    CrystalBall=par[4]*exp(-t*t/2);
  }
  else {
    double A=pow(par[3]/absAlpha,par[3])*exp(-absAlpha*absAlpha/2.0);
    double B=par[3]/absAlpha-absAlpha;
    if (pow((B-t),par[3])!=0){
      CrystalBall=par[4]*A/pow((B-t),par[3]);
    }
  }

  return CrystalBall; // return SPE distribution
  }*/
/*
void minuitFunction(int& nDim, double* gout, double& result, double par[], int flg) {
  result = spe_CrystalBall(par[0], par[1]);
}
 

 */


