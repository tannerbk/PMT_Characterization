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
  //const unsigned int dark_count = 3; 
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

      //PMT straight on from LED is PMT2
      //PMT at an angle from LED is PMT1

      //Define ROOT histograms here 
      TH1F *voltages = new TH1F("Voltages over Pedestal Window for PMT1","",5000,-1.0,1.0);
      TH1F *noise_voltages = new TH1F("Voltages over Noise Window for PMT1","",5000,-1.0,1.0);
      TH1F *window_voltages = new TH1F("Voltages over Entire Window for PMT1","",5000,-1.0,1.0);
      TH1F *variances = new TH1F("Variance over Pedestal Window","",500000,-0.1,1.0);
      TH1F *noise_variances = new TH1F("Variance over Signal Window","",500000,-0.1,1.0);
      TH1F *noise_pedestals = new TH1F("Pedestal for PMT1","",100000,-1,1); 
      TH1F *corrected_noise_pedestals = new TH1F("Corrected Pedestal for PMT1","",100000,-1,1);
      TH1F *charges_noise = new TH1F("Signal Charge for PMT","",20000,-150,200);  
      TH1F *total_charges = new TH1F("Charge over Entire Window for PMT1","",5000,-150,200); 
      TH1F *acceptable_waveform = new TH1F("Acceptable Waveform","", 1252, 0, 252); 
      TH1F *pulse_waveform = new TH1F("Unacceptable Waveform","", 1252, 0, 252); 
      TH1F *average_waveform = new TH1F("Average Waveform PMT1","", 1252, 0, 252); 
      TH1F *waveform_hist = new TH1F("Avg Waveform Histogram","",10000,-0.05,0.05);     

      //Dark rate stuff commented out 
      //TH1F *dark_voltages = new TH1F("Voltages over Noise Window for Dark Pulses/Background on PMT1","",5000,-1.0,1.0);
      //TH1F *dark_charges = new TH1F("Dark/Background Charge for PMT1","",5000,-150,200); 
      //TGraphErrors *Dark_Graph = new TGraphErrors;
      
      TH1F *big_noise_pedestals = new TH1F("Pedestal for Trigger PMT","",5000,-1.0,1.0); 
      TH1F *big_noise_voltages =  new TH1F("Voltages over Signal Window for Trigger PMT","",5000,-1.0,1.0);
      TH1F *big_noise_charges =  new TH1F("Signal Charge for Trigger PMT","",5000,-150,200); 
      TH1F *big_average_waveform = new TH1F("Average Waveform for Trigger PMT","", 1252, 0, 252); 
      // PMT 2 data commented out

      unsigned long i; 
      unsigned long j; 
      unsigned int window_count = 0.0;
      unsigned int new_window_count = 0.0; 
      //int dark_file_count = 0.0; 
      string filename; 
      H5File file; 
      DataSet dataset;
      //unsigned int plus = 0.0; 
      //Float_t dark_rate[dark_count]; 
      //Float_t dark_rate_errorbars[dark_count]; 
      //Float_t fake_time[dark_count]; 
      const unsigned long length_trace = 2502; 
      

      //for (plus=0.0; plus<dark_count; plus++){ // this is used for dark rate graph
      //	fake_time[plus] = plus;
      //}

      // seperates file input 
      int g = 0; 
      int p = 0;

      // average waveform  
      float waveform_voltage[length_trace] = {0};
      float big_waveform_voltage[length_trace] = {0}; 

      ofstream chargefile;
      chargefile.open("info.text");
	  
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
	  double total_charge = 0.0; 
	  //double dark_charge = 0.0; 
	  float voltage = 0.0; 
	  float noise_voltage = 0.0; 
	  float window_voltage = 0.0; 
	  //float dark_voltage = 0.0; 
	  float variance = 0.0; 
	  float noise_variance = 0.0;  
	  float noise_pedestal = 0.0; 
	  float corrected_noise_pedestal = 0.0; 
	  //int darkrate = 0.0; 

	  for (j = 0; j < datacluster->n_traces; j++){
	    cout << "Analyzing trace " << j+1 << " out of " 
		 <<datacluster->n_traces<< " total traces " << "for channel 3 ZN0112"<<'\r';
	    cout.flush(); 

	    Read_Trace(datacluster,j); 	    	    
	    noise_pedestal = TMath::Mean (window_width, datacluster->data_out)*dy;
	  
	    // initialize 
	    ncharge = 0.0; 
	    total_charge = 0.0; 
	    //dark_charge = 0.0; 
	    variance = 0.0; 
	    noise_variance = 0.0; 
	    noise_voltage = 0.0; 	   
	    
	    for(i=0; i < window_width; i++){
	      voltage = ((float)datacluster->data_out[i]*dy-noise_pedestal);
	      voltages->Fill(voltage);
	      variance = variance + voltage * voltage; 
	    }
	    
	    if(variance < 0.00022) {
	      corrected_noise_pedestal = TMath::Mean (window_width, datacluster->data_out)*dy; 
	      for(i = window_width; i < window_width + integrate_length ; i++){ // noise window
		noise_voltage = ((float)datacluster->data_out[i]*dy-corrected_noise_pedestal);
		noise_voltages->Fill(noise_voltage); 
		noise_variance = noise_variance + noise_voltage * noise_voltage; 
		ncharge = ncharge+(noise_voltage*((-1000.0*dx*1e9)/termination_ohms));
	      }  
	      for(i=0; i < window_length; i++){ // entire window 
		window_voltage = ((float)datacluster->data_out[i]*dy-corrected_noise_pedestal);  
		window_voltages->Fill(window_voltage);
		total_charge = total_charge + (window_voltage*((-1000.0*dx*1e9)/termination_ohms));
		waveform_voltage[i] = waveform_voltage[i] + window_voltage; 
	      }	    
	      for(i=0; i < window_length; i++){ // draws an acceptable trace 
		acceptable_waveform->SetBinContent(i+1,(float)datacluster->data_out[i]*dy-corrected_noise_pedestal); 
	      }

	      window_count++; // number of traces 
	  
	      charges_noise->Fill(ncharge); 
	      total_charges->Fill(total_charge); 
	      noise_variances->Fill(noise_variance); 
	      corrected_noise_pedestals->Fill(corrected_noise_pedestal); 
	    }
	
	    /*
	    if(variance < 0.00022) {
	      for(i = window_width; i < window_width + integrate_length; i++){ // noise window
		dark_voltage = ((float)datacluster->data_out[i]*dy-corrected_noise_pedestal);
		dark_voltages->Fill(dark_voltage);
		dark_charge = dark_charge + (dark_voltage*((-1000.0*dx*1e9)/termination_ohms));
	      }
	      dark_charges->Fill(dark_charge);
	      darkrate++;
	    }
	    */

	    if(variance > 0.0004){
	      for(i=0; i<window_length; i++){
		pulse_waveform->SetBinContent(i+1,(float)datacluster->data_out[i]*dy-noise_pedestal);
	      }
	    }

	    //store histograms 
	    variances->Fill(variance); 
	    noise_pedestals->Fill(noise_pedestal);
	  }

	  /*
	  dark_rate[dark_file_count] = darkrate/0.04; 
	  dark_rate_errorbars[dark_file_count] = sqrt(darkrate)/0.04; 
	  chargefile << "File " << dark_file_count << " has a dark rate of " << dark_rate[dark_file_count] << " with error size " << dark_rate_errorbars[dark_file_count] << endl; 
	  dark_file_count++; 
	  */

	  file.close(); 
	}
	ifs.close(); 
      }

      /*
      // dark rate graph
      TCanvas *c1 = new TCanvas("c1","Charge Cut Dark Rate",200,10,700,500);
      Dark_Graph = new TGraphErrors(dark_count, fake_time, dark_rate, 0 , dark_rate_errorbars);
      Dark_Graph->SetTitle("Charge Cut Dark Rate");
      Dark_Graph->SetMarkerStyle(0);
      Dark_Graph->SetMarkerColor(4);
      Dark_Graph->Draw("AP*");
      c1->Update(); 
      */
 
      Float_t avg_waveform_voltage[length_trace];  
      for(i=0; i<length_trace; i++){ 
	avg_waveform_voltage[i] = waveform_voltage[i]/window_count;
	average_waveform->SetBinContent(i+1, avg_waveform_voltage[i]);  
	waveform_hist->Fill(avg_waveform_voltage[i]); 
      }

      cout  << "\n The number of Traces accepted is " << window_count << " for the ZN0112 PMT." << endl;

      for(p=0;p<1;p++){

	ifstream ifs ( argv[1] , ifstream::in ); // Open File List Stream
	ifs >> filename;
	
	while (ifs.good()){ // While Stream Is Open, Analyze New Files.

	  cout<<filename<< " trig pmt" << endl;
	  file.openFile(filename, H5F_ACC_RDONLY); //Open HDF5 File 
	  ifs >> filename;

	  dataset = file.openDataSet("channel1"); // open channel 3
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
	  
	  // for PMT2 
	  double big_noise_charge = 0.0;
	  float big_noise_voltage = 0.0; 
	  float big_noise_pedestal = 0.0; 
	  float big_window_voltage = 0.0; 
	
	  for (j = 0; j < datacluster->n_traces; j++){
	    cout << "Analyzing trace " << j+1 << " out of " 
	    	 <<datacluster->n_traces<< " total traces " << " for channel 1 trigger PMT."<<'\r';
	    cout.flush(); 

	    Read_Trace(datacluster,j); 
	    big_noise_pedestal = TMath::Mean (500, datacluster->data_out)*dy; // pedestal correction
	    
	    big_noise_charge = 0.0; 

	    for(i = 500; i < 500 + integrate_length; i++){ // signal window for trigger pmt 
	      big_noise_voltage = ((float)datacluster->data_out[i]*dy-big_noise_pedestal);
	      big_noise_voltages->Fill(big_noise_voltage); 
	      big_noise_charge = big_noise_charge+(big_noise_voltage*((-1000.0*dx*1e9)/termination_ohms)); 
	    }

	    for(i = 0; i < window_length ; i++){
	      big_window_voltage = ((float)datacluster->data_out[i]*dy-big_noise_pedestal);
	      big_waveform_voltage[i] = big_waveform_voltage[i] + big_window_voltage; 
	    }

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
      chargefile << first_bin << " " << last_bin << endl;      
      while(first_bin < last_bin){
	charge_finder[f] = charges_noise->GetBinContent(first_bin);
	chargefile << charge_finder[f] << endl; 
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

      // noise fit 
      Int_t zero_bin = axis->FindBin(0.0); 
      Double_t entries_zero_bin = charges_noise->GetBinContent(zero_bin); 
      Double_t entrie_checker = entries_zero_bin;
      while(entrie_checker >= entries_zero_bin/2.0){
	zero_bin++;
	entrie_checker = charges_noise->GetBinContent(zero_bin); 
      }      
      Double_t charge_half_noise_height = charges_noise->GetBinCenter(zero_bin); 
 
      // draw the fits 
      TCanvas *c1 = new TCanvas("c1","Charge",200,10,700,500);
      charges_noise->Fit("gaus","","",bottom_charge_fit, top_charge_fit);
      TF1 *myfunc=charges_noise->GetFunction("gaus");
      Double_t p0=myfunc->GetParameter(0); cout << "Overall normalization for SPE Gaussian " << p0 << endl; 
      Double_t p1=myfunc->GetParameter(1); 
      Double_t p2=myfunc->GetParameter(2); 
      charges_noise->Fit("gaus","","",-charge_half_noise_height, charge_half_noise_height);
      TF1 *myfunc2=charges_noise->GetFunction("gaus");
      //Double_t p00=myfunc2->GetParameter(0);
      //Double_t p11=myfunc2->GetParameter(1); 
      Double_t p22=myfunc2->GetParameter(2);
      charges_noise->Fit("pol2","","",6*p22, bottom_charge_fit);
      TF1 *myfunc3=charges_noise->GetFunction("pol2");
      Double_t p000=myfunc3->GetParameter(0); 
      Double_t p111=myfunc3->GetParameter(1); 
      Double_t p222=myfunc3->GetParameter(2);
      charges_noise->Draw(); 
      c1->Update(); 

      Double_t Low_charge_counter = axis->FindBin(p22); // count number of entries above the electronic noise width
      Double_t High_charge_counter = axis->FindBin(p1 + 3*p2); // count number of entries 3 sigma above the charge peak 
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
	  
      cout << "Electronic Noise Width is " << p22 << endl; 
      cout << "The Charge Peak is " << p1 << endl; 
      cout << "The Charge FWHM is " << p2*2*sqrt(2*log(2)) << endl; 
      cout << "High Charge Tail " << high_charge_entry * 100 / low_charge_entry << "%" << endl; 
      
      chargefile.close(); 

      // Output Histograms to File
      if (argc > 2){
	TFile f(argv[2],"new");
	//PMT1
	noise_pedestals->Write();
	corrected_noise_pedestals->Write();
	voltages->Write(); 
	noise_voltages->Write(); 
	window_voltages->Write();
	//dark_voltages->Write(); 
	variances->Write();
	noise_variances->Write();
	charges_noise->Write();  
	total_charges->Write(); 
	//dark_charges->Write(); 
	acceptable_waveform->Write();
	pulse_waveform->Write();
	average_waveform->Write(); 
	waveform_hist->Write();
	//Dark_Graph->Write();
	  
	//PMT2
	big_noise_pedestals->Write(); 
	big_noise_voltages->Write();
	big_noise_charges->Write(); 
	big_average_waveform->Write(); 
	
      }
  
      // clean-up 
      delete noise_pedestals;
      delete corrected_noise_pedestals; 
      delete voltages;
      delete noise_voltages; 
      delete window_voltages;
      //delete dark_voltages; 
      delete variances; 
      delete noise_variances; 
      delete charges_noise;
      delete total_charges;
      //delete dark_charges;  
      delete acceptable_waveform; 
      delete pulse_waveform;
      delete average_waveform; 
      delete waveform_hist; 
      //delete Dark_Graph; 
      
      delete big_noise_pedestals; 
      delete big_noise_voltages; 
      delete big_noise_charges; 
      delete big_average_waveform; 

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
