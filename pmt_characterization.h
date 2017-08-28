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
#include <TMinuit.h>
#include <TVirtualFitter.h>

// HDF5 Library
#include "H5Cpp.h"

using namespace std;

#ifndef PMT_CHARACTERIZATION
#define PMT_CHARACTERIZATION

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

using namespace std;


/* Stuct for HDF5 information */
const int RANK_OUT = 2;

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


/* Find the start of the window by looking at average waveforms */
int GetWindowStart(char* datafile);
double GetSampling(char* datafile);
size_t GetLength(char* datafile);


/* Lognormal fits */
double Lognormal(double t, double tau, double sigma, double mag);
double SingleLognormal(Double_t* x, Double_t* par);
double DoubleLognormal(Double_t* x, Double_t* par);
double TripleLognormal(Double_t* x, Double_t* par);


/* Datacluser attributes */
Attribute horiz_interval;
Attribute vertical_gain;
Attribute vertical_gain_trigger;


/* Attribute Value Variables and datafiles */
double dx,dy,dy2;
H5File file;
DataSet dataset;
DataSet dataset_trigger;
char* datafile;


/* Scope Settings */
const string trigger_channel = "channel2";
const string pmt_channel = "channel3";
unsigned long window_width;


/* Plot Settings */
const int high_voltage = 1840;
TString box_title = Form("R5912-MOD");


/* Various modes */
int mode;
const char* converged = "CONVERGED ";


/* Analysis settings */
const int length_late = 60;
const float variance_cut = 0.005;
// Charge cuts to fill timing histogram
const float charge_cut_low = 0.4;
const float charge_cut_high = 3.0;
// Charge cuts to integrate above to count coincidence
const float charge_cut_integral = 0.3;
const int size_signal_window = 400;
const int termination_ohms = 50;
string filename;


/* Vectors to store timing info */
vector<double> trigger_time;
vector<double> pmt_time;

/* Average waveform vectors */
vector<double> waveform_voltage;
vector<double> waveform_voltage_trigger;


/* Vectors for double pulsing */
vector<double> prompt_pulse;
vector<double> late_pulse;

#endif
