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
#include <TLegend.h>

// HDF5 Library
#include "H5Cpp.h"

#ifndef PMT_CHARACTERIZATION
#define PMT_CHARACTERIZATION

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

using namespace std;


/* Find the start of the window by looking at average waveforms */
int GetWindowStart(char* datafile);
double GetSampling(char* datafile);
size_t GetLength(char* datafile);


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
const string trigger_channel = "channel1";
const string pmt_channel = "channel2";
unsigned long window_width;

/* Plot Settings */
const int high_voltage = 1840;
TString box_title = Form("R5912-MOD");

/* Various modes */
int mode;
const char* converged = "CONVERGED ";

/* Analysis settings */
const int length_late = 45;
const float variance_cut = 0.005;
const int prompt_peak_fraction = 0.2;
// Charge cuts to fill timing histogram
const float charge_cut_low = 0.4;
const float charge_cut_high = 3.0;
// Charge cuts to integrate above to count coincidence
const float charge_cut_integral = 0.4;
const int size_signal_window = 300;
const int termination_ohms = 75; /* SNO PMTs are 75 ohms terminated */
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
vector<double> prompt_q;
vector<double> late_q;

/* 2D charge time vector */
vector<double> kTime;
vector<double> kCharge;

// Histograms for both channels (trigger and analysis channel)
TH1F *pedestals = new TH1F("Pedestal","",100000,-1,1);
TH1F *variances = new TH1F("Variance","",500000,-0.1,1.0);
TH1F *charges_signal = new TH1F("Charge","",800,-4.0,12.0);
TH1F *peaks = new TH1F("Peaks","",5000,-1.0,0.1);
TH1F *trigger_pedestals = new TH1F("Pedestal_Trigger","",100000,-1.0,1.0);

TH1F *prompt_hits_charge = new TH1F("promptq","",800,-4.0,12.0);
TH1F *late_hits_charge = new TH1F("lateq","",800,-4.0,12.0);
TH1F *double_hits_prompt_charge = new TH1F("double_promptq","",800,-4.0,12.0);
TH1F *double_hits_late_charge = new TH1F("double_lateq","",800,-4.0,12.0);
TH1F *dark_hits_charge = new TH1F("darkq","",800,-4.0,12.0);

TH1F *prompt_charge = new TH1F("prompt_charge","",800,-4.0,12.0);

#endif
