#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include <TFile.h>
#include <TTree.h>

#include "H5Cpp.h"

#include "pmt_characterization.hh"
#include "calibration.hh"
#include "read_data.hh"
#include "tools.hh"
#include "const.hh"
#include "version.hh"

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

using namespace std;

typedef struct DataCluster DataCluster;
pmt_data pdata;
pmt_meta_data meta;

void help(){

    printf("Run the code as:\n");
    printf("./pmt_characterization input.txt output fast gr ch ped_window\n");
}

int main(int argc, char* argv[]){

     if(argc != 7){
         printf("Incorrect number of arguments\n.");
         help();
         return 0;
     }

     char* datafile = argv[1];
     char* output = argv[2];
     char* dig = argv[3];
     char* gr = argv[4];
     char* ch = argv[5];
     float ped_window = atof(argv[6]);

     PMTChar fPMTChar;

     fPMTChar.pmt_characterization(datafile, output, dig, gr, ch, ped_window);

     return 0;
}

void PMTChar::data_structure(TTree* output){
    /*
    The output data structure
    */
    output->Branch("time", &pdata.time);

    output->Branch("pedestal", &pdata.pedestal);

    output->Branch("stddev", &pdata.stddev);

    output->Branch("charge", &pdata.charge);

    output->Branch("samples_above_threshold", &pdata.samples_above_threshold);
    output->Branch("ncrossings", &pdata.ncrossings);

    output->Branch("peak_voltage", &pdata.peak_voltage);
}

void PMTChar::meta_structure(TTree* output){
    /*
    The meta structure
    */
    output->Branch("nfiles", &meta.nfiles);
    output->Branch("nwaveforms", &meta.nwaveforms);
    output->Branch("pedestal_window", &meta.pedestal_window);
    output->Branch("coincidence_count", &meta.coincidence_count);
    output->Branch("analysis_code_version", &meta.version);
}

void PMTChar::pmt_characterization(char* datafile,
                                   char* outname,
                                   char* dig,
                                   char* gr,
                                   char* ch,
                                   double pedestal_window){
    /*
    Generate the charge and timing information 
    */
    Calibration fCal;
    ReadData fReadData;
    Tools fTools;

    // ns per sample
    Attribute _dx;
    double dx = 0;
    // ADC bits
    Attribute _bits;
    double bits = 0;

    // For average waveform
    std::vector<double> waveform_voltage;

    // Open output filename and write data to a TTree
    char fname[256];
    sprintf(fname, "%s_%s_%s_%s.root", outname, dig, gr, ch);
    TFile *fout = new TFile(fname, "RECREATE");
    TTree *output = new TTree("output", "output");
    TTree *meta_output = new TTree("meta", "meta");
    data_structure(output);
    meta_structure(meta_output);

    TH1D* wfm = new TH1D("wfm","wfm",1024,0,204.8);

    // Load calibration information, maps samples->time
    std::vector<double> cal_pmt;
    std::vector<double> cal_trigger; // not used
    if(fCal.get_cal(1, gr, gr, cal_pmt, cal_trigger)){
        std::cout << "Invalid calibration file." << std::endl;
        return;
    }

    // Load in list of data files
    string filename;
    ifstream ifs (datafile, ifstream::in);
    ifs >> filename;

    H5File file;
    meta.nfiles = 0;
    meta.nwaveforms = 0;
    meta.coincidence_count = 0;

    while(ifs.good()){

        // Open HDF5 File, then HDF5 DataSet and Read in Attributes
        file.openFile(filename.c_str(), H5F_ACC_RDONLY);
        ifs >> filename;

        Group digitizer = file.openGroup(dig);
        Group digitizer_group = digitizer.openGroup(gr);
        Group channel = digitizer_group.openGroup(ch);
        DataSet dataset = channel.openDataSet("samples");

        // Initialize the dataclusters
        DataCluster *datacluster = fReadData.Init_Data(&dataset);

        std::vector<DataCluster*> dataclusters;
        dataclusters.push_back(datacluster);

        // Specific to digitizer, so in principle could be different for trigger PMT
        _dx = digitizer.openAttribute("ns_sample");
        _bits = digitizer.openAttribute("bits");
        _dx.read(PredType::NATIVE_DOUBLE, &dx);
        _bits.read(PredType::NATIVE_DOUBLE, &bits);

        // Waveform length and resolution
        unsigned long window_length = datacluster->trace_length;
        double dy = DYNAMIC_RANGE/pow(2, bits);

        if(meta.nfiles == 0){
            cout << " * * * * * * * * * * * * * * * * * * * " << endl;
            cout << "Time/sample:      " << dx << " ns" << endl;
            cout << "Resolution:       " << dy << " V" << endl;
            cout << "Trace length:     " << window_length << " samples or " << window_length*dx << " ns." << endl;
            cout << "Pedestal window:  " << sample_offset*dx << " - " << pedestal_window*dx << " ns "
                 << "or " << sample_offset << " - " << pedestal_window << " samples." << endl;
            cout << "Scanning:         " << pedestal_window*dx << " - " << window_length - sample_offset_last 
                 << " ns or " << pedestal_window << " - " << window_length - sample_offset << " samples." << endl;
            cout << " * * * * * * * * * * * * * * * * * * * " << endl;
            cout << "PMT threshold: " << voltage_threshold << " mV" << endl;
            cout << " * * * * * * * * * * * * * * * * * * * " << endl;
        }
        else if(meta.nfiles % 200 == 0){
            cout << "File: " << filename << " (" << meta.nfiles << ")" << endl;
        }

        // Total number of digitized waveforms. Should be the same for all channels.
        const int ntraces = datacluster->n_traces;
        // First indexed that is digitized
        std::vector<int> start_index(ntraces, 0);
        DataSet dataset_start_index = digitizer_group.openDataSet("start_index");
        fReadData.Read_Trace_1D(dataset_start_index, ntraces, start_index);

        for(int j = 0; j < ntraces; j++){

            meta.nwaveforms += 1;

            if(j == 0) waveform_voltage.resize(window_length);

            // Read the waveforms for the CHESS and trigger PMT and the associated
            // raw triggers for the corresponding digitizer.
            for(size_t idc = 0; idc < dataclusters.size(); idc++){
                fReadData.Read_Trace(dataclusters[idc], j);
            }

            // Baseline over the selected pedestal window
            std::vector<double> pedestals = fTools.calculate_pedestals(dataclusters,
                                                                       sample_offset,
                                                                       pedestal_window,
                                                                       dy);

            // Baseline calculation for each channel
            pdata.pedestal = pedestals[0];

            pdata.stddev = fTools.calculate_stddev(dataclusters[0],
                                                   sample_offset,
                                                   pedestal_window,
                                                   pdata.pedestal, dy);

            // Location of the peak in samples
            peak_bin = 0;

            // Peak (in mV) for each waveform
            pdata.peak_voltage = 0;

            // N samples above threshold for CHESS PMT
            pdata.samples_above_threshold = 0;
            double nsamples_above = 0;

            bool crossed = false;
            pdata.ncrossings = 0;

            // Loop through window after pedestal looking for pulses that crossed threshold.
            for(unsigned int i = pedestal_window; i < (window_length - sample_offset_last); i++){

                double voltage = fTools.get_voltage(i, datacluster, pdata.pedestal, dy);

                waveform_voltage.at(i) += voltage;

                if(simple_write_waveforms){
                    wfm->SetBinContent(i, voltage);
                }

                // Peak voltages (minimum) for CHESS PMT
                if(voltage < pdata.peak_voltage){
                    pdata.peak_voltage = voltage;
                    peak_bin = i;
                }

                // Count the number of threshold crossings
                if(voltage < voltage_threshold){
                    nsamples_above += 1;
                    if(!crossed){
                        pdata.ncrossings+=1;
                    }
                    crossed = true;
                }
                else{
                    nsamples_above = 0;
                    crossed = false;
                }
                if(nsamples_above > pdata.samples_above_threshold){
                    pdata.samples_above_threshold = nsamples_above;
                }
            }

            // Defines a pulse that crosses threshold
            if(pdata.peak_voltage >= voltage_threshold) continue;

            // Constant fraction discriminator applied to CHESS PMT waveform
            int sample = fTools.const_frac(peak_bin, lookback,
                                           datacluster, dy,
                                           pdata.pedestal,
                                           pdata.peak_voltage,
                                           const_frac_thresh);

            // Index of the start index + threshold crossing sample
            int index = fmod((start_index[j]+sample),window_length);
            // Apply a linear interpolation between two samples around threshold crossing
            double dt = fTools.interpolate(sample, datacluster,
                                           pdata.pedestal, dy,
                                           pdata.peak_voltage, cal_pmt,
                                           index, const_frac_thresh);
            // Convert from sample to time using calibrated information
            pdata.time = fCal.get_time(window_length, start_index[j],
                                       sample-1, cal_pmt);
            pdata.time += dt;

            pdata.charge = fTools.get_charge(peak_bin-integration_samples_back,
                                             peak_bin+integration_samples_forward,
                                             datacluster, pdata.pedestal, dy, dx);

            meta.coincidence_count += 1;
            output->Fill();
        }
        file.close();
        meta.nfiles += 1;
    }
    ifs.close();

    cout << "Total coincidence rate: " << double(meta.coincidence_count)*100/meta.nwaveforms << "%" << endl;
    meta.pedestal_window = pedestal_window;
    meta.version = VERSION;
    meta_output->Fill();
    meta_output->Write();
    output->Write();
    fout->Close();
}

