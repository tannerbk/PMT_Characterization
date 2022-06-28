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
pmt_data data;
pmt_meta_data meta;

void help(){

    printf("Run the code as:\n");
    printf("./pmt_characterization input.txt output fast gr ch gr_trig ch_trig gr_empty ch_empty ped_window\n");
}

int main(int argc, char* argv[]){

     if(argc != 12){
         printf("Incorrect number of arguments\n.");
         help();
         return 0;
     }

     char* datafile = argv[1];
     char* output = argv[2];
     char* dig = argv[3];
     char* gr = argv[4];
     char* ch = argv[5];
     char* gr_trig = argv[6];
     char* ch_trig = argv[7];
     char* gr_empty = argv[8];
     char* ch_empty = argv[9];
     float ped_window = atof(argv[10]);
     int led = atoi(argv[11]);

     PMTChar fPMTChar;

     fPMTChar.pmt_characterization(datafile, output, dig, gr, ch, gr_trig, ch_trig, gr_empty, ch_empty, ped_window, led);

     return 0;
}

void PMTChar::data_structure(TTree* output){
    /*
    The output data structure
    */
    output->Branch("time", &data.time);
    output->Branch("time_ttl", &data.time_ttl);
    output->Branch("dt_pmt", &data.dt);
    output->Branch("time_trigger", &data.time_trigger);
    output->Branch("time_trigger_ttl", &data.time_trigger_ttl);
    output->Branch("dt_trigger", &data.dt_trigger);
    output->Branch("deltat", &data.deltat);

    output->Branch("pedestal", &data.pedestal);
    output->Branch("pedestal_trigger", &data.pedestal_trigger);
    output->Branch("pedestal_tr", &data.pedestal_tr);
    output->Branch("pedestal_trigger_tr", &data.pedestal_trigger_tr);
    output->Branch("pedestal_empty", &data.pedestal_empty);

    output->Branch("stddev", &data.stddev);
    output->Branch("stddev_trigger", &data.stddev_trigger);

    output->Branch("charge", &data.charge);
    output->Branch("charge_empty", &data.charge_empty);
    output->Branch("trigger_charge", &data.trigger_charge);
    output->Branch("trigger_charge_empty", &data.trigger_charge_empty);

    output->Branch("samples_above_threshold", &data.samples_above_threshold);
    output->Branch("ncrossings", &data.ncrossings);

    output->Branch("peak_voltage", &data.peak_voltage);
    output->Branch("peak_voltage_trigger", &data.peak_voltage_trigger);
    output->Branch("peak_tr_voltage", &data.peak_tr_voltage);
    output->Branch("peak_tr_voltage_trigger", &data.peak_tr_voltage_trigger);
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
                                   char* gr_trig,
                                   char* ch_trig,
                                   char* gr_empty,
                                   char* ch_empty,
                                   double pedestal_window,
                                   int led){
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

    // Open output filename and write data to a TTree
    char fname[256];
    sprintf(fname, "%s_%s_%s_%s.root", outname, dig, gr, ch);
    TFile *fout = new TFile(fname, "RECREATE");
    TTree *output = new TTree("output", "output");
    TTree *meta_output = new TTree("meta", "meta");
    data_structure(output);
    meta_structure(meta_output);

    TH1D* wfm = new TH1D("wfm","wfm",1024,0,204.8);
    TH1D* wfm_trig = new TH1D("wfm_trig","wfm_trig",1024,0,204.8);

    // Load calibration information, maps samples->time
    std::vector<double> cal_pmt;
    std::vector<double> cal_trigger;
    if(fCal.get_cal(1, gr, gr_trig, cal_pmt, cal_trigger)){
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

        // Grab the CHESS PMT information
        Group digitizer = file.openGroup(dig);
        Group digitizer_group = digitizer.openGroup(gr);

        Group channel = digitizer_group.openGroup(ch);
        Group channel_tr = digitizer_group.openGroup("tr");

        DataSet dataset = channel.openDataSet("samples");
        DataSet dataset_tr = channel_tr.openDataSet("samples");

        // Grab the trigger PMT information
        Group trigger_digitizer_group = digitizer.openGroup(gr_trig);
        Group trigger_channel = trigger_digitizer_group.openGroup(ch_trig);
        Group trigger_channel_tr = trigger_digitizer_group.openGroup("tr");

        DataSet trigger_dataset = trigger_channel.openDataSet("samples");
        DataSet trigger_dataset_tr = trigger_channel_tr.openDataSet("samples");

        // Grab the empty PMT information
        Group empty_digitizer_group = digitizer.openGroup(gr_empty);
        Group empty_channel = empty_digitizer_group.openGroup(ch_empty);
        DataSet empty_dataset = empty_channel.openDataSet("samples");

        // Initialize the dataclusters
        DataCluster *datacluster = fReadData.Init_Data(&dataset);
        DataCluster *datacluster_trigger = fReadData.Init_Data(&trigger_dataset);

        DataCluster *datacluster_tr = fReadData.Init_Data(&dataset_tr);
        DataCluster *datacluster_trigger_tr = fReadData.Init_Data(&trigger_dataset_tr);

        DataCluster *datacluster_empty = fReadData.Init_Data(&empty_dataset);

        std::vector<DataCluster*> dataclusters;
        dataclusters.push_back(datacluster);
        dataclusters.push_back(datacluster_trigger);
        dataclusters.push_back(datacluster_tr);
        dataclusters.push_back(datacluster_trigger_tr);
        dataclusters.push_back(datacluster_empty);

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
            //cout << "Pedestal window:  " << pedestal_window*dx << " " << (window_length - sample_offset_last)*dx << " ns "
                 << "or " << sample_offset << " - " << pedestal_window << " samples." << endl;
                 //<< "or " << pedestal_window << " - " << (window_length - sample_offset_last) << " samples." << endl;
            cout << "Scanning:         " << pedestal_window*dx << " - " << window_length - sample_offset_last 
                 << " ns or " << pedestal_window << " - " << window_length - sample_offset << " samples." << endl;
            cout << " * * * * * * * * * * * * * * * * * * * " << endl;
            cout << "Using constant fraction for trigger: " << CONST_FRAC_TRIGGER << endl; 
            cout << "Trigger threshold: " << trigger_voltage_threshold << " mV" << endl;
            cout << "Array PMT threshold: " << voltage_threshold << " mV" << endl;
            cout << " * * * * * * * * * * * * * * * * * * * " << endl;
        }
        else if(meta.nfiles % 200 == 0){
            cout << "File: " << filename << " (" << meta.nfiles << ")" << endl;
        }

        // Total number of digitized waveforms. Should be the same for all channels.
        const int ntraces = datacluster->n_traces;
        // First indexed that is digitized
        std::vector<int> start_index(ntraces, 0);
        std::vector<int> start_index_trigger(ntraces, 0);
        DataSet dataset_start_index = digitizer_group.openDataSet("start_index");
        DataSet dataset_trigger_start_index = trigger_digitizer_group.openDataSet("start_index");
        fReadData.Read_Trace_1D(dataset_start_index, ntraces, start_index);
        fReadData.Read_Trace_1D(dataset_trigger_start_index, ntraces, start_index_trigger);

        std::vector<long long> trigger_times(ntraces, 0);
        DataSet dataset_trigger_time = digitizer_group.openDataSet("trigger_time");
        fReadData.Read_Trace_1D_Long(dataset_trigger_time, ntraces, trigger_times);

        for(int j = 0; j < ntraces; j++){

            meta.nwaveforms += 1;

            // Read the waveforms for the CHESS and trigger PMT and the associated
            // raw triggers for the corresponding digitizer.
            for(size_t idc = 0; idc < dataclusters.size(); idc++){
                fReadData.Read_Trace(dataclusters[idc], j);
            }

            // Baseline over the selected pedestal window
            std::vector<double> pedestals = fTools.calculate_pedestals(dataclusters,
                                                                       sample_offset,
                                                                       pedestal_window,
                                                                       //window_length - sample_offset_last,
                                                                       dy);

            // Baseline calculation for each channel
            data.pedestal = pedestals[0];
            data.pedestal_trigger = pedestals[1];
            data.pedestal_tr = pedestals[2];
            data.pedestal_trigger_tr = pedestals[3];
            data.pedestal_empty = pedestals[4];

            data.stddev = fTools.calculate_stddev(dataclusters[0],
                                                  sample_offset,
                                                  pedestal_window,
                                                  //window_length - sample_offset_last,
                                                  data.pedestal, dy);

            data.stddev_trigger = fTools.calculate_stddev(dataclusters[1],
                                                          sample_offset,
                                                          pedestal_window,
                                                          //window_length - sample_offset_last,
                                                          data.pedestal_trigger, dy);

            // Location of the peak in samples
            peak_bin = 0;
            peak_bin_trigger = 0;
            peak_tr_bin = 0;
            peak_tr_bin_trigger = 0;

            // Peak (in mV) for each waveform
            data.peak_voltage = 0;
            data.peak_tr_voltage = 0;
            data.peak_voltage_trigger = 0;
            data.peak_tr_voltage_trigger = 0;

            // N samples above threshold for CHESS PMT
            data.samples_above_threshold = 0;
            double nsamples_above = 0;

            bool crossed = false;
            data.ncrossings = 0;

            // Loop through full window looking for pulses that crossed threshold.
            //for(unsigned int i = sample_offset; i < (window_length - sample_offset_last); i++){
            // Loop through window after pedestal looking for pulses that crossed threshold.
            for(unsigned int i = pedestal_window; i < (window_length - sample_offset_last); i++){
            //for(unsigned int i = sample_offset; i < pedestal_window; i++){

                double voltage = fTools.get_voltage(i, datacluster, data.pedestal, dy);
                double trigger_voltage = fTools.get_voltage(i, datacluster_trigger, data.pedestal_trigger, dy);
                double tr_voltage = fTools.get_voltage(i, datacluster_tr, data.pedestal_tr, dy);
                double tr_trigger_voltage = fTools.get_voltage(i, datacluster_trigger_tr, data.pedestal_trigger_tr, dy);

                if(simple_write_waveforms){
                    wfm->SetBinContent(i, voltage);
                    //wfm->SetBinError(i, data.stddev);
                    wfm_trig->SetBinContent(i, trigger_voltage);
                    //wfm_trig->SetBinError(i, data.stddev_trigger);
                }

                // Peak voltage for upward going TTL pulse, CHESS PMT
                if(tr_voltage > data.peak_tr_voltage){
                    data.peak_tr_voltage = tr_voltage;
                    peak_tr_bin = i;
                }

                // Peak voltage for upward going TTL pulse, trigger PMT
                if(tr_trigger_voltage > data.peak_tr_voltage_trigger){
                    data.peak_tr_voltage_trigger = tr_trigger_voltage;
                    peak_tr_bin_trigger = i;
                }

                // Peak voltages (minimum) for CHESS PMT
                if(voltage < data.peak_voltage){
                    data.peak_voltage = voltage;
                    peak_bin = i;
                }

                // Count the number of threshold crossings
                if(voltage < voltage_threshold){
                    nsamples_above += 1;
                    if(!crossed){
                        data.ncrossings+=1;
                    }
                    crossed = true;
                }
                else{
                    nsamples_above = 0;
                    crossed = false;
                }
                if(nsamples_above > data.samples_above_threshold){
                    data.samples_above_threshold = nsamples_above;
                }

                // Peak voltage for trigger PMT
                if(trigger_voltage < data.peak_voltage_trigger){
                    data.peak_voltage_trigger = trigger_voltage;
                    peak_bin_trigger = i;
                }
            }

            if(!led){

              // Defines a pulse that crosses threshold
              //if(data.peak_voltage >= voltage_threshold) continue;
              if(data.peak_voltage_trigger >= trigger_voltage_threshold) continue;

              int trigger_sample = 0;
              if(CONST_FRAC_TRIGGER){
                  // Constant fraction discriminator applied to trigger PMT
                  trigger_sample = fTools.const_frac(peak_bin_trigger, lookback_trigger,
                                                     datacluster_trigger, dy,
                                                     data.pedestal_trigger, data.peak_voltage_trigger,
                                                     const_frac_thresh_trigger);
              }
              else{
                  // Constant threshold discriminator applied to trigger PMT
                  trigger_sample = fTools.const_threshold(peak_bin_trigger, lookback_trigger,
                                                          datacluster_trigger, dy,
                                                          data.pedestal_trigger,
                                                          const_thresh_trigger);
              }

              // Constant fraction discriminator applied to CHESS PMT waveform
              int sample = fTools.const_frac(peak_bin, lookback,
                                             datacluster, dy,
                                             data.pedestal,
                                             data.peak_voltage,
                                             const_frac_thresh);

              // Constant fraction discriminiator applied to TTL pulse for CHESS PMT
              int tr_sample = fTools.const_frac_ttl(peak_tr_bin, lookback_trigger,
                                                    datacluster_tr, dy,
                                                    data.pedestal_tr,
                                                    data.peak_tr_voltage,
                                                    const_frac_thresh);

              // Constant fraction discriminiator applied to TTL pulse for trigger PMT
              int tr_trigger_sample = fTools.const_frac_ttl(peak_tr_bin_trigger, lookback_trigger,
                                                            datacluster_trigger_tr, dy,
                                                            data.pedestal_trigger_tr,
                                                            data.peak_tr_voltage_trigger,
                                                            const_frac_thresh);

              // Index of the start index + threshold crossing sample
              int index = fmod((start_index[j]+sample),window_length);
              // Apply a linear interpolation between two samples around threshold crossing
              double dt = fTools.interpolate(sample, datacluster,
                                             data.pedestal, dy,
                                             data.peak_voltage, cal_pmt,
                                             index, const_frac_thresh);
              // Convert from sample to time using calibrated information
              data.time = fCal.get_time(window_length, start_index[j],
                                        sample-1, cal_pmt);
              data.time += dt;

              // Repeat for the CHESS PMT group trigger signal
              int tr_index = fmod((start_index[j]+tr_sample),window_length);
              double dt_tr = fTools.interpolate(tr_sample, datacluster_tr,
                                                data.pedestal_tr, dy,
                                                data.peak_tr_voltage, cal_pmt,
                                                tr_index, const_frac_thresh);
              data.time_ttl = fCal.get_time(window_length, start_index[j],
                                            tr_sample-1, cal_pmt);
              data.time_ttl += dt_tr;

              // Time difference between CHESS PMT time and TTL pulse
              data.dt = data.time - data.time_ttl;

              // Index of the start index + threshold crossing sample (trigger PMT)
              int trigger_index = fmod((start_index_trigger[j]+trigger_sample),window_length);

              double dt_trig = 0;
              // Apply a linear interpolation between two samples around threshold crossing
              if(CONST_FRAC_TRIGGER){
                  dt_trig = fTools.interpolate(trigger_sample, datacluster_trigger,
                                               data.pedestal_trigger, dy,
                                               data.peak_voltage_trigger, cal_trigger,
                                               trigger_index, const_frac_thresh_trigger);
              }
              else{
                  dt_trig = fTools.interpolate_const(trigger_sample, datacluster_trigger,
                                                     data.pedestal_trigger, dy,
                                                     data.peak_voltage_trigger, cal_trigger,
                                                     trigger_index, const_thresh_trigger);
              }

              // Convert from sample to time using calibrated information
              data.time_trigger = fCal.get_time(window_length, start_index_trigger[j],
                                                trigger_sample-1, cal_trigger);
              data.time_trigger += dt_trig;

              // Repeat for the trigger PMT group trigger signal
              int tr_trigger_index = fmod((start_index_trigger[j]+tr_trigger_sample),window_length);
              double dt_tr_trig = fTools.interpolate(tr_trigger_sample, datacluster_trigger_tr,
                                                     data.pedestal_trigger_tr, dy,
                                                     data.peak_tr_voltage_trigger, cal_trigger,
                                                     tr_trigger_index, const_frac_thresh);
              // Convert from sample to time using calibrated information
              data.time_trigger_ttl = fCal.get_time(window_length, start_index_trigger[j],
                                                    tr_trigger_sample-1, cal_trigger);
              data.time_trigger_ttl += dt_tr_trig;

              // Time difference between trigger time and TTL pulse
              data.dt_trigger = data.time_trigger - data.time_trigger_ttl;

              // Ultimately it is the difference in time differences that we use
              data.deltat = (data.dt - data.dt_trigger);

              // Integrate the trigger waveform
              data.trigger_charge = fTools.get_charge(trigger_sample-integration_samples_back,
                                                      trigger_sample+integration_samples_forward,
                                                      datacluster_trigger, data.pedestal_trigger, dy, dx);

              data.trigger_charge_empty = fTools.get_charge(trigger_sample-integration_samples_back,
                                                            trigger_sample+integration_samples_forward,
                                                            datacluster_empty, data.pedestal_empty, dy, dx);

              data.charge = fTools.get_charge(peak_bin-integration_samples_back,
                                              peak_bin+integration_samples_forward,
                                              datacluster, data.pedestal, dy, dx);

              data.charge_empty = fTools.get_charge(peak_bin-integration_samples_back,
                                                    peak_bin+integration_samples_forward,
                                                    datacluster_empty, data.pedestal_empty, dy, dx);

            }

            else{
              data.charge = fTools.get_charge(peak_bin-integration_samples_back_led,
                                              peak_bin+integration_samples_forward_led,
                                              datacluster, data.pedestal, dy, dx);

              data.charge_empty = fTools.get_charge(peak_bin-integration_samples_back_led,
                                                    peak_bin+integration_samples_forward_led,
                                                    datacluster_empty, data.pedestal_empty, dy, dx);
            }

            if(simple_write_waveforms){
                wfm->Write();
                wfm_trig->Write();
            }

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

