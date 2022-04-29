#include <TTree.h>
#include <TH1D.h>
#include <TFile.h>

#ifndef _SIMPLE_TIMING_H_
#define _SIMPLE_TIMING_H_

struct pmt_data {

    // Pedestals for each channel
    double pedestal;
    double pedestal_trigger;
    double pedestal_tr;
    double pedestal_trigger_tr;
    double pedestal_empty;

    double stddev;
    double stddev_trigger;

    // CHESS and trigger timing
    double time;
    double time_ttl;
    double dt;
    double time_trigger;
    double time_trigger_ttl;
    double dt_trigger;
    double deltat;

    // Charge
    double charge;
    double charge_empty;
    double trigger_charge;
    double trigger_charge_empty;

    int samples_above_threshold;
    int ncrossings;

    // Peak voltages
    double peak_voltage;
    double peak_tr_voltage;
    double peak_voltage_trigger;
    double peak_tr_voltage_trigger;
};

struct pmt_meta_data {

    int nfiles;
    int nwaveforms;
    int pedestal_window;
    int coincidence_count;
    double version;
};

class PMTChar {

private:

    // Bin corresponding to peak
    int peak_bin;
    int peak_bin_trigger;
    int peak_tr_bin;
    int peak_tr_bin_trigger;

    void data_structure(TTree *output);
    void meta_structure(TTree *output);

public:

    void pmt_characterization(char* datafile, // input txt file
                              char* output,   // output root file
                              char* dig,      // digitizer
                              char* gr,       // digitizer group
                              char* ch,       // digitizer channel
                              char* gr_trig,  // digitizer group for trigger PMT
                              char* ch_trig,  // digitizer channel for trigger PMT
                              char* gr_empty, // digitizer group for empty channel
                              char* ch_empty, // digitizer channel for empty channel
                              double pedestal_window); // Length of the pedestal window (in samples)

};
#endif

