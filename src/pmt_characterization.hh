#include <TTree.h>
#include <TH1D.h>
#include <TFile.h>

#ifndef _SIMPLE_TIMING_H_
#define _SIMPLE_TIMING_H_

struct pmt_data {

    // Pedestals for each channel
    double pedestal;

    double stddev;

    // CHESS and trigger timing
    double time;

    // Charge
    double charge;
    double charge_empty;

    int samples_above_threshold;
    int ncrossings;

    // Peak voltages
    double peak_voltage;
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

    void data_structure(TTree *output);
    void meta_structure(TTree *output);

public:

    void pmt_characterization(char* datafile, // input txt file
                              char* output,   // output root file
                              char* dig,      // digitizer
                              char* gr,       // digitizer group
                              char* ch,       // digitizer channel
                              double pedestal_window); // Length of the pedestal window (in samples)
};
#endif

