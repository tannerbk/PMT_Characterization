#include <stdio.h>
#include <iostream>
#include <string.h>
#include <math.h>

#include "H5Cpp.h"

#include "calibration.hh"

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

using namespace std;

const std::string cal0 = "/data/snoplus/home/tannerbk/chess_timing/digitizer_calibration/fast.hdf5";
const std::string lappd0 = "/data/snoplus/home/tannerbk/chess_timing/digitizer_calibration/lappd_0.hdf5";
const std::string lappd1 = "/data/snoplus/home/tannerbk/chess_timing/digitizer_calibration/lappd_1.hdf5";

std::vector<double> Calibration::open_calibration_file(std::string filename,
                                                       std::string freq,
                                                       std::string dig_group){
    /*
    Open and parse the map that takes samples to time.
    */
    H5File file;
    file.openFile(filename.c_str(), H5F_ACC_RDONLY);

    Group digitizer_freq = file.openGroup(freq.c_str());
    Group digitizer_group = digitizer_freq.openGroup(dig_group.c_str());

    DataSet dataset = digitizer_group.openDataSet("cell_delay");

    float rdata[1024];
    dataset.read(rdata, PredType::NATIVE_FLOAT);
    std::vector<double> calibration_data;
    for(int i = 0; i < 1024; i++){
        calibration_data.push_back(rdata[i]);
    }

    return calibration_data;
}

int Calibration::get_cal(int index,
                         std::string gr1,
                         std::string gr2,
                         std::vector<double> &cal_pmt,
                         std::vector<double> &cal_trigger){
    /*
    Grab the digitizer calibration information that maps from
    sample number to time. In other words, each sample is not
    exactly 0.2 ns (for a 5 GHz sampling rate).

    TO-DO: Add sampling rate as argument.
    */
    std::string calibration_file;
    if(index == 0) calibration_file = cal0;
    if(index == 1) calibration_file = lappd0;
    if(index == 2) calibration_file = lappd1;

    cal_pmt = open_calibration_file(calibration_file, "5GHz", gr1);
    cal_trigger = open_calibration_file(calibration_file, "5GHz", gr2);

    return 0;
}

double Calibration::get_time(int window_length,
                             int first_sample,
                             int sample,
                             std::vector<double> calibration){
    /*
    Calculates the associated time of 'sample' using the first
    digitized sample and the length of the waeform. 
    */
    int index = fmod((first_sample+sample), window_length);
    double tdiff = (calibration[index] - calibration[first_sample]);
    double time = fmod((fmod(tdiff, 204.8) + 204.8),204.8);

    return time;
}

double Calibration::get_livetime(std::vector<long long> trigger_times){
    /*
    Calculated the livetime of the dataset using the trigger timestamps
    */
    double dt = 0;
    for(size_t i = 0; i < trigger_times.size(); i++){

        if(i == 0) continue;

        double deltat = ((trigger_times[i] - trigger_times[i-1]) & 0x7FFFFFFF)*20.0/1e9;
        dt+=deltat;
    }

    return dt;
}

