#include "tools.hh"
#include "read_data.hh"
#include "const.hh"

#include <cmath>

#include <iostream>
using namespace std;

double Tools::get_charge(int start_window,
                         int end_window,
                         DataCluster *datacluster,
                         double pedestal,
                         double dy,
                         double dx){
    /*
    Integrate the waveform between (start - end window).
    Convert to charge in pC.
    */
    double charge = 0.0;
    for(int i = start_window; i < end_window; i++){
        double voltage = get_voltage(i, datacluster, pedestal, dy);
        charge+=(-voltage*dx)/termination_ohms; // in pC
    }
    return charge;
}

double Tools::const_threshold(int sample,
                              int back,
                              DataCluster *datacluster,
                              double dy,
                              double pedestal,
                              double const_thresh){
    /*
    Given the sample associated with the peak (minimum)
    of the waveform, step back and find the const. thresh
    crossing. Return the associated sample.
    */
    // Start at the peak and work backwards
    for(int i = sample; i > sample - back; i--){
        double voltage = get_voltage(i, datacluster, pedestal, dy);
        if(voltage > const_thresh){
            // Return invalid if the first sample is already above threshold
            if(i == (sample - back)) return INVALID;
            // Otherwise we found the CFD crossing
            return i + 1;
        }
    }

    return INVALID; // no value identified
}

double Tools::const_threshold_ttl(int sample,
                                  int back,
                                  DataCluster *datacluster,
                                  double dy,
                                  double pedestal,
                                  double const_thresh){
    /*
    Given the sample associated with the peak (minimum)
    of the waveform, step back and find the const. thresh
    crossing. Return the associated sample.
    */
    // Start at the peak and work backwards
    for(int i = sample; i > sample - back; i--){
        double voltage = get_voltage(i, datacluster, pedestal, dy);
        if(voltage < const_thresh){
            // Return invalid if the first sample is already above threshold
            if(i == (sample - back)) return INVALID;
            // Otherwise we found the CFD crossing
            return i + 1;
        }
    }

    return INVALID; // no value identified
}


double Tools::const_frac(int sample,
                         int back,
                         DataCluster *datacluster,
                         double dy,
                         double pedestal,
                         double peak_voltage,
                         double threshold_fraction){
    /*
    Given the sample associated with the peak (minimum)
    of the waveform, step back and find the 60% fraction
    crossing. Return the associated sample.
    */
    // Start at the peak and work backwards
    for(int i = sample; i > sample - back; i--){
        double voltage = get_voltage(i, datacluster, pedestal, dy);
        if(voltage > threshold_fraction * peak_voltage){
            // Return invalid if the first sample is already above threshold
            if(i == (sample - back)) return INVALID;
            // Otherwise we found the CFD crossing
            return i + 1;
        }
    }

    return INVALID; // no value identified
}

double Tools::const_frac_ttl(int sample,
                             int back,
                             DataCluster *datacluster,
                             double dy,
                             double pedestal,
                             double peak_voltage,
                             double threshold_fraction){
    /*
    Given the sample associated with the peak (minimum)
    of the waveform, step back and find the 60% fraction
    crossing. Return the associated sample.
    */
    for(int i = sample - back; i < sample; i++){
        double voltage = get_voltage(i, datacluster, pedestal, dy);
        if(voltage > threshold_fraction * peak_voltage){
            // Return invalid if the first sample is already above threshold
            if(i == (sample - back)) return INVALID;
            // Otherwise we found the CFD crossing
            return i;
        }
    }

    return INVALID; // no value identified
}

double Tools::interpolate_no_cal_const(int sample,
                                       DataCluster* datacluster,
                                       double pedestal,
                                       double dy,
                                       double ns_sample,
                                       double const_value){
    /*
    Interpolate between two samples, no calibration applied
    */
    double voltage = get_voltage(sample, datacluster, pedestal, dy);
    double prev_voltage = get_voltage(sample - 1, datacluster, pedestal, dy);
    double deltav = (voltage - prev_voltage);
    double dx = (const_value - prev_voltage)/deltav;
    double dt = ns_sample*dx;

    return dt;
}

double Tools::interpolate_no_cal(int sample,
                                 DataCluster* datacluster,
                                 double pedestal,
                                 double dy,
                                 double ns_sample,
                                 double peak_voltage,
                                 double threshold_fraction){
    /*
    Interpolate between two samples, no calibration applied
    */
    double voltage = get_voltage(sample, datacluster, pedestal, dy);
    double prev_voltage = get_voltage(sample - 1, datacluster, pedestal, dy);
    double deltav = (voltage - prev_voltage);
    double peak_frac = peak_voltage*threshold_fraction;
    double dx = (peak_frac - prev_voltage)/deltav;
    double dt = ns_sample*dx;

    return dt;
}

double Tools::interpolate_const(int sample,
                                DataCluster* datacluster,
                                double pedestal,
                                double dy,
                                double peak_voltage,
                                std::vector<double> cal,
                                int index,
                                double const_thresh){
    /*
    Interpolate between two samples (also convert to time)
    */
    double voltage = get_voltage(sample, datacluster, pedestal, dy);
    double prev_voltage = get_voltage(sample - 1, datacluster, pedestal, dy);
    double deltav = (voltage - prev_voltage);
    double dx = (const_thresh - prev_voltage)/deltav;

    int prev_index = index - 1;
    if(index == 0){
        prev_index = 1023;
    }
    double tdiff = (cal[index] - cal[prev_index]);
    // Calculate time difference with potential rollover
    double time = fmod((fmod(tdiff, 204.8) + 204.8),204.8);
    double dt = time*dx;

    return dt;
}

double Tools::interpolate(int sample,
                          DataCluster* datacluster,
                          double pedestal,
                          double dy,
                          double peak_voltage,
                          std::vector<double> cal,
                          int index,
                          double threshold_fraction){
    /*
    Interpolate between two samples (also convert to time)
    */
    double voltage = get_voltage(sample, datacluster, pedestal, dy);
    double prev_voltage = get_voltage(sample - 1, datacluster, pedestal, dy);
    double deltav = (voltage - prev_voltage);
    double peak_frac = peak_voltage*threshold_fraction;
    double dx = (peak_frac - prev_voltage)/deltav;

    int prev_index = index - 1;
    if(index == 0){
        prev_index = 1023;
    }
    double tdiff = (cal[index] - cal[prev_index]);
    // Calculate time difference with potential rollover
    double time = fmod((fmod(tdiff, 204.8) + 204.8),204.8);
    double dt = time*dx;

    return dt;
}

double Tools::get_voltage(int sample,
                          DataCluster* datacluster,
                          double pedestal,
                          double dy){
    /*
    Convert ADC counts to voltage of a sample, in mV
    */
    return (datacluster->data_out[sample]*dy - pedestal)*1000;
}

std::vector<double> Tools::calculate_pedestals(std::vector<DataCluster*> dataclusters,
                                               int pedestal_window_low,
                                               int pedestal_window_high,
                                               double dy){
    /*
    Calculate the baseline 'pedestals' over the selected pedestal window,
    for each channel in the dataset
    */
    int nchannels = dataclusters.size();
    std::vector<double> pedestals(nchannels, 0);

    for(int i = pedestal_window_low; i < pedestal_window_high; i++){
        for(size_t idc = 0; idc < dataclusters.size(); idc++){
            pedestals[idc] += dataclusters[idc]->data_out[i]*dy;
        }
    }

    for(size_t idc = 0; idc < dataclusters.size(); idc++){
        pedestals[idc] /= (pedestal_window_high - pedestal_window_low);
    }

    return pedestals;
}

double Tools::calculate_stddev(DataCluster* datacluster,
                               int pedestal_window_low,
                               int pedestal_window_high,
                               double pedestal,
                               double dy){
    /*
    Calculate the baseline stddev over the selected pedestal window,
    for each channel in the dataset
    */
    double stddev = 0;
    for(int i = pedestal_window_low; i < pedestal_window_high; i++){
        double voltage = (datacluster->data_out[i]*dy);
        stddev += pow((voltage - pedestal), 2);
    }
    stddev *= 1000; // convert to mV
    stddev /= (pedestal_window_high - pedestal_window_low);
    stddev = sqrt(stddev);

    return stddev;
}

