#include <vector>
#include "H5Cpp.h"
#include "read_data.hh"

#ifndef _TOOLS_
#define _TOOLS_

typedef struct DataCluster DataCluster;

class Tools{

public:

    double get_charge(int start_window,
                      int end_window,
                      DataCluster *datacluster,
                      double pedestal,
                      double dy,
                      double dx);

    double const_threshold(int sample,
                           int lookback,
                           DataCluster *datacluster,
                           double dy,
                           double pedestal,
                           double const_thresh);

    double const_threshold_ttl(int sample,
                               int lookback,
                               DataCluster *datacluster,
                               double dy,
                               double pedestal,
                               double const_thresh);

    double const_frac(int sample,
                      int lookback,
                      DataCluster *datacluster,
                      double dy,
                      double pedestal,
                      double peak_voltage,
                      double threshold_fraction);

    double const_frac_ttl(int sample,
                          int lookback,
                          DataCluster *datacluster,
                          double dy,
                          double pedestal,
                          double peak_voltage,
                          double threshold_fraction);

    double interpolate_no_cal_const(int sample,
                                    DataCluster* datacluster,
                                    double pedestal,
                                    double dy,
                                    double ns_sample,
                                    double const_value);

    double interpolate_no_cal(int sample,
                              DataCluster* datacluster,
                              double pedestal,
                              double dy,
                              double ns_samples,
                              double peak_voltage,
                              double threshold_fraction);

    double interpolate_const(int sample,
                             DataCluster* datacluster,
                             double pedestal,
                             double dy,
                             double peak_voltage,
                             std::vector<double> cal,
                             int index,
                             double const_value);

    double interpolate(int sample,
                       DataCluster* datacluster,
                       double pedestal,
                       double dy,
                       double peak_voltage,
                       std::vector<double> cal,
                       int index,
                       double threshold_fraction);

    double get_voltage(int sample,
                       DataCluster* datacluster,
                       double pedestal,
                       double dy);

    std::vector<double> calculate_pedestals(std::vector<DataCluster*> datacluster,
                                            int pedestal_window_low,
                                            int pedestal_window_high,
                                            double dy);

    double calculate_stddev(DataCluster* datacluster,
                            int pedestal_window_low,
                            int pedestal_window_high,
                            double pedestal,
                            double dy);

};
#endif

