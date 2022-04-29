#include <string>
#include <vector>

#ifndef _CALIBRATION_H_
#define _CALIBRATION_H_

class Calibration {

private:

    std::vector<double> open_calibration_file(std::string filename,
                                              std::string freq,
                                              std::string dig_group);

public:

    int get_cal(int index,
                std::string gr1,
                std::string gr2,
                std::vector<double> &cal_pmt,
                std::vector<double> &cal_trigger);

    double get_time(int window_length,
                    int first_sample,
                    int sample,
                    std::vector<double> calibration);

    double get_livetime(std::vector<long long> trigger_times);

};
#endif

