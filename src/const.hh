// Digitizer dynamic range in volts
const double DYNAMIC_RANGE = 1.0;

// Termination resistance
const double termination_ohms = 50.0;

// Threshold in mV for an SPE pulse
const double voltage_threshold = -5.0;

// Threshold in mV for the trigger pulse
const double trigger_voltage_threshold = -10.0;

// Name of samples that pulse stays above threshold
const int nsamples_above_threshold_threshold = 5;

// Number of samples at the beginning of the waveform to skip
const int sample_offset = 25;

// Name of samples a the end of the waveform to skip
const int sample_offset_last = 124;

// Num. of samples to look back at CHESS PMT
const int lookback = 50;

// Num. of samples to look back after trigger PMT pulse
const int lookback_trigger = 50;

// Number of samples to integrate for charge
const int integration_samples_back = 30;

// Number of samples to integrate for charge
const int integration_samples_forward = 60;

// Const. fraction discriminator
const double const_frac_thresh = 0.6;

// Const. fraction discriminator
const double const_frac_thresh_trigger = 0.2;

// Const. thresh on trigger, in mV
const double const_thresh_trigger = -2.0;

// Write the fit waveforms to file
const bool simple_write_waveforms = false;

const bool CONST_FRAC_TRIGGER = false;

// Invalid time identified
const int INVALID = -9999;
