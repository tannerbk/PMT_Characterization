CREATE TABLE pmt_information(
    key                      bigserial PRIMARY KEY,
    timestamp                timestamp with time zone default now(),
    entries                  bigint NOT NULL,
    source                   text NOT NULL,
    pmt_type                 text NOT NULL,
    pmt_id                   text NOT NULL,
    high_voltage             smallint NOT NULL, -- V
    tts_sigma                real, -- ns
    tts_sigma_err            real, -- ns
    late_pulsing_pct         real,
    after_pulsing_pct        real,
    pre_pulsing_pct          real,
    dark_rate                real, -- Hz
    charge_peak              real, -- pC
    charge_width             real, -- pC
    charge_peak_to_valley    real,
    high_charge_pct          real,
    threshold                real,
    coincidence_rate         real,     
    settling_time            real, -- hr
    trigger_threshold        real, -- mV
    trigger_q_cut            real, -- pC
    magnetic_compensation    text NOT NULL,
    comment                  text
)

