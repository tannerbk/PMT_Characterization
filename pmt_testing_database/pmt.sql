CREATE TABLE pmt_information(
    key                     bigserial PRIMARY KEY,
    timestamp               timestamp with time zone default now(),
    source                  text, -- NOT NULL
    pmt_type                text, -- NOT NULL
    pmt_id                  text, -- NOT NULL
    high_voltage            smallint, -- NOT NULL
    tts_sigma               real,
    late_pulsing_pct        real,
    after_pulsing_pct       real,
    pre_pulsing_pct         real,
    dark_rate               real,
    charge_peak             real,
    charge_width            real,
    charge_peak_to_valley   real,
    high_charge_pct         real
    -- entries              bigint
    -- magnetic compensation real
)
