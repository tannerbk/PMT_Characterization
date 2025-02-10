Start by compiling the code by doing:
```
cd src/
make
```

The PMT characterization code assumes at three channels have been digitized in the data; (1) the measurement PMT (2) the trigger PMT and (3) an empty channel.

We run the PMT characterization code as:
```
./src/run_pmt_characterization data.txt output lappd_0 gr0 ch1 gr0 ch0 gr0 ch2 100 0
``` 
where `data.txt` lists the h5 files you want to analyze, output is the name of the produce `.root` file, `lappd_0` is the name of the digitzer, and the list of groups and channels keeps track of the measurement PMT, the trigger PMT, and the empty channel. The last two argument set (1) the pedestal window and (2) the configuration of the setup [0=Cherenkov source, 1=LED source].

For the standard configuration we can run the code as:
```
python3 analysis.py -d path_to_directory -v high_voltage -i pmt_id -k channel
```
