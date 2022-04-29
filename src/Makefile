LDFLAGS=-g -L$(ROOTSYS)/lib
LDLIBS=`root-config --libs --cflags` -Wall -lcurl -pthread

.PHONY: all run_pmt_characterization

PMT_OBJS=pmt_characterization.cpp calibration.cpp read_data.cpp tools.cpp

all: run_pmt_characterization

run_pmt_characterization: $(OBJS)
	h5c++ $(LDFLAGS) -o run_pmt_characterization $(PMT_OBJS) $(LDLIBS)

clean:
	rm -f run_pmt_characterization *.o

