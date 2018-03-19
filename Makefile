ReadData_C:
	h5c++ -o pmtanalysis pmt_characterization.cpp tools.cpp -I$(ROOTSYS)/include -L$(ROOTSYS)/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lMinuit -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -lhdf5_cpp -g3 -std=c++0x -fpic -mcmodel=medium -Wall
