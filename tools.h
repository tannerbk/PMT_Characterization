#include <cmath>
#include <TMath.h>

#include "H5Cpp.h"

#ifndef TOOLS
#define TOOLS

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

using namespace std;

/* Stuct for HDF5 information */
const int RANK_OUT = 2;

struct DataCluster{
  DataSet *dataset; // Dataset pointer
  DataSpace dataspace; // DataSet's DataSpace
  DataSpace memspace; // MemSpace Object for Data Extraction
  hsize_t offset[RANK_OUT]; // Data Extraction Parameters...
  hsize_t count[RANK_OUT];
  hsize_t offset_out[RANK_OUT];
  hsize_t count_out[RANK_OUT];
  unsigned long trace_length; // Length of a Scope Trace
  unsigned long n_traces; // Number of traces in DataSet
  char * data_out; // Pointer to Data Buffer
};

typedef struct DataCluster DataCluster;

/* DataCluster Methods */
DataCluster * Init_Data(DataSet *dataset);
int Read_Trace(DataCluster *datacluster, unsigned long trace_index);

/* Lognormal fits */
double Lognormal(double t, double tau, double sigma, double mag);
double SingleLognormal(double* x, double* par);
double DoubleLognormal(double* x, double* par);
double TripleLognormal(double* x, double* par);

#endif
