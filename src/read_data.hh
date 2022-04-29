#include <vector>
#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

#ifndef _READ_DATA_H_
#define _READ_DATA_H_

const int RANK_OUT = 2; // Data Rank

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
    unsigned int * data_out; // Pointer to Data Buffer
};

class ReadData {

public:

    void Read_Trace(DataCluster *datacluster, unsigned long trace_index);
    DataCluster * Init_Data(DataSet *dataset);
    void Read_Trace_1D(DataSet dataset, const int ntraces, std::vector<int> &index);
    void Read_Trace_1D_Long(DataSet dataset, const int ntraces, std::vector<long long> &index);
    void Read_Trace_2D_Float(DataSet dataset, std::vector<std::vector<float> > &index);

};
#endif

