#include "read_data.hh"
#include <string.h>

void ReadData::Read_Trace(DataCluster *datacluster,
                         unsigned long trace_index){
    /* 
    Updates a DataCluster datacluster so that its
    buffer contains trace number trace_index 
    */
    datacluster->offset[0] = (hsize_t)trace_index;

    datacluster->dataspace.selectHyperslab(H5S_SELECT_SET,
                                           datacluster->count,
                                           datacluster->offset);

    datacluster->memspace.selectHyperslab(H5S_SELECT_SET,
                                          datacluster->count_out,
                                          datacluster->offset_out);

    datacluster->dataset->read(datacluster->data_out,
                               PredType::NATIVE_UINT32,
                               datacluster->memspace,
                               datacluster->dataspace);
}

DataCluster * ReadData::Init_Data(DataSet *dataset){
    /*
    Initialize the DataCluster
    */
    DataCluster * datacluster = new DataCluster[1];

    // Get dataspace of the dataset
    datacluster->dataset = dataset;
    datacluster->dataspace = datacluster->dataset->getSpace();

    hsize_t dims_out[2];
    datacluster->dataspace.getSimpleExtentDims(dims_out, NULL);
    datacluster->trace_length = (unsigned long)(dims_out[1]);
    datacluster->n_traces = (unsigned long)(dims_out[0]);

    datacluster->data_out = new unsigned int[datacluster->trace_length];
    for (unsigned long i = 0; i < datacluster->trace_length; i++){
        datacluster->data_out[i]= 0;
    }

    // Define hyperslab in the dataset.
    datacluster->offset[0] = 0;
    datacluster->offset[1] = 0;
    datacluster->count[0] = 1;
    datacluster->count[1] = datacluster->trace_length;
    datacluster->dataspace.selectHyperslab(H5S_SELECT_SET, datacluster->count, datacluster->offset);

    // Define the memory dataspace.
    hsize_t dimsm[2];
    dimsm[0] = dims_out[0];
    dimsm[1] = dims_out[1];
    datacluster->memspace = DataSpace(RANK_OUT, dimsm);

    // Define memory hyperslab.
    datacluster->offset_out[0] = 0;
    datacluster->offset_out[1] = 0;
    datacluster->count_out[0] = 1;
    datacluster->count_out[1] = datacluster->trace_length;
    datacluster->memspace.selectHyperslab(H5S_SELECT_SET, datacluster->count_out, datacluster->offset_out);

    return datacluster;
}

void ReadData::Read_Trace_1D(DataSet dataset,
                             const int ntraces,
                             std::vector<int> &index){
    /*
    Read in the 1D array
    */
    H5::DataSpace dataspace = dataset.getSpace();

    const int rank = dataspace.getSimpleExtentNdims();
    hsize_t offset[rank];
    hsize_t count[rank];
    offset[0] = 0;
    count[0]  = ntraces;

    dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);
    hsize_t dimsm[rank];
    dimsm[0] = ntraces;
    DataSpace memspace(rank, dimsm);

    hsize_t offset_out[rank];
    hsize_t count_out[rank];
    offset_out[0] = 0;
    count_out[0]  = ntraces;
    memspace.selectHyperslab(H5S_SELECT_SET, count_out, offset_out);

    int data[ntraces];
    dataset.read(data, PredType::NATIVE_INT, memspace, dataspace);
    for(int i = 0; i < ntraces; i++){
        index[i] = data[i];
    }
}

void ReadData::Read_Trace_1D_Long(DataSet dataset,
                                  const int ntraces,
                                  std::vector<long long> &index){
    /*
    Read in the 1D array
    */
    H5::DataSpace dataspace = dataset.getSpace();

    const int rank = dataspace.getSimpleExtentNdims();
    hsize_t offset[rank];
    hsize_t count[rank];
    offset[0] = 0;
    count[0]  = ntraces;

    dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);
    hsize_t dimsm[rank];
    dimsm[0] = ntraces;
    DataSpace memspace(rank, dimsm);

    hsize_t offset_out[rank];
    hsize_t count_out[rank];
    offset_out[0] = 0;
    count_out[0]  = ntraces;
    memspace.selectHyperslab(H5S_SELECT_SET, count_out, offset_out);

    long long data[ntraces];
    dataset.read(data, PredType::NATIVE_LLONG, memspace, dataspace);
    for(int i = 0; i < ntraces; i++){
        index[i] = data[i];
    }
}

void ReadData::Read_Trace_2D_Float(DataSet dataset,
                                   std::vector<std::vector<float> > &index){

    H5::DataSpace dataspace = dataset.getSpace();

    hsize_t dims_out[2];
    dataspace.getSimpleExtentDims(dims_out, NULL);

    // Define hyperslab
    hsize_t offset[2];
    hsize_t count[2];
    offset[0] = 0;
    offset[1] = 0;
    count[0] = dims_out[0];
    count[1] = dims_out[1];
    dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);

    // Define the memory dataspace
    hsize_t dimsm[2];
    dimsm[0] = dims_out[0];
    dimsm[1] = dims_out[1];
    DataSpace memspace(2, dimsm);

    // Read the data
    hsize_t offset_out[2];
    hsize_t count_out[2];
    offset_out[0] = 0;
    offset_out[1] = 0;
    count_out[0] = dims_out[0];
    count_out[1] = dims_out[1];
    memspace.selectHyperslab(H5S_SELECT_SET, count_out, offset_out);

    float data[dims_out[0]][dims_out[1]];
    memset(data, 0, dims_out[0]*dims_out[1]*sizeof(float));
    dataset.read(data, PredType::NATIVE_FLOAT, memspace, dataspace);

    for(size_t i = 0; i < dims_out[0]; i++){
        for(size_t j = 0; j < dims_out[1]; j++){
            index[i][j] = data[i][j];
        }
    }
}

