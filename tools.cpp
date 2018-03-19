#include "tools.h"

DataCluster * Init_Data(DataSet *dataset){

  DataCluster * datacluster = new DataCluster[1];

  /*
   * Get dataspace of the dataset.
   */
  datacluster->dataset = dataset;
  datacluster->dataspace = datacluster->dataset->getSpace();

  /*
   * Get the dimension size of each dimension in the dataspace and
   * display them.
   */

  hsize_t dims_out[2];
  datacluster->dataspace.getSimpleExtentDims( dims_out, NULL);
  datacluster->trace_length = (unsigned long)(dims_out[1]);
  datacluster->n_traces = (unsigned long)(dims_out[0]);

  // Data Buffer
  datacluster->data_out = new char[datacluster->trace_length]; // Scope data is size char
  for (unsigned long i = 0; i < datacluster->trace_length; i++) datacluster->data_out[i]= 0;

  /*
   * Define hyperslab in the dataset.
   */

  datacluster->offset[0] = 0;
  datacluster->offset[1] = 0;
  datacluster->count[0] = 1;
  datacluster->count[1] = datacluster->trace_length;
  datacluster->dataspace.selectHyperslab( H5S_SELECT_SET, datacluster->count, datacluster->offset );

  /*
   * Define the memory dataspace.
   */

  hsize_t dimsm[2]; /* memory space dimensions */
  dimsm[0] = dims_out[0];
  dimsm[1] = dims_out[1];
  datacluster->memspace = DataSpace( RANK_OUT, dimsm );

  /*
   * Define memory hyperslab.
   */

  datacluster->offset_out[0] = 0;
  datacluster->offset_out[1] = 0;
  datacluster->count_out[0] = 1;
  datacluster->count_out[1] = datacluster->trace_length;
  datacluster->memspace.selectHyperslab( H5S_SELECT_SET, datacluster->count_out, datacluster->offset_out );

  return datacluster;
}


/* Method: Read_Trace(DataCluster *datacluster, unsigned long trace_index)
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Updates a DataCluster datacluster so that its buffer contains trace number trace_index
 */

int Read_Trace(DataCluster *datacluster, unsigned long trace_index){
  datacluster->offset[0]= (hsize_t)trace_index;
  datacluster->dataspace.selectHyperslab( H5S_SELECT_SET, datacluster->count, datacluster->offset );
  datacluster->memspace.selectHyperslab( H5S_SELECT_SET, datacluster->count_out, datacluster->offset_out );
  datacluster->dataset->read( datacluster->data_out, PredType::NATIVE_CHAR, datacluster->memspace, datacluster->dataspace );
  return 0; // No protection...
}

double Lognormal(double t, double tau, double sigma, double mag){
    double q = exp(-0.5 * pow(log(t / tau) / sigma, 2.0));
    q *= mag / (t * sigma * sqrt(2 * TMath::Pi()));
    return -q;
}

double SingleLognormal(double* x, double* par){
    return Lognormal(x[0], par[0], par[1], par[2]);
}

double DoubleLognormal(double* x, double* par){
    double q = Lognormal(x[0], par[0], par[1], par[2]);
    q += Lognormal(x[0], par[3], par[4], par[5]);
    return q;
}

double TripleLognormal(double* x, double* par){
    double q = Lognormal(x[0], par[0], par[1], par[2]);
    q += Lognormal(x[0], par[3], par[4], par[5]);
    q += Lognormal(x[0], par[6], par[7], par[8]);
    return q;
}
