// #ifndef __LARCV3THREADIO_BATCHFILLERSPARSETENSOR_CXX__
// #define __LARCV3THREADIO_BATCHFILLERSPARSETENSOR_CXX__

// #include "BatchFillerSparseTensor.h"

// #include <random>

// namespace larcv3 {

// template <size_t dimension>
// BatchFillerSparseTensor<dimension>::BatchFillerSparseTensor(const std::string name)
//     : BatchFillerTemplate<float>(name) {}


// template <size_t dimension>
// void BatchFillerSparseTensor<dimension>::configure(const PSet& cfg) {
//   LARCV_DEBUG() << "start" << std::endl;
//   _tensor2d_producer = cfg.get<std::string>("TensorProducer");

//   _unfilled_voxel_value = cfg.get<float>("UnfilledVoxelValue", 0.0);
//   _slice_v = cfg.get<std::vector<size_t> >("Channels", _slice_v);

//   LARCV_DEBUG() << "done" << std::endl;
// }

// template <size_t dimension>
// void BatchFillerSparseTensor<dimension>::initialize() {}

// template <size_t dimension>
// void BatchFillerSparseTensor<dimension>::_batch_begin_() {
//   if(!batch_data().dim().empty() && (int)(batch_size()) != batch_data().dim().front()) {
//     LARCV_INFO() << "Batch size changed " << batch_data().dim().front() << "=>" << batch_size() << std::endl;
//     auto dim = batch_data().dim();
//     dim[0] = batch_size();
//     this->set_dim(dim);
//   }
// }

// template <size_t dimension>
// void BatchFillerSparseTensor<dimension>::_batch_end_() {
//   if (logger().level() <= msg::kINFO)
//     LARCV_INFO() << "Total data size: " << batch_data().data_size()
//                  << std::endl;
// }

// template <size_t dimension>
// void BatchFillerSparseTensor<dimension>::finalize() { _entry_data.clear(); }

// template <size_t dimension>
// bool BatchFillerSparseTensor<dimension>::process(IOManager& mgr) {

//   /*

//   We read in sparse tensor (2D, 3D, doesn't matter).  We determine the image size from the meta, 
//   and we set the output datasize to match.

//   */


//   LARCV_DEBUG() << "start" << std::endl;

//   auto const& voxel_data =
//       mgr.get_data<larcv3::EventSparseTensor<dimension> >(_tensor_producer);
//   if (voxel_data.as_vector().empty()) {
//     LARCV_CRITICAL()
//         << "Could not locate non-empty sparse tensor data w/ producer name "
//         << _tensor_producer << std::endl;
//     throw larbys();
//   }

//   _num_channels = _slice_v.size();


//   if (_num_channels != voxel_data.as_vector().size()){
//     LARCV_CRITICAL() << "Number of requested channels does not match number of available channels." << std::endl;
//     throw larbys();
//   }


//   // Next, determine the dimensions.
//   // There is a requirement that all projection_id have the same meta shape
//   auto const& voxel_meta = voxel_data.as_vector().front().meta();


//   std::vector<int> dimensions;

//   for


//   std::vector<int> dense_dim;
//   dense_dim.resize(4);
//   dense_dim[0] = batch_size();
//   dense_dim[1] = voxel_meta.number_of_voxels(0);
//   dense_dim[2] = voxel_meta.number_of_voxels(1);
//   dense_dim[3] = _num_channels;
//   this->set_dense_dim(dense_dim);


//   if (_entry_data.size() != batch_data().entry_data_size())
//     _entry_data.resize(batch_data().entry_data_size(), 0.);


//   // Reset all values to 0.0 (or whatever is specified)
//   for (auto& v : _entry_data) v = _unfilled_voxel_value;

//   // Get the random x/y/z flipping
//   bool flip_cols = false;
//   bool flip_rows = false;
//   if (_augment){
//     flip_cols = bool(rand() % 2);
//     flip_rows = bool(rand() % 2);
//   }



//   for ( auto const& voxel_set : voxel_data.as_vector()){
//     auto & meta = voxel_set.meta();
//     auto projection_id = meta.id();

//     // Check that this projection ID is in the lists of channels:
//     bool found = false;
//     int count = 0;
//     for (auto & channel : _slice_v){
//       if (channel == projection_id){
//         found = true;
//         break;
//       }
//       count ++;
//     }
//     if (!found) continue;
//     size_t max_voxel(voxel_set.size());
//     if (max_voxel > _max_voxels) {
//       max_voxel = _max_voxels;
//       LARCV_INFO() << "Truncating the number of voxels!" << std::endl;
//     }

//     size_t index;
//     int row(0), col(0);
//     int row_mult(1), row_add(0);
//     int col_mult(1), col_add(0);

//     if (flip_rows){
//       row_mult = -1;
//       row_add = meta.rows() - 1;
//     }
//     if (flip_cols){
//       col_mult = -1;
//       col_add = meta.cols() - 1;
//     }

//     // Get all of the indexes:
//     std::vector<size_t> indexes = voxel_set.indexes_vec();

//     // Convert them all to coordinates:
//     std::vector<size_t> coordinates;
//     meta.coordinates(indexes, coordinates);


//     if (_include_values){
//       std::vector<float>  values  = voxel_set.values_vec();
// #ifdef LARCV_OPENMP
//       #pragma omp parallel private(row, col, index) shared(_entry_data, row_mult, row_add, col_mult, col_add)
// #endif
//       for (size_t i_voxel = 0; i_voxel < max_voxel; i_voxel ++) {
//       // for (auto const& voxel : voxel_set.as_vector()) {
//         row = coordinates[i_voxel*2];
//         col = coordinates[i_voxel*2 + 1];

//         col = col_mult*col + col_add;
//         row = row_mult*row + row_add;


//         index = count*(_max_voxels*point_dim) + i_voxel*point_dim;
//         _entry_data.at(index) = row;
//         _entry_data.at(index + 1) = col;
//         _entry_data.at(index + 2) = values[i_voxel];
//       }
//     }
//     else{
// #ifdef LARCV_OPENMP
//       #pragma omp parallel private(row, col, index) shared(_entry_data, row_mult, row_add, col_mult, col_add)
// #endif
//       for (size_t i_voxel = 0; i_voxel < max_voxel; i_voxel ++) {
//       // for (auto const& voxel : voxel_set.as_vector()) {
//         row = coordinates[i_voxel*2];
//         col = coordinates[i_voxel*2 + 1];

//         col = col_mult*col + col_add;
//         row = row_mult*row + row_add;

//         index = count*(_max_voxels*point_dim) + i_voxel*point_dim;
//         _entry_data.at(index) = row;
//         _entry_data.at(index + 1) = col;
//       }
//     }




//   }

//   // record the entry data
//   LARCV_INFO() << "Inserting entry data of size " << _entry_data.size()
//                << std::endl;

//   set_entry_data(_entry_data);

//   return true;
// }

// static BatchFillerSparseTensor2DProcessFactory
//     __global_BatchFillerSparseTensor2DProcessFactory__;


// }
// #endif
