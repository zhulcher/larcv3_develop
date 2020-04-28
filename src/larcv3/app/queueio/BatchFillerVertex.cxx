#ifndef __LARCV3THREADIO_BATCHFILLERVERTEX_CXX__
#define __LARCV3THREADIO_BATCHFILLERVERTEX_CXX__

#include "BatchFillerVertex.h"
#include "larcv3/core/dataformat/EventTensor.h"
#include "larcv3/core/dataformat/EventParticle.h"
#include <random>

namespace larcv3 {

  static BatchFillerVertexProcessFactory __global_BatchFillerVertexProcessFactory__;

  BatchFillerVertex::BatchFillerVertex(const std::string name)
    : BatchFillerTemplate<float>(name)
  {}

  void BatchFillerVertex::configure(const PSet& cfg)
  {
    _part_producer = cfg.get<std::string>("ParticleProducer");
  }

  void BatchFillerVertex::initialize()
  {}

  void BatchFillerVertex::_batch_begin_()
  {
    if (!batch_data().dim().empty() && (int)(batch_size()) != batch_data().dim().front()) {
      auto dim = batch_data().dim();
      dim[0] = batch_size();
      this->set_dim(dim);
    }
  }

  void BatchFillerVertex::_batch_end_()
  {
    if (logger().level() <= msg::kINFO) {
      LARCV_INFO() << "Total data size: " << batch_data().data_size() << std::endl;

      std::vector<size_t> ctr_v(_vtx_dim, 0);
      auto const& data = batch_data().data();
      for (size_t i = 0; i < data.size(); ++i) {
        if (data[i] < 1.) continue;
        ctr_v[i % _vtx_dim] += 1;
      }
      std::stringstream ss;
      ss << "Class fractions (0";
      for (size_t i = 1; i < _vtx_dim; ++i) ss << "," << i;
      ss << ") ... (" << ctr_v[0];
      for (size_t i = 1; i < _vtx_dim; ++i) ss << "," << ctr_v[i];
      ss << ")";
      LARCV_INFO() << ss.str() << std::endl;
    }
  }

  void BatchFillerVertex::finalize()
  {}

  bool BatchFillerVertex::process(IOManager & mgr)
  {

    auto const& event_part = mgr.get_data<larcv3::EventParticle>(_part_producer);

    // Refresh the dimension:
    std::vector<int> dim(2);
    dim[0] = batch_size();
    dim[1] = _vtx_dim;
    set_dim(dim);
    // In this case, set the dense dim as the same:
    set_dense_dim(dim);

    // labels
    auto const& part_v = event_part.as_vector();
    
    // vertexes
    _entry_data.resize(_vtx_dim, 0);
    for (auto const& part : part_v) {
      // _entry_data = {part.x(), part.y(), part.z()};
      _entry_data[0] = part.x();
      _entry_data[1] = part.y();
      _entry_data[2] = part.z();
      LARCV_DEBUG() << "Vertex for particle with PDG" << part.pdg_code() << " is " 
                    << "[" << part.x() << "," << part.y() << "," << part.z() << "]" << std::endl;
    }

    set_entry_data(_entry_data);

    return true;
  }

}
#endif