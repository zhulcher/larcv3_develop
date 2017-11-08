#ifndef __THRESHOLDVOXEL3D_CXX__
#define __THRESHOLDVOXEL3D_CXX__

#include "ThresholdVoxel3D.h"
#include "larcv/core/DataFormat/EventVoxel3D.h"
namespace larcv {

  static ThresholdVoxel3DProcessFactory __global_ThresholdVoxel3DProcessFactory__;

  ThresholdVoxel3D::ThresholdVoxel3D(const std::string name)
    : ProcessBase(name)
  {}
    
  void ThresholdVoxel3D::configure(const PSet& cfg)
  {
    _target_producer = cfg.get<std::string>("TargetProducer");
    _output_producer = cfg.get<std::string>("OutputProducer");

    _voxel_value_min = cfg.get<float>("MinThreshold",0.);
    _voxel_value_max = cfg.get<float>("MaxThreshold");
  }

  void ThresholdVoxel3D::initialize()
  {}

  bool ThresholdVoxel3D::process(IOManager& mgr)
  {
    auto const& ev_tensor3d = mgr.get_data<larcv::EventSparseTensor3D>(_target_producer);
    auto ev_output = (larcv::EventSparseTensor3D*)(mgr.get_data("sparse3d",_output_producer));

    ev_output->meta(ev_tensor3d.meta());

    for(auto const& vox : ev_tensor3d.as_vector()) {
      if(vox.value() < _voxel_value_min) continue;
      ((VoxelSet*)ev_output)->emplace(vox.id(), std::min(vox.value(),_voxel_value_max), false);
    }
    return true;
  }

  void ThresholdVoxel3D::finalize()
  {}

}
#endif