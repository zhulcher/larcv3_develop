/**
 * \file EventVoxel2D.h
 *
 * \ingroup core_DataFormat
 *
 * \brief Class def header for a class larcv3::EventSparseTensor and larcv3::EventClusterPixel2D
 *
 * @author cadams
 */

/** \addtogroup core_DataFormat

    @{*/
#ifndef __LARCV3DATAFORMAT_VOXEL_SERIALIZATION_HELPER_H
#define __LARCV3DATAFORMAT_VOXEL_SERIALIZATION_HELPER_H

#include <iostream>
#include "larcv3/core/dataformat/DataProductFactory.h"
#include "larcv3/core/dataformat/Voxel.h"
#include "larcv3/core/dataformat/ImageMeta.h"


namespace larcv3 {

  // Create a suite of functions to manage the shared IO of sparse tensors and clusters

  
  // This class keeps a copy of the meta for an object, so it is managed here.
  
  template<size_t dimension>
  class VoxelSerializationHelper {
  
  public:

    VoxelSerializationHelper();

    // Write the projection meta to the group
    void write_projection_meta(H5::Group *group);

    // Read the projection meta from the group
    void read_projection_meta(H5::Group *group);

    // Write a set of voxels, packaged in the way proscribed by meta, to the group:
    void write_voxels(H5::Group * group, const std::vector<std::vector<larcv3::VoxelSet> > & voxels );

    // Read an entry of voxels from a group
    void read_voxels(H5::Group * group, size_t entry, std::vector<std::vector<larcv3::VoxelSet> > & voxels);

    // Initialization is the same for sparse image and cluster set
    void initialize_voxel_group(H5::Group * group);

    // Function to write a image meta as an attribute
    void write_image_meta(H5::H5Object * obj);

    // Function to read all meta from an object
    void read_image_meta(H5::H5Object *obj);

    bool initialized(){return _initialized;}

    void initialize_for_read(H5::Group * obj);
    void initialize_for_write(H5::Group * obj);

    std::vector<larcv3::ImageMeta<dimension> > image_meta;
    ImageMeta2D projection_meta;

  private:
    bool _initialized;
  };




}

#endif
/** @} */ // end of doxygen group
