%module dataformat
%{
/* Put header files here or function declarations like below */

// Basic Types:
#include "larcv3/core/dataformat/DataFormatTypes.h"
#include "larcv3/core/dataformat/EventID.h"
#include "larcv3/core/dataformat/Point.h"
#include "larcv3/core/dataformat/BBox.h"
#include "larcv3/core/dataformat/Vertex.h"
#include "larcv3/core/dataformat/Particle.h"
#include "larcv3/core/dataformat/ImageMeta.h"
#include "larcv3/core/dataformat/Image2D.h"
#include "larcv3/core/dataformat/Voxel.h"

// IO Interface classes:
#include "larcv3/core/dataformat/DataProductFactory.h"
#include "larcv3/core/dataformat/EventBase.h"
#include "larcv3/core/dataformat/EventParticle.h"
#include "larcv3/core/dataformat/EventImage2D.h"
#include "larcv3/core/dataformat/IOManager.h"
#include "larcv3/core/dataformat/EventSparseTensor.h"
#include "larcv3/core/dataformat/EventSparseCluster.h"


using namespace larcv3;
%}


%include "std_string.i"
%include "std_vector.i"

// Basic Types:
%include "larcv3/core/dataformat/DataFormatTypes.h"
%include "larcv3/core/dataformat/EventID.h"
%include "larcv3/core/dataformat/Point.h"
%include "larcv3/core/dataformat/BBox.h"
%include "larcv3/core/dataformat/Vertex.h"
%include "larcv3/core/dataformat/Particle.h"
%include "larcv3/core/dataformat/ImageMeta.h"
%include "larcv3/core/dataformat/Image2D.h"
%include "larcv3/core/dataformat/Voxel.h"

// IO Interface classes:
%include "larcv3/core/dataformat/DataProductFactory.h"
%include "larcv3/core/dataformat/EventBase.h"
%include "larcv3/core/dataformat/EventParticle.h"
%include "larcv3/core/dataformat/EventImage2D.h"
%include "larcv3/core/dataformat/IOManager.h"
%include "larcv3/core/dataformat/EventSparseTensor.h"
%include "larcv3/core/dataformat/EventSparseCluster.h"


// Define instantiations of LARCV3 templates:
%template(Point2D) larcv3::Point<2>;
%template(Point3D) larcv3::Point<3>;
%template(BBox2D) larcv3::BBox<2>;
%template(BBox3D) larcv3::BBox<3>;
%template(ImageMeta2D) larcv3::ImageMeta<2>;
%template(ImageMeta3D) larcv3::ImageMeta<3>;
%template(SparseTensor2D) larcv3::SparseTensor<2>;
%template(SparseTensor3D) larcv3::SparseTensor<3>;
%template(EventSparseTensor2D) larcv3::EventSparseTensor<2>;
%template(EventSparseTensor3D) larcv3::EventSparseTensor<3>;
%template(SparseCluster2D) larcv3::SparseCluster<2>;
%template(SparseCluster3D) larcv3::SparseCluster<3>;
%template(EventSparseCluster2D) larcv3::EventSparseCluster<2>;
%template(EventSparseCluster3D) larcv3::EventSparseCluster<3>;


// Add templates for the dataproducts that get wrapped vectors:
%template(VectorOfDouble) std::vector<double>;
%template(VectorOfFloat) std::vector<float>;
%template(VectorOfSizet) std::vector<size_t>;
%template(VectorOfVoxel) std::vector<larcv3::Voxel>;
%template(VectorOfVoxelSet) std::vector<larcv3::VoxelSet>;
%template(VectorOfSparse2D) std::vector<larcv3::SparseTensor2D>;
%template(VectorOfSparse3D) std::vector<larcv3::SparseTensor3D>;
%template(VectorOfCluster2D) std::vector<larcv3::SparseCluster2D>;
%template(VectorOfCluster3D) std::vector<larcv3::SparseCluster3D>;
%template(VectorOfParticle) std::vector<larcv3::Particle>;
%template(VectorOfImage2D) std::vector<larcv3::Image2D>;



