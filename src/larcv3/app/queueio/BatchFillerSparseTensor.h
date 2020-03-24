/**
 * \file BatchFillerSparseTensor.h
 *
 * \ingroup ThreadIO
 *
 * \brief Class def header for a class BatchFillerSparseTensor
 *
 * @author cadams
 */

/** \addtogroup ThreadIO

    @{*/
#ifndef __LARCV3THREADIO_BATCHFILLERSPARSETENSOR_H__
#define __LARCV3THREADIO_BATCHFILLERSPARSETENSOR_H__

#include "larcv3/core/processor/ProcessFactory.h"
#include "BatchFillerTemplate.h"

#include "larcv3/core/dataformat/EventSparseTensor.h"

namespace larcv3 {

  /**     
     This class reads sparse tensors from file (or after conversions), and stores them
     in dense format.
  */
  template <size_t dimension>
  class BatchFillerSparseTensor : public BatchFillerTemplate<float> {

  public:

    /// Default constructor
    BatchFillerSparseTensor(const std::string name="BatchFillerSparseTensor");

    /// Default destructor
    ~BatchFillerSparseTensor(){}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();

  protected:

    void _batch_begin_();
    void _batch_end_();

  private:


    std::string _sparse_tensor_producer;
    float _unfilled_voxel_value;
    std::vector<size_t> _slice_v;
    std::vector<size_t> _compression_v;


    std::vector<float>  _entry_data;
    size_t _num_channels;
    bool _allow_empty;
    bool _augment;
  };

  /**
     \class larcv3::BatchFillerSparseTensorFactory
     \brief A concrete factory class for larcv3::BatchFillerSparseTensor
  */
  class BatchFillerSparseTensorProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    BatchFillerSparseTensorProcessFactory() { ProcessFactory::get().add_factory("BatchFillerSparseTensor",this); }
    /// dtor
    ~BatchFillerSparseTensorProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) {
      return new BatchFillerSparseTensor(instance_name);
    }
  };

}

#endif
/** @} */ // end of doxygen group

