#ifndef __LARCV3DATAFORMAT_MPIIOMANAGER_CXX__
#define __LARCV3DATAFORMAT_MPIIOMANAGER_CXX__

#include "larcv3/app/mpi_io/MPIIOManager.h"


#ifndef LARCV_MPI
fgdfgsf
#endif

namespace larcv3{

    MPIIOManager::MPIIOManager(IOMode_t mode, std::string name) :
        IOManager(mode, name)
    {
    }

    /// Configuration PSet construction so you don't have to call setter functions
    MPIIOManager::MPIIOManager(const PSet& cfg) : 
        IOManager(cfg)
    {
    }

    /// Configuration PSet file construction so you don't have to call setter functions
    MPIIOManager::MPIIOManager(std::string config_file, std::string name) : 
        IOManager(config_file, name)
    {

    }

    /// Destructor
    MPIIOManager::~MPIIOManager(){
    }

bool MPIIOManager::initialize(){
    
    std::cout << "Entered initialize " << std::endl;
    int ierr = MPI_Init(NULL, NULL) ;

    // int * initialized(nullptr);
    // int error_val =  MPI_Initialized(initialized);
    // std::cout << " * Initialized is  " << * initialized << std::endl;
    // std::cout << " Initialized is  " << initialized << std::endl;
    // std::cout << " error_val is  " << error_val << std::endl;

    // if (! initialized ){
    //     // Initialize the MPI environment
    // }


    // // Get the number of processes
    // int world_size;
    // MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // // Get the rank of the process
    // int world_rank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // std::cout << "world size: " << world_size << std::endl;
    // std::cout << "world rank: " << world_rank << std::endl;
    // return true;
    // MPI_Finalize();
    return true;
}

void MPIIOManager::open_new_input_file(std::string filename){
  


  // try{
  //   _in_open_file = H5::H5File(filename.c_str(), H5F_ACC_RDONLY, H5::FileCreatPropList::DEFAULT, _fapl);
  // }
  // catch ( ... ) {
  //   LARCV_CRITICAL() << "Open attempt failed for a file: " << filename
  //                    << std::endl;
  //   throw larbys();
  // }

  // _active_in_event_id_dataset   = _in_open_file.openGroup("Events").openDataSet("event_id");
  // _active_in_event_id_dataspace = _active_in_event_id_dataset.getSpace();

}


}  // namespace larcv3
#endif
