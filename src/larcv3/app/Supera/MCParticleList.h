/**
* \file MCParticleList.h
*
* \ingroup Supera
*
* \brief Class def header for a class MCParticleList
*
* @author kterao
*/

/** \addtogroup Supera

@{*/
#ifndef SUPERA_MCPARTICLELIST_H
#define SUPERA_MCPARTICLELIST_H

#include <iostream>
#include <vector>
#include <algorithm>

//#include "EDepSim/TG4Trajectory.h" Zach

#ifdef __has_include
#if __has_include("larcv3/core/dataformat/Particle.h")
#define larcv larcv3
#endif
#endif

namespace larcv
{
  class Particle;
}

namespace supera {

  /**
  \class MCParticleList
  User defined class MCParticleList ... these comments are used to generate
  doxygen documentation!
  */
  class MCParticleList{

  public:

    /// Default constructor
    MCParticleList()
      : _run(-1), _event(-1)
    {}

    /// Default destructor
    ~MCParticleList(){}

    void Update(const std::vector<larcv::Particle> & particles, int run, int event);

    const std::vector<int>& PdgCode()        const { return _pdgcode_v;        }
    const std::vector<int>& ParentIndex()    const { return _parent_index_v;   }
    const std::vector<int>& ParentTrackId()  const { return _parent_trackid_v; }
    const std::vector<int>& ParentPdgCode()  const { return _parent_pdg_v;     }
    const std::vector<int>& TrackIdToIndex() const { return _trackid2index;    }

  private:
    int _run;   ///< a unique event identifier (together with _event)
    int _event; ///< a unique event identifier (together with _run)
    std::vector<int> _trackid_v;        ///< Track ID, index = std::vector<simb::MCParticle> index
    std::vector<int> _pdgcode_v;        ///< PDG code, index = std::vector<simb::MCParticle> index
    std::vector<int> _parent_index_v;   ///< Parent index in std::vector<simb::MCParticle> index, index = std::vector<simb::MCParticle> index
    std::vector<int> _parent_trackid_v; ///< Parent track ID, index = std::vector<simb::MCParticle> index
    std::vector<int> _parent_pdg_v;     ///< Parent PDG, index = std::vector<simb::MCParticle> index
    std::vector<int> _trackid2index;    ///< TrackID => std::vector<simb::MCParticle> index converter
  };
}

#endif
/** @} */ // end of doxygen group
