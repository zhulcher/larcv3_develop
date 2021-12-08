#ifndef TG4PrimaryVertex_hxx_seen
#define TG4PrimaryVertex_hxx_seen

#include <string>
#include <vector>

#include "SuperaLorentz.h"

#if __has_include("larcv3/core/dataformat/Particle.h")
#define larcv larcv3
#endif

class TG4PrimaryVertex;
class TG4PrimaryParticle;

namespace EDepSim {class PersistencyManager;}
typedef std::vector<TG4PrimaryVertex> TG4PrimaryVertexContainer;

/// A class to save a G4 primary vertex into a rooooooooooot output file without linking
/// to geant.
class TG4PrimaryVertex {
    friend class EDepSim::PersistencyManager;
public:
    typedef std::vector<TG4PrimaryParticle> PrimaryParticles;

    TG4PrimaryVertex(void)
        : Position(0,0,0,0), GeneratorName("none"),
          InteractionNumber(0), CrossSection(0.0), DiffCrossSection(0.0),
          Weight(0.0), Probability(0.0) {}
    virtual ~TG4PrimaryVertex();

    /// The initial position of the particle.
    const larcv::SLorentzvector& GetPosition() const {return Position ;}

    /// The name of the generator that created this vertex.
    const char* GetGeneratorName() const {return GeneratorName.c_str();}

    /// The reaction that created this vertex.
    const char* GetReaction() const {return Reaction.c_str();}

    /// The name of the input file.
    const char* GetFilename() const {return Filename.c_str();}

    /// The index (or identifier) of the interaction in the kinematics file.
    int GetInteractionNumber() const {return InteractionNumber;}

    /// The cross section for the reaction that created this vertex.
    double GetCrossSection() const {return CrossSection;}

    /// The differential cross section for the kinematics of the reaction that
    /// created this vertex.
    double GetDiffCrossSection() const {return DiffCrossSection;}

    /// The weight of the interaction.  This will be set to one if the
    /// interaction is not reweighted.  If the vertex is oversampled, this
    /// will be less than one.
    double GetWeight() const {return Weight;}

    /// The overall probability of the interaction that created this vertex.
    /// This includes the effect of the cross section, path length through the
    /// material, etc.  This should be one if it is not filled.
    double GetProbability() const {return Probability;}

    /// The PrimaryVertex points for this PrimaryVertex.
    PrimaryParticles Particles;

    /// The informational vertices associated with this vertex.
    TG4PrimaryVertexContainer Informational;

// The public fields are deprecated but still supported by default in the
// current version.
#define EDEPSIM_USE_PUBLIC_FIELDS
    
#if defined(EDEPSIM_USE_PUBLIC_FIELDS)&&!defined(EDEPSIM_FORCE_PRIVATE_FIELDS)&&!defined(__CINT__)
public:
#ifdef EDEPSIM_WARN_PUBLIC_FIELDS
#warning Using deprecated public fields.  Please consider using the accessor.  For example, to access PrimaryId, use GetPrimaryId().
#endif
#else
private:
#endif
    
    /// The initial position of the particle.
    larcv::SLorentzvector Position;

    /// The name of the generator that created this vertex.
    std::string GeneratorName;

    /// The reaction that created this vertex.
    std::string Reaction;

    /// The name of the input file.
    std::string Filename;

    /// The index (or identifier) of the interaction in the kinematics file.
    int InteractionNumber;

    /// The cross section for the reaction that created this vertex.
    float_t CrossSection;

    /// The differential cross section for the kinematics of the reaction that
    /// created this vertex.
    float_t DiffCrossSection;

    /// The weight of the interaction.  This will be set to one if the
    /// interaction is not reweighted.  If the vertex is oversampled, this
    /// will be less than one.
    float_t Weight;

    /// The overall probability of the interaction that created this vertex.
    /// This includes the effect of the cross section, path length through the
    /// material, etc.  This should be one if it is not filled.
    float_t Probability;

};

/// A class to save a G4 primary particle into a rooooooooooot output file without
/// linking to geant.
class TG4PrimaryParticle {
    friend class EDepSim::PersistencyManager;
public:
    TG4PrimaryParticle(void)
        : TrackId(-1), PDGCode(0), Momentum(0,0,0,0) {}
    virtual ~TG4PrimaryParticle();

    /// The Track Id of the matching trajectory.  Particles that are not
    /// tracked will have negative track id values.
    int GetTrackId() const {return TrackId;}

    /// The name of the particle.
    const char* GetName() const {return Name.c_str();}
    
    /// The PDG code of the particle.
    int GetPDGCode() const {return PDGCode;}

    /// The initial momentum of the particle
    const larcv::SLorentzvector& GetMomentum() const {return Momentum;}

#if defined(EDEPSIM_USE_PUBLIC_FIELDS)&&!defined(EDEPSIM_FORCE_PRIVATE_FIELDS)&&!defined(__CINT__)
public:
#ifdef EDEPSIM_WARN_PUBLIC_FIELDS
#warning Using deprecated public fields.  Please consider using the accessor.  For example, to access PrimaryId, use GetPrimaryId().
#endif
#else
private:
#endif
    
    /// The Track Id of the matching trajectory.  Particles that are not
    /// tracked will have negative track id values.
    int TrackId;

    /// The name of the particle.
    std::string Name;
    
    /// The PDG code of the particle.
    int PDGCode;

    /// The initial momentum of the particle
    larcv::SLorentzvector Momentum;

};
#endif
