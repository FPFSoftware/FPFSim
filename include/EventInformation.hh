#ifndef EventInformation_HH
#define EventInformation_HH

#include "G4VUserEventInformation.hh"
#include "generators/GeneratorVertexMetadata.hh"
#include "globals.hh"
#include <vector>

// Event information
class EventInformation : public G4VUserEventInformation
{
  public:
    
    EventInformation();
    EventInformation(std::vector<GeneratorVertexMetadata> genMetadata);
    virtual ~EventInformation();

    /// Gets vertex metadata full vector
    inline std::vector<GeneratorVertexMetadata> GetEventMetadata() const { return fGenMetadata; }

    /// Gets metadata per vertex
    inline GeneratorVertexMetadata GetMetadataPerVertex(int i) const { return fGenMetadata.at(i); }

    /// Prints the information about the event.
    virtual void Print() const;
    virtual void Print(int i) const;

  private:
    /// Set of vertex metadata
    std::vector<GeneratorVertexMetadata> fGenMetadata;
};

#endif
