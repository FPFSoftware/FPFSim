#ifndef GENERATOR_BASE_HH
#define GENERATOR_BASE_HH

#include <vector>
#include "G4Event.hh"
#include "G4UImessenger.hh"
#include "generators/GeneratorVertexMetadata.hh"

class GeneratorBase
{
  public:

    GeneratorBase() : fGeneratorName("unknown"), fMessenger(nullptr) {}
    virtual ~GeneratorBase() {}

    // Called once (by PrimaryGeneratorAction) to load the external file,
    // open files, set branch addresses, or load histograms, etc..
    virtual void LoadData() = 0;

    // Called for each event to generate primaries
    virtual void GeneratePrimaries(G4Event *event) = 0;

    // return name of current generator
    const G4String GetGeneratorName() { return fGeneratorName; }

    // reset the list of metadata
    void ResetEventMetadata() { fVertexMetadata.clear(); }

    // return full event vertex metadata
    const std::vector<GeneratorVertexMetadata> GetEventMetadata() { return fVertexMetadata; }
    // return single vertex metadata
    const GeneratorVertexMetadata GetEventMetadataPerVertex(G4int i) { return fVertexMetadata.at(i); }

  protected : 

    G4String fGeneratorName; 
    G4UImessenger* fMessenger;
    std::vector<GeneratorVertexMetadata> fVertexMetadata;
};

#endif
