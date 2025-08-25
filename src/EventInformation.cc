#include "EventInformation.hh"
#include "generators/GeneratorVertexMetadata.hh"
#include <vector>

EventInformation::EventInformation()
{}

EventInformation::EventInformation(std::vector<GeneratorVertexMetadata> genMetadata)
{
  fGenMetadata = genMetadata;
}

EventInformation::~EventInformation()
{}

void EventInformation::Print() const 
{
  G4cout << G4endl;
  G4cout << "EventInformation: " << fGenMetadata.size() << " vertex(es)" << G4endl;
  for(int i=0; i<fGenMetadata.size(); i++ )
  {
    G4cout << "Vertex Information: " << i << G4endl
      << "Generator " << fGenMetadata[i].generatorType << G4endl
      << "Process name : " << fGenMetadata[i].processName << G4endl
      << "Initiator PDG : " << fGenMetadata[i].pdg << G4endl
      << "Initiator p4 : (" << fGenMetadata[i].p4.X() << ", " << fGenMetadata[i].p4.Y() << ", " 
                            << fGenMetadata[i].p4.Z() << ", " << fGenMetadata[i].p4.E() << ")" << G4endl
      << "Initiator x4 : (" << fGenMetadata[i].x4.X() << ", " << fGenMetadata[i].x4.Y() << ", " 
                            << fGenMetadata[i].x4.Z() << ", " << fGenMetadata[i].x4.T() << ")" << G4endl;
  }
}

