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
    G4cout << G4endl;
    G4cout << "Vertex : " << i << G4endl
      << "Generator : " << fGenMetadata[i].generatorType << G4endl
      << "Process name : " << fGenMetadata[i].processName << G4endl
      << "Initiator PDG : " << fGenMetadata[i].pdg << G4endl
      << "Initiator p4 : (" << fGenMetadata[i].p4.x() << ", " << fGenMetadata[i].p4.y() << ", " 
                            << fGenMetadata[i].p4.z() << ", " << fGenMetadata[i].p4.e() << ")" << G4endl
      << "Initiator x4 : (" << fGenMetadata[i].x4.x() << ", " << fGenMetadata[i].x4.y() << ", " 
                            << fGenMetadata[i].x4.z() << ", " << fGenMetadata[i].x4.t() << ")" << G4endl;
  }
}

