#include "PrimaryParticleInformation.hh"

PrimaryParticleInformation::PrimaryParticleInformation(G4int aID, G4int aPDG, G4double aMass, G4double aCharge,
    G4ThreeVector aMomentum, G4ThreeVector aVertex) :
  fPartID(aID), fPDG(aPDG), fMass(aMass), fCharge{aCharge},
  fMomentumMC(aMomentum), fVertexMC(aVertex)
{}

PrimaryParticleInformation::~PrimaryParticleInformation()
{}

void PrimaryParticleInformation::Print() const {
  G4cout<<G4endl;
  G4cout<<"PrimaryParticleInformation: PDG code "<<fPDG<<G4endl
    <<"Particle unique ID : "<<fPartID<<G4endl
    <<"MC momentum : "<<fMomentumMC<<" MeV"<<G4endl
    <<"MC vertex : "<<fVertexMC<<" mm"<<G4endl;
}

