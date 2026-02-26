#include "generators/GeneratorBase.hh"
#include "generators/GPSGenerator.hh"
#include "generators/GeneratorVertexMetadata.hh"

#include "geometry/GeometricalParameters.hh"

#include "G4GeneralParticleSource.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4LorentzVector.hh"
#include "Randomize.hh"

GPSGenerator::GPSGenerator()
{
  fGeneratorName = "gun";
  fMessenger = nullptr; // rely on standard /gps/ directory

  fGPS = new G4GeneralParticleSource();

  // set GPS default values, 
  // this are usually overriden in macro
  G4ParticleDefinition* myParticle;
  myParticle = G4ParticleTable::GetParticleTable()->FindParticle("e-");
  fGPS->SetParticleDefinition(myParticle);
  fGPS->GetCurrentSource()->GetEneDist()->SetMonoEnergy(5*GeV);  // kinetic energy
  fGPS->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
  G4long seeds[2];
  time_t systime = time(NULL);
  seeds[0] = (G4long) systime;
  seeds[1] = (G4long) (systime*G4UniformRand());
  G4Random::setTheSeeds(seeds);
  G4Random::showEngineStatus();
  G4double x0 = G4UniformRand();
  G4double y0 = G4UniformRand();
  G4double z0 = G4UniformRand();
  x0 = GeometricalParameters::Get()->GetFLArEPosition().x() +
       (x0-0.5)*GeometricalParameters::Get()->GetFLArEFidVolSize().x();
  y0 = GeometricalParameters::Get()->GetFLArEPosition().y() +
       (y0-0.5)*GeometricalParameters::Get()->GetFLArEFidVolSize().y();
  z0 = GeometricalParameters::Get()->GetFLArEPosition().z() +
       (z0-0.5)*GeometricalParameters::Get()->GetFLArEFidVolSize().z();;
  fGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Point");
  fGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(x0, y0, z0));

}

GPSGenerator::~GPSGenerator()
{
  delete fGPS;
  delete fMessenger;
}

void GPSGenerator::LoadData()
{ // nothing to load!
} 

void GPSGenerator::GeneratePrimaries(G4Event* anEvent)
{
  // complete line from PrimaryGeneratorAction..
  G4cout << "): General Particle Source ===oooOOOooo===" << G4endl;

  // Generate primary vertex first
  fGPS->GeneratePrimaryVertex(anEvent);

  // preparing to ship metadata
  GeneratorVertexMetadata metadata;
  metadata.generatorType = fGeneratorName;
  metadata.processName = "Gun";
  metadata.weight = 1.0;
  metadata.pdg = fGPS->GetParticleDefinition()->GetPDGEncoding();
  metadata.mass = fGPS->GetParticleDefinition()->GetPDGMass();
  metadata.charge = fGPS->GetParticleDefinition()->GetPDGCharge();
  G4PrimaryVertex* vtx = anEvent->GetPrimaryVertex();
  G4PrimaryParticle* pp = vtx->GetPrimary(0);
  metadata.x4 = G4LorentzVector(vtx->GetX0(),vtx->GetY0(),vtx->GetZ0(),vtx->GetT0());
  metadata.p4 = G4LorentzVector(pp->GetPx(),pp->GetPy(),pp->GetPz(),pp->GetTotalEnergy());
  fVertexMetadata.push_back(metadata);
}



