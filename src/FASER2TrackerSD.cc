#include "FASER2TrackerSD.hh"
#include "G4SystemOfUnits.hh"

FASER2TrackerSD::FASER2TrackerSD(G4String name) :
  G4VSensitiveDetector(name) {
  G4cout << "creating a sensitive detector with name: " << name << G4endl;
  collectionName.insert(name);
}

FASER2TrackerSD::~FASER2TrackerSD(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void FASER2TrackerSD::Initialize(G4HCofThisEvent *HCE) {
  fHitCollection = new FASER2TrackerHitsCollection(GetName(), collectionName[0]);

  if (fHCID < 0)
    fHCID = GetCollectionID(0);
  HCE->AddHitsCollection(fHCID, fHitCollection);

  fTmpHits.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FASER2TrackerSD::EndOfEvent(G4HCofThisEvent *) {
  for (auto it = fTmpHits.begin(); it != fTmpHits.end(); ++it)
    fHitCollection->insert(*it);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool FASER2TrackerSD::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist){

  G4Track* track = aStep->GetTrack();
  //track->SetTrackStatus(fStopAndKill);
  G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
  G4StepPoint *postStepPoint = aStep->GetPostStepPoint();

  G4ThreeVector posHit = preStepPoint->GetPosition();
  G4int pdgid = track->GetParticleDefinition()->GetPDGEncoding();
  G4double energy = track->GetDynamicParticle()->Get4Momentum().e();

  if (energy/MeV < 2)
  {
    G4cout << "Hit energy too low: " << energy/MeV << " MeV, skipping hit." << G4endl;
    return 0;
  }

  G4double px = track->GetDynamicParticle()->Get4Momentum().px();
  G4double py = track->GetDynamicParticle()->Get4Momentum().py();
  G4double pz = track->GetDynamicParticle()->Get4Momentum().pz();
  G4double m = track->GetDynamicParticle()->Get4Momentum().m();
  G4double charge = track->GetDynamicParticle()->GetCharge();
  G4double time = track->GetDynamicParticle()->Get4Momentum().t();
  G4ThreeVector delta_momentum = aStep->GetDeltaMomentum();
  G4double delta_energy = aStep->GetDeltaEnergy();
  G4int sensor_id = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();

  FASER2TrackerHit* tmpHit = new FASER2TrackerHit();

  G4TouchableHandle touchable = preStepPoint->GetTouchableHandle();
  G4ThreeVector sensorCenterGlobal = touchable->GetTranslation();
  G4double sensorCentreZ = sensorCenterGlobal.z();

  // tmpHit->SetPosition(posHit[0]/mm, posHit[1]/mm, posHit[2]/mm); // in mm
  // fix the hit z-position to be the centre of the sensor - this way every hit on the same sensor has the same z-pos

  tmpHit->SetPosition(posHit[0]/mm, posHit[1]/mm, sensorCentreZ/mm); // in mm
  tmpHit->SetPDGID(pdgid);
  tmpHit->SetEnergy(energy/GeV);
  tmpHit->SetCharge(charge);
  tmpHit->SetPx(px/GeV);
  tmpHit->SetPy(py/GeV);
  tmpHit->SetPz(pz/GeV);
  tmpHit->SetDeltaPx(delta_momentum.x()/GeV);
  tmpHit->SetDeltaPy(delta_momentum.y()/GeV);
  tmpHit->SetDeltaPz(delta_momentum.z()/GeV);
  tmpHit->SetDeltaE(delta_energy/GeV);
  tmpHit->SetMass(m/GeV);
  tmpHit->SetTrackID(track->GetTrackID());
  tmpHit->SetParentID(track->GetParentID());
  tmpHit->SetCopyNumSensor(sensor_id);
  tmpHit->SetT(time/ns);
  tmpHit->SetTrackVertex(track->GetVertexPosition()/mm);
  tmpHit->SetTrackP4(track->GetDynamicParticle()->Get4Momentum()/GeV);
  if (track->GetParentID() == 0) {
    tmpHit->SetIsPrimaryTrack(1);
    tmpHit->SetIsSecondaryTrack(0);
  }
  else {
    tmpHit->SetIsPrimaryTrack(0);
    tmpHit->SetIsSecondaryTrack(1);
  }
  fTmpHits.push_back(tmpHit); 
  
  return 0;
}