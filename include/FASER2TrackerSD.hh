#ifndef FASER2TRACKERSD_HH
#define FASER2TRACKERSD_HH

#include "FASER2TrackerHit.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"
#include "G4ThreeVector.hh"
#include  "G4LorentzVector.hh"

typedef G4THitsCollection<FASER2TrackerHit> FASER2TrackerHitsCollection;

class FASER2TrackerSD : public G4VSensitiveDetector
{
public:
  FASER2TrackerSD(G4String);
  ~FASER2TrackerSD();

  void Initialize(G4HCofThisEvent *HCE);
  /// Temporary map of hits is stored in hit collection, to be retrieved
  /// for analysis by the event action
  // void EndOfEvent(G4HCofThisEvent *HCE);

  G4bool ProcessHits(G4Step*, G4TouchableHistory*);

private:
  /// Hit collection stored in the event, filled in at the end of event based
  /// on temporary hits
  FASER2TrackerHitsCollection *fHitCollection = nullptr;
  /// ID of hit collection
  G4int fHCID = -1;
  // Container to store the track IDs of the tracks which have hit this SD
  // Structure is map<sensorID, vector<trackID>>
  // We have this so that each particle hits each tracker only once - subsequent hits to the same tracker will not be recorded
  std::map<G4int, std::vector<G4int>> fTrackIDRecord; 
};

#endif