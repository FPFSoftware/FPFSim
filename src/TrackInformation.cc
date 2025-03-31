#include "TrackInformation.hh"

G4ThreadLocal G4Allocator<TrackInformation> *
                                   aTrackInformationAllocator = 0;

TrackInformation::TrackInformation()
  : G4VUserTrackInformation()
{
  fFromPrimaryPizero = 0;
  fFromFSLPizero = 0;
  fFromPrimaryLepton = 0;
}

TrackInformation::TrackInformation(const G4Track* aTrack) 
  : G4VUserTrackInformation()
{
  fFromPrimaryPizero = 0;
  fFromFSLPizero = 0;
  fFromPrimaryLepton = 0;
}

TrackInformation::~TrackInformation()
{;}

TrackInformation& TrackInformation
::operator =(const TrackInformation& aTrackInfo)
{
  fFromPrimaryPizero = aTrackInfo.fFromPrimaryPizero;
  fFromFSLPizero = aTrackInfo.fFromFSLPizero;
  fFromPrimaryLepton = aTrackInfo.fFromPrimaryLepton;

  return *this;
}

void TrackInformation::Print() const
{
    G4cout << "Is from primary pizero " << fFromPrimaryPizero << G4endl;
    G4cout << "Is from final state lepton decay pizero " << fFromFSLPizero << G4endl;
    G4cout << "Is from primary lepton (tau or muon)    " << fFromPrimaryLepton << G4endl;
}

