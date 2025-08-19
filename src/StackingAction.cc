#include "StackingAction.hh"
#include "G4TrackStatus.hh"
#include "G4VProcess.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "AnalysisManager.hh"
#include "G4TrackingManager.hh"

StackingAction::StackingAction(RunAction* aRunAction, EventAction* aEventAction) :
  G4UserStackingAction(), fRunAction(aRunAction), fEventAction(aEventAction)
{;}

G4ClassificationOfNewTrack StackingAction::ClassifyNewTrack (const G4Track* aTrack)
{
  // for each track, build track ID to primary ancestor
  // primaries have themselves as ancestor
  // go up the tree until original primary for everything else
  G4int trackID = aTrack->GetTrackID();
  G4int parentID = aTrack->GetParentID();
  G4int ancestorID = -1;

  // Register primary tracks
  if (parentID==0) 
  {
    fEventAction->AddPrimaryTrack();
    ancestorID = trackID; // primary is its own ancestor
  }

  // Register only secondaries, i.e. tracks having ParentID > 0
  if (parentID > 0) 
  {
    fEventAction->AddSecondaryTrack();
    if (aTrack->GetParticleDefinition()->GetPDGEncoding() != 22)
    {
      fEventAction->AddSecondaryTrackNotGamma();
    }

    // if you are a secondary, ancestor is your parent's ancestor
    ancestorID = AnalysisManager::GetInstance()->GetTrackPrimaryAncestor(parentID);
  }

  // add track with its ancestor!!!
  AnalysisManager::GetInstance()->SetTrackPrimaryAncestor(trackID,ancestorID);

  // Do not affect track classification. Just return what would have
  // been returned by the base class
  return G4UserStackingAction::ClassifyNewTrack(aTrack);
}
