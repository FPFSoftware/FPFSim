#include "TrackingAction.hh"
#include "TrackInformation.hh"
#include "AnalysisManager.hh"
#include "FPFTrajectory.hh"

#include "G4TrackingManager.hh"
#include "G4Track.hh"

TrackingAction::TrackingAction() : G4UserTrackingAction() {;}

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  // find out whether saving the full trajectory points or not
   bool storeTrajectoryPoints = AnalysisManager::GetInstance()->GetSaveTrajectories();

  // initialize our custom trajectory class
  auto *traj = new FPFTrajectory(aTrack, storeTrajectoryPoints);

  // hand it over to the manager
  fpTrackingManager->SetTrajectory(traj);
  fpTrackingManager->SetStoreTrajectory(true);
}

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  /// the following blocks store additional information into the track
  /// via the TrackInformation object (G4VTrackUserInformation)
  /// it is currently used in some FLARE SD volume to change how to save things
  /// TODO: consider for removal to simplify structure!
  if (aTrack->GetParentID()==0) 
  {
    if (aTrack->GetParticleDefinition()->GetPDGEncoding()==111) 
    {
      // in case of pizero in the list of primary track
      // its decay products are also counted as primary particles, mostly 2 gammas
      G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
      if (secondaries) 
      {
        size_t nSeco = secondaries->size();
        if (nSeco>0) 
        {
          for (size_t i=0; i<nSeco; ++i) 
          {
            if ((*secondaries)[i]->GetCreatorProcess()->GetProcessName()=="Decay") 
            {
              TrackInformation* info =  new TrackInformation();
              info->SetTrackIsFromPrimaryPizero(1);
              (*secondaries)[i]->SetUserInformation(info);
            }
          }
        }
      }
    }
  }

  if (aTrack->GetTrackID()==1 &&
      (abs(aTrack->GetParticleDefinition()->GetPDGEncoding())==15 ||
       abs(aTrack->GetParticleDefinition()->GetPDGEncoding())==13)) 
  {
    // in case of the lepton decays, the decay products are counted as primary particles
    // * tau- decay (dominant)
    // * mu- decay
    G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
    if (secondaries) 
    {
      size_t nSeco = secondaries->size();
      if (nSeco>0) 
      {
        for (size_t i=0; i<nSeco; ++i) 
        {
          if ((*secondaries)[i]->GetCreatorProcess()->GetProcessName()=="Decay") 
          {
            TrackInformation* info =  new TrackInformation();
            info->SetTrackIsFromPrimaryLepton(1);
            (*secondaries)[i]->SetUserInformation(info);
          }
        }
      }
    }
  }

  if (aTrack->GetParentID()==1 && aTrack->GetCreatorProcess()->GetProcessName()=="Decay") 
  {
    // in case of tau decay pizero
    // decay products of this pizero are also counted as primary particles, mostly 2 gammas
    if (aTrack->GetParticleDefinition()->GetPDGEncoding()==111) 
    {
      G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
      if (secondaries) 
      {
        size_t nSeco = secondaries->size();
        if (nSeco>0) 
        {
          for (size_t i=0; i<nSeco; ++i) 
          {
            if ((*secondaries)[i]->GetCreatorProcess()->GetProcessName()=="Decay") 
            {
              TrackInformation* info =  new TrackInformation();
              info->SetTrackIsFromFSLPizero(1);
              (*secondaries)[i]->SetUserInformation(info);
            }
          }
        }
      }
    }
  }
}
