#include "TrackingAction.hh"
#include "TrackInformation.hh"
#include "AnalysisManager.hh"

#include "G4TrackingManager.hh"
#include "G4Track.hh"

TrackingAction::TrackingAction() : G4UserTrackingAction() {;}

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
}

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  if (aTrack->GetParentID()==0) {
    AnalysisManager::GetInstance()->AddOnePrimaryTrack();
  }
  if (aTrack->GetParentID()==0) {
    if (aTrack->GetParticleDefinition()->GetPDGEncoding()==111) {
      // in case of pizero in the list of primary track
      // its decay products are also counted as primary particles, mostly 2 gammas
      G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
      if (secondaries) {
        size_t nSeco = secondaries->size();
        if (nSeco>0) {
          for (size_t i=0; i<nSeco; ++i) {
            if ((*secondaries)[i]->GetCreatorProcess()->GetProcessName()=="Decay") {
              TrackInformation* info =  new TrackInformation();
              info->SetTrackIsFromPrimaryPizero(1);
              (*secondaries)[i]->SetUserInformation(info);
              AnalysisManager::GetInstance()->AddOnePrimaryTrack();
            }
          }
        }
      } 
    }
  }
  if (aTrack->GetTrackID()==1 &&
      (abs(aTrack->GetParticleDefinition()->GetPDGEncoding())==15 ||
       abs(aTrack->GetParticleDefinition()->GetPDGEncoding())==13)) {
    // in case of the lepton decays, the decay products are counted as primary particles
    // * tau- decay (dominant)
    // * mu- decay
    G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
    if (secondaries) {
      size_t nSeco = secondaries->size();
      if (nSeco>0) {
        for (size_t i=0; i<nSeco; ++i) {
          if ((*secondaries)[i]->GetCreatorProcess()->GetProcessName()=="Decay") {
            TrackInformation* info =  new TrackInformation();
            info->SetTrackIsFromPrimaryLepton(1);
            (*secondaries)[i]->SetUserInformation(info);
            AnalysisManager::GetInstance()->AddOnePrimaryTrack();
          }
        }
      }
    }
  }
  if (aTrack->GetParentID()==1 && aTrack->GetCreatorProcess()->GetProcessName()=="Decay") {
    // in case of tau decay pizero
    // decay products of this pizero are also counted as primary particles, mostly 2 gammas
    if (aTrack->GetParticleDefinition()->GetPDGEncoding()==111) {
      G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
      if (secondaries) {
        size_t nSeco = secondaries->size();
        if (nSeco>0) {
          for (size_t i=0; i<nSeco; ++i) {
            if ((*secondaries)[i]->GetCreatorProcess()->GetProcessName()=="Decay") {
              TrackInformation* info =  new TrackInformation();
              info->SetTrackIsFromFSLPizero(1);
              (*secondaries)[i]->SetUserInformation(info);
              AnalysisManager::GetInstance()->AddOnePrimaryTrack();
            }
          }
        }
      }
    }
  }
  //TrackInformation* aTrackInfo = (TrackInformation*)(aTrack->GetUserInformation());
  //if (aTrackInfo) {
  //  if (aTrackInfo->IsTrackFromPrimaryTau() | aTrackInfo->IsTrackFromPrimaryPizero()) {
  //    std::cout<<aTrack->GetParentID()<<" "<<aTrack->GetParticleDefinition()->GetPDGEncoding()<<std::endl;
  //    aTrackInfo->Print();
  //  }
  //}
}
