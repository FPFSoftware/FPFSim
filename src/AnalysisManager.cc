#include <vector>
#include <functional>
#include <iostream>
#include <string>
#include <map>
#include <iomanip>
#include <random>

#include <G4Event.hh>
#include <G4SDManager.hh>
#include <G4SystemOfUnits.hh>
#include <Randomize.hh>
#include <G4Poisson.hh>
#include <G4Trajectory.hh>

#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <THnSparse.h>
#include <TString.h>
#include <Math/ProbFunc.h>

#include "AnalysisManager.hh"
#include "geometry/GeometricalParameters.hh"
#include "LArBoxSD.hh"
#include "FASER2TrackerSD.hh"
#include "LArBoxHit.hh"
#include "PrimaryParticleInformation.hh"
#include "reco/PCAAnalysis3D.hh"
#include "reco/Cluster3D.hh"
#include "reco/LinearFit.hh"
#include "reco/ShowerLID.hh"
#include "reco/Barcode.hh"
#include "FPFParticle.hh"
#include "FPFNeutrino.hh"

//---------------------------------------------------------------------
//---------------------------------------------------------------------
// AnalysisManager "singleton" instance
// once initialized, can be used to point to AnalysisManager
// from anywhere else in the codebase
AnalysisManager *AnalysisManager::fInstance = 0;

AnalysisManager *AnalysisManager::GetInstance()
{
  if (!fInstance)
  {
    G4cout << "AnalysisManager: Re-initialization" << G4endl;
    fInstance = new AnalysisManager();
  }
  return fInstance;
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------
AnalysisManager::AnalysisManager()
{
  fFile = nullptr;
  fFilename = "test.root";
  fH5Filename = "test.h5";

  fMessenger = new AnalysisManagerMessenger(this);

  fEvt = nullptr;
  fTrk = nullptr;
  fPrim = nullptr;
  fFLArEHits = nullptr;
  fFLArEHCALHits = nullptr;
  fActsHitsTree = nullptr;
  fActsParticlesTree = nullptr;

  fSaveTrack = false;
  fSave3DEvd = false;
  fSave2DEvd = false;
  fSaveActs = true;
  fAddDiffusion = "false";
}

AnalysisManager::~AnalysisManager() {}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::bookEvtTree()
{
  fEvt = new TTree("event", "event info");
  fEvt->Branch("evtID", &evtID, "evtID/I");

  /* TODO: REVIEW, CLEAN-UP or MOVE TO FLARE-SPECIFIC TREE
  fEvt->Branch("FPFNeutrino"                , &neutrino, 96000, 0);
  fEvt->Branch("FPFParticle"                , &primaries, 96000, 0);
  fEvt->Branch("TotalDedxLongitudinal"     , TotalDedxLongitudinal      , "TotalDedxLongitudinal[3000]/D");
  fEvt->Branch("TrueTotalDedxLongitudinal" , TrueTotalDedxLongitudinal  , "TrueTotalDedxLongitudinal[3000]/D");
  fEvt->Branch("nPrimaryParticle"          , &nPrimaryParticle          , "nPrimaryParticle/I");
  fEvt->Branch("primaryParentPDG"          , primaryParentPDG           , "primaryParentPDG[nPrimaryParticle]/I");
  fEvt->Branch("primaryTrackLength"        , primaryTrackLength         , "primaryTrackLength[nPrimaryParticle]/D");
  fEvt->Branch("primaryTrackLengthInTPC"   , primaryTrackLengthInTPC    , "primaryTrackLengthInTPC[nPrimaryParticle]/D");
  fEvt->Branch("ProngEInDetector"          , ProngEInDetector           , "ProngEInDetector[nPrimaryParticle]/D");
  fEvt->Branch("ProngEInLAr"               , ProngEInLAr                , "ProngEInLAr[nPrimaryParticle]/D");
  fEvt->Branch("ProngEInHadCal"            , ProngEInHadCal             , "ProngEInHadCal[nPrimaryParticle]/D");
  fEvt->Branch("ProngEInMuonFinder"        , ProngEInMuonFinder         , "ProngEInMuonFinder[nPrimaryParticle]/D");
  fEvt->Branch("ProngEInMuonFinderLayer1X" , ProngEInMuonFinderLayer1X  , "ProngEInMuonFinderLayer1X[nPrimaryParticle]/D");
  fEvt->Branch("ProngEInMuonFinderLayer1Y" , ProngEInMuonFinderLayer1Y  , "ProngEInMuonFinderLayer1Y[nPrimaryParticle]/D");
  fEvt->Branch("ProngEInMuonFinderLayer2X" , ProngEInMuonFinderLayer2X  , "ProngEInMuonFinderLayer2X[nPrimaryParticle]/D");
  fEvt->Branch("ProngEInMuonFinderLayer2Y" , ProngEInMuonFinderLayer2Y  , "ProngEInMuonFinderLayer2Y[nPrimaryParticle]/D");
  fEvt->Branch("ProngAngleToBeamDir"       , ProngAngleToBeamDir        , "ProngAngleToBeamDir[nPrimaryParticle]/D");
  fEvt->Branch("ShowerLength"              , ShowerLength               , "ShowerLength[nPrimaryParticle]/D");
  fEvt->Branch("ShowerLengthInLAr"         , ShowerLengthInLAr          , "ShowerLengthInLAr[nPrimaryParticle]/D");
  fEvt->Branch("ShowerWidth"               , ShowerWidth                , "ShowerWidth[nPrimaryParticle]/D");
  fEvt->Branch("ShowerWidthInLAr"          , ShowerWidthInLAr           , "ShowerWidthInLAr[nPrimaryParticle]/D");
  fEvt->Branch("ProngAvgdEdx"              , ProngAvgdEdx               , "ProngAvgdEdx[nPrimaryParticle]/D");
  fEvt->Branch("ProngAvgdEdxInLAr"         , ProngAvgdEdxInLAr          , "ProngAvgdEdxInLAr[nPrimaryParticle]/D");
  fEvt->Branch("ProngdEdxAlongTrack"       , ProngdEdxAlongTrack        , "ProngdEdxAlongTrack[nPrimaryParticle][100]/D");
  fEvt->Branch("ProngdEdxTrackLength"      , ProngdEdxTrackLength       , "ProngdEdxTrackLength[nPrimaryParticle][100]/I");
  fEvt->Branch("dir_pol_x"                 , dir_pol_x                  , "dir_pol_x[nPrimaryParticle]/D");
  fEvt->Branch("dir_pol_y"                 , dir_pol_y                  , "dir_pol_y[nPrimaryParticle]/D");
  fEvt->Branch("dir_pol_z"                 , dir_pol_z                  , "dir_pol_z[nPrimaryParticle]/D");
  fEvt->Branch("dir_coc_x"                 , dir_coc_x                  , "dir_coc_x[nPrimaryParticle]/D");
  fEvt->Branch("dir_coc_y"                 , dir_coc_y                  , "dir_coc_y[nPrimaryParticle]/D");
  fEvt->Branch("dir_coc_z"                 , dir_coc_z                  , "dir_coc_z[nPrimaryParticle]/D");

  fEvt->Branch("sparseFractionMem"         , &sparseFractionMem         , "sparseFractionMem/D");
  fEvt->Branch("sparseFractionBins"        , &sparseFractionBins        , "sparseFractionBins/D");

  fEvt->Branch("edepInLAr"                 , &edepInLAr                 , "edepInLAr/D");
  fEvt->Branch("edepInHadCalX"             , &edepInHadCalX             , "edepInHadCalX/D");
  fEvt->Branch("edepInHadCalY"             , &edepInHadCalY             , "edepInHadCalY/D");
  fEvt->Branch("edepInMuonFinderX"         , &edepInMuonFinderX         , "edepInMuonFinderX/D");
  fEvt->Branch("edepInMuonFinderY"         , &edepInMuonFinderY         , "edepInMuonFinderY/D");
  fEvt->Branch("edepInHadAborb"            , &edepInHadAborb            , "edepInHadAborb/D");
  fEvt->Branch("edepInMuonFinderAbsorb"    , &edepInMuonFinderAbsorb    , "edepInMuonFinderAbsorb/D");
  fEvt->Branch("missCountedEnergy"         , &missCountedEnergy         , "missCountedEnergy/D");

  fEvt->Branch("nFromFSLParticles"         , &nFromFSLParticles         , "nFromFSLParticles/I");
  */
}

void AnalysisManager::bookPrimTree()
{
  fPrim = new TTree("primaries", "primaries info");
  fPrim->Branch("evtID", &evtID, "evtID/I");
  fPrim->Branch("vtxID", &primVtxID, "vtxID/I");
  fPrim->Branch("PDG", &primPDG, "PDG/I");
  fPrim->Branch("trackID", &primTrackID, "trackID/I");
  fPrim->Branch("barcode", &primParticleID, "bardcode/I");
  fPrim->Branch("mass", &primM, "mass/F");
  fPrim->Branch("charge", &primQ, "charge/F");
  fPrim->Branch("Vx", &primVx, "Vx/F"); // position
  fPrim->Branch("Vy", &primVy, "Vy/F");
  fPrim->Branch("Vz", &primVz, "Vz/F");
  fPrim->Branch("Vt", &primVt, "Vt/F");
  fPrim->Branch("Px", &primPx, "Px/F"); // momentum
  fPrim->Branch("Py", &primPy, "Py/F");
  fPrim->Branch("Pz", &primPz, "Pz/F");
  fPrim->Branch("E", &primE, "E/F");    // initial total energy
  fPrim->Branch("KE", &primKE, "KE/F"); // initial kinetic energy
  fPrim->Branch("Eta", &primEta, "Eta/F");
  fPrim->Branch("Phi", &primPhi, "Phi/F");
  fPrim->Branch("Pt", &primPt, "Pt/F");
  fPrim->Branch("P", &primP, "P/F");
}

void AnalysisManager::bookTrkTree()
{
  fTrk = new TTree("trajectories", "trajectories info");
  fTrk->Branch("evtID", &evtID, "evtID/I");
  fTrk->Branch("trackTID", &trackTID, "trackTID/I");
  fTrk->Branch("trackPID", &trackPID, "trackPID/I");
  fTrk->Branch("trackPDG", &trackPDG, "trackPDG/I");
  fTrk->Branch("trackKinE", &trackKinE, "trackKinE/D");
  fTrk->Branch("trackNPoints", &trackNPoints, "trackNPoints/I");
  fTrk->Branch("trackPointX", &trackPointX);
  fTrk->Branch("trackPointY", &trackPointY);
  fTrk->Branch("trackPointZ", &trackPointZ);
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::bookFLArETrees()
{
  // create subdirectory in file
  fFLArEDir = fFile->mkdir("flare","FLArE output",kTRUE);
  fFile->cd(fFLArEDir->GetName());

  // FLArE hits tree
  fFLArEHits = new TTree("flare_hits", "FLArE hits info");
  fFLArEHits->Branch("flareEvtID", &evtID, "flareEvtID/I");
  fFLArEHits->Branch("flareTrackID", &flareTrackID, "flareTrackID/I");
  fFLArEHits->Branch("flareBarcode", &flareParticleID, "flareParticleID/I");
  fFLArEHits->Branch("flareParentID", &flareParentID, "flareParentID/I");
  fFLArEHits->Branch("flarePDG", &flarePDG, "flarePDG/I");
  fFLArEHits->Branch("flareCopyNum", &flareCopyNum, "flareCopyNum/I");
  fFLArEHits->Branch("flareT", &flareT, "flareT/I");
  fFLArEHits->Branch("flareX", &flareX, "flareX/F"); // Pre-position
  fFLArEHits->Branch("flareY", &flareY, "flareY/F");
  fFLArEHits->Branch("flareZ", &flareZ, "flareZ/F");
  fFLArEHits->Branch("flarePx", &flarePx, "flarePx/F"); // momentum
  fFLArEHits->Branch("flarePy", &flarePy, "flarePy/F");
  fFLArEHits->Branch("flarePz", &flarePz, "flarePz/F");
  fFLArEHits->Branch("flareDeltaPx", &flareDeltaPx, "flareDeltaPx/F");
  fFLArEHits->Branch("flareDeltaPy", &flareDeltaPy, "flareDeltaPy/F");
  fFLArEHits->Branch("flareDeltaPz", &flareDeltaPz, "flareDeltaPz/F");
  fFLArEHits->Branch("flareEdep", &flareEdep, "flareEdep/F");

  // FLArE HCAL hits
  fFLArEHCALHits = new TTree("hcal_hits", "FLArE HCAL hits info");

  fFLArEHCALHits->Branch("hadEvtID", &evtID, "hadEvtID/I");
  fFLArEHCALHits->Branch("hadTrackID", &flareTrackID, "hadTrackID/I");
  fFLArEHCALHits->Branch("hadBarcode", &flareParticleID, "hadParticleID/I");
  fFLArEHCALHits->Branch("hadParentID", &flareParentID, "hadParentID/I");
  fFLArEHCALHits->Branch("hadPDG", &flarePDG, "hadPDG/I");
  fFLArEHCALHits->Branch("hadCopyNum", &flareCopyNum, "hadCopyNum/I");
  fFLArEHCALHits->Branch("hadT", &flareT, "hadT/I");
  fFLArEHCALHits->Branch("hadX", &flareX, "hadX/F"); // Pre-position
  fFLArEHCALHits->Branch("hadY", &flareY, "hadY/F");
  fFLArEHCALHits->Branch("hadZ", &flareZ, "hadZ/F");
  fFLArEHCALHits->Branch("hadPx", &flarePx, "hadPx/F"); // momentum
  fFLArEHCALHits->Branch("hadPy", &flarePy, "hadPy/F");
  fFLArEHCALHits->Branch("hadPz", &flarePz, "hadPz/F");
  fFLArEHCALHits->Branch("hadDeltaPx", &flareDeltaPx, "hadDeltaPx/F");
  fFLArEHCALHits->Branch("hadDeltaPy", &flareDeltaPy, "hadDeltaPy/F");
  fFLArEHCALHits->Branch("hadDeltaPz", &flareDeltaPz, "hadDeltaPz/F");
  fFLArEHCALHits->Branch("hadEdep", &flareEdep, "hadEdep/F");
  fFLArEHCALHits->Branch("hadIsZX", &flareIsZX, "hadIsZX/B");

  fFile->cd();
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::bookFASER2Trees()
{
  // create subdirectory in file
  fFASER2Dir = fFile->mkdir("faser2","FASER2 output",kTRUE);
  fFile->cd(fFASER2Dir->GetName());

  //* Acts Hits Tree [i == unsigned int; F == float; l == Long unsigned 64 int]
  fActsHitsTree = new TTree("hits", "ActsHitsTree");
  fActsHitsTree->Branch("event_id", &ActsHitsEventID, "event_id/i");
  fActsHitsTree->Branch("geometry_id", &ActsHitsGeometryID, "geometryid/l");
  fActsHitsTree->Branch("particle_id", &ActsHitsParticleID, "particle_id/l");
  fActsHitsTree->Branch("tx", &ActsHitsX, "tx/F");
  fActsHitsTree->Branch("ty", &ActsHitsY, "ty/F");
  fActsHitsTree->Branch("tz", &ActsHitsZ, "tz/F");
  fActsHitsTree->Branch("tt", &ActsHitsT, "tt/F");
  fActsHitsTree->Branch("tpx", &ActsHitsPx, "tpx/F");
  fActsHitsTree->Branch("tpy", &ActsHitsPy, "tpy/F");
  fActsHitsTree->Branch("tpz", &ActsHitsPz, "tpz/F");
  fActsHitsTree->Branch("te", &ActsHitsE, "tpe/F");
  fActsHitsTree->Branch("deltapx", &ActsHitsDeltaPx, "deltapx/F");
  fActsHitsTree->Branch("deltapy", &ActsHitsDeltaPy, "deltapy/F");
  fActsHitsTree->Branch("deltapz", &ActsHitsDeltaPz, "deltapz/F");
  fActsHitsTree->Branch("deltae", &ActsHitsDeltaE, "deltae/F");
  fActsHitsTree->Branch("index", &ActsHitsIndex, "index/I");
  fActsHitsTree->Branch("volume_id", &ActsHitsVolumeID, "volume_id/i");
  fActsHitsTree->Branch("boundary_id", &ActsHitsBoundaryID, "boundary_id/i");
  fActsHitsTree->Branch("layer_id", &ActsHitsLayerID, "layer_id/i");
  fActsHitsTree->Branch("approach_id", &ActsHitsApproachID, "approach_id/i");
  fActsHitsTree->Branch("sensitive_id", &ActsHitsSensitiveID, "sensitive_id/i");

  //* Acts truth particle tree
  fActsParticlesTree = new TTree("particles", "ActsParticlesTree");
  fActsParticlesTree->Branch("event_id", &ActsHitsEventID, "event_id/i");
  fActsParticlesTree->Branch("particle_id", &ActsParticlesParticleId);
  fActsParticlesTree->Branch("particle_type", &ActsParticlesParticleType);
  fActsParticlesTree->Branch("process", &ActsParticlesProcess);
  fActsParticlesTree->Branch("vx", &ActsParticlesVx);
  fActsParticlesTree->Branch("vy", &ActsParticlesVy);
  fActsParticlesTree->Branch("vz", &ActsParticlesVz);
  fActsParticlesTree->Branch("vt", &ActsParticlesVt);
  fActsParticlesTree->Branch("px", &ActsParticlesPx);
  fActsParticlesTree->Branch("py", &ActsParticlesPy);
  fActsParticlesTree->Branch("pz", &ActsParticlesPz);
  fActsParticlesTree->Branch("m", &ActsParticlesM);
  fActsParticlesTree->Branch("q", &ActsParticlesQ);
  fActsParticlesTree->Branch("eta", &ActsParticlesEta);
  fActsParticlesTree->Branch("phi", &ActsParticlesPhi);
  fActsParticlesTree->Branch("pt", &ActsParticlesPt);
  fActsParticlesTree->Branch("p", &ActsParticlesP);
  fActsParticlesTree->Branch("vertex_primary", &ActsParticlesVertexPrimary);
  fActsParticlesTree->Branch("vertex_secondary", &ActsParticlesVertexSecondary);
  fActsParticlesTree->Branch("particle", &ActsParticlesParticle);
  fActsParticlesTree->Branch("generation", &ActsParticlesGeneration);
  fActsParticlesTree->Branch("sub_particle", &ActsParticlesSubParticle);
  fActsParticlesTree->Branch("e_loss", &ActsParticlesELoss);
  fActsParticlesTree->Branch("total_x0", &ActsParticlesPathInX0);
  fActsParticlesTree->Branch("total_l0", &ActsParticlesPathInL0);
  fActsParticlesTree->Branch("number_of_hits", &ActsParticlesNumberOfHits);
  fActsParticlesTree->Branch("outcome", &ActsParticlesOutcome);

  fFile->cd();
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::BeginOfRun()
{
  G4cout << "Run has been started, preparing output" << G4endl;

  if (fFile)
    delete fFile;

  // Preparing output files
  fFile = new TFile(fFilename.c_str(), "RECREATE");
  fH5Filename = fFilename;
  if (fH5Filename.find(".root") != std::string::npos)
  {
    const size_t pos = fH5Filename.find(".root");
    fH5Filename.resize(pos);
  }
  fH5Filename += ".h5";
  fH5file = hep_hpc::hdf5::File(fH5Filename, H5F_ACC_TRUNC);

  // Booking common output trees
  bookEvtTree();
  bookPrimTree();
  if (fSaveTrack) bookTrkTree();

  fSDNamelist = GeometricalParameters::Get()->GetSDNamelist();
  G4cout << "Number of SDs : " << fSDNamelist.size() << G4endl;

  // Loop through active sentive detectors
  // book trees as necessary
  for (auto sdname : fSDNamelist)
  {
    G4cout << sdname.first << " " << sdname.second << G4endl;

    // if FLArE is enabled, book trees
    if ( sdname.second.find("FLArE") != std::string::npos )
    {
      fFlareSDs.push_back(sdname.first);
      bookFLArETrees();
    }

    // if FASER2 is enabled, book ACTS trees
    // but give the option to disable output if required
    else if (fSaveActs && sdname.second.find("FASER2") != std::string::npos )
    {
      fFaser2SDs.push_back(sdname.first);
      bookFASER2Trees();
    }
  }
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::EndOfRun()
{
  // save common trees at the top of the output file
  fFile->cd();
  fEvt->Write();
  fPrim->Write();
  if (fSaveTrack) fTrk->Write();

  // save detector-specif trees in their directories
  if (fFlareSDs.size()>0)
  {
    fFile->cd(fFLArEDir->GetName());
    fFLArEHits->Write();
    fFLArEHCALHits->Write();
    fFile->cd(); // go back to top
  }
  if (fFaser2SDs.size()>0)
  {
    fFile->cd(fFASER2Dir->GetName());
    fActsHitsTree->Write();
    fActsParticlesTree->Write();
    fFile->cd(); // go back to top
  }

  fFile->Close();
  fH5file.close();
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::BeginOfEvent()
{
  // TODO: REVIEW, CLEAN-UP UNUSED VARIABLES
  // TODO: change arrays into std::vectors
  // reset vectors that need to be cleared for a new event
  // only reset arrays or vectors, no need for other defaults??

  //neutrino = FPFNeutrino();
  primaries.clear();
  primaryIDs.clear();

  /*
  for (G4int j = 0; j < 3000; ++j)
  {
    TotalDedxLongitudinal[j] = 0;
    TrueTotalDedxLongitudinal[j] = 0;
  }
  for (G4int i = 0; i < 1000; ++i)
  {
    primaryParentPDG[i] = 0;
    primaryTrackLength[i] = 0;
    primaryTrackLengthInTPC[i] = 0;
    ProngEInDetector[i] = 0;
    ProngEInLAr[i] = 0;
    ProngEInHadCal[i] = 0;
    ProngEInMuonFinder[i] = 0;
    ProngEInMuonFinderLayer1X[i] = 0;
    ProngEInMuonFinderLayer1Y[i] = 0;
    ProngEInMuonFinderLayer2X[i] = 0;
    ProngEInMuonFinderLayer2Y[i] = 0;
    ProngAngleToBeamDir[i] = -1;
    ShowerLength[i] = -1;
    ShowerLengthInLAr[i] = -1;
    ShowerWidth[i] = 0;
    ShowerWidthInLAr[i] = 0;
    ProngAvgdEdx[i] = -1;
    ProngAvgdEdxInLAr[i] = -1;
    for (G4int j = 0; j < 100; ++j)
    {
      ProngdEdxAlongTrack[i][j] = 0;
      ProngdEdxTrackLength[i][j] = -1;
    }
    dir_pol_x[i] = -999;
    dir_pol_y[i] = -999;
    dir_pol_z[i] = -999;
    dir_coc_x[i] = -999;
    dir_coc_y[i] = -999;
    dir_coc_z[i] = -999;
    fromFSLParticlePDG[i] = 0;
  }

  allTracksPTPair.clear();
  trackClusters.clear();
  tracksFromFSL.clear();
  tracksFromFSLSecondary.clear();
  tracksFromFSPizeroSecondary.clear();
  tracksFromFSLDecayPizeroSecondary.clear();
  */
  trackPointX.clear();
  trackPointY.clear();
  trackPointZ.clear();

  ActsParticlesParticleId.clear();
  ActsParticlesParticleType.clear();
  ActsParticlesProcess.clear();
  ActsParticlesVx.clear();
  ActsParticlesVy.clear();
  ActsParticlesVz.clear();
  ActsParticlesVt.clear();
  ActsParticlesPx.clear();
  ActsParticlesPy.clear();
  ActsParticlesPz.clear();
  ActsParticlesM.clear();
  ActsParticlesQ.clear();
  ActsParticlesEta.clear();
  ActsParticlesPhi.clear();
  ActsParticlesPt.clear();
  ActsParticlesP.clear();
  ActsParticlesVertexPrimary.clear();
  ActsParticlesVertexSecondary.clear();
  ActsParticlesParticle.clear();
  ActsParticlesGeneration.clear();
  ActsParticlesSubParticle.clear();
  ActsParticlesELoss.clear();
  ActsParticlesPathInX0.clear();
  ActsParticlesPathInL0.clear();
  ActsParticlesNumberOfHits.clear();
  ActsParticlesOutcome.clear();
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::EndOfEvent(const G4Event *event)
{

  /// evtID
  evtID = event->GetEventID();

  // FILL EVENT TREE
  // nothing else except evtID for now...
  // TODO: need to find good way to propagate generator-level info
  fEvt->Fill();

  //-----------------------------------------------------------

  // FILL PRIMARIES/TRAJECTORIES TREE
  FillPrimariesTree(event);
  if(fSaveTrack) FillTrajectoriesTree(event);

  //-----------------------------------------------------------

  // Get the hit collections
  // If there is no hit collection, there is nothing to be done
  fHCofEvent = event->GetHCofThisEvent();
  if (!fHCofEvent)
  {
    G4cout << "No hits recorded in any sensitive volume --> nothing to save!" << G4endl;
    return;
  }

  //-----------------------------------------------------------

  // FILL DETECTOR HITS TREES
  if( fFlareSDs.size() > 0 )  FillFLArEOutput();
  if( fFaser2SDs.size() > 0 ) FillFASER2Output();

  /*
  // update number of primary particles after FillPrimaryTruthTree
  // including decay products from primary tau and pizero
  nPrimaryParticle = countPrimaryParticle; // TODO: they're equivalent, should only store one of them
  nFromFSLParticles = tracksFromFSLSecondary.size();
  nFromFSPizeroParticles = tracksFromFSPizeroSecondary.size();
  nFromFSLDecayPizeroParticles = tracksFromFSLDecayPizeroSecondary.size();
  // G4cout<<"passed9"<<G4endl;

  // find all the tracks originate from the final state lepton, include FSL itself (TID=1)
  // should only work with neutrino interaction generator
  // exception to single particle generator: tau, mu
  /* if (neutrino.NuPDG()!=0 || abs(primaries[0].PDG())==15 || abs(primaries[0].PDG())==13) {
  //this part causing errors, so commenting out */
  /*   tracksFromFSL.insert(1); */
  /*   for (auto x : allTracksPTPair) { */
  /*     if (tracksFromFSL.find(x.first) != tracksFromFSL.end()) { */
  /*       tracksFromFSL.insert(x.second); */
  /*     } */
  /*   } */
  /* } */
  // G4cout<<"passed10"<<G4endl;

  // tracksFromFSL includes all the tracks orginating from the fsl
  // tracksFromFSLSecondary only inclues the tracks directly decayed from the fsl
  /* std::cout<<"Recorded tracks       : "<<allTracksPTPair.size()<<std::endl; */
  /* std::cout<<"Tracks from FSL       : "<<tracksFromFSL.size()<<std::endl; */
  /* std::cout<<"Tracks from FSL (2nd) : "<<tracksFromFSLSecondary.size()<<std::endl; */
  /* std::cout<<"number of primary particles : "<<nPrimaryParticle */
  /*   <<" , in which number of particles from fsl : "<<nFromFSLParticles<<std::endl; */
  // std::cout<<"Test nTestNPrimaryTrack : "<<nTestNPrimaryTrack<<std::endl;
  //  return;
  //- cluster all tracks to primary particles
  //- mark the index of the final state lepton from the neutrino interaction
  /* trackClusters.resize(nPrimaryParticle); */
  /* for (int iPrim= 0; iPrim< nPrimaryParticle; ++iPrim) { */
  /*   trackClusters[iPrim].insert(primaries[iPrim].TID()); */
  /*   if (primaries[iPrim].TID()==1 && neutrino.NuPDG()!=0) { */
  /*     fPrimIdxFSL = iPrim; */
  /*   } */
  /* } */
  /* if (fPrimIdxFSL>=0) primaries[fPrimIdxFSL].SetProngType(0); */

  /* for (auto x : allTracksPTPair) { */
  /*   // if this track is the fsl (TID=1) and it decays (nFromFSLParticles>0), */
  /*   // then it forms a single cluster by itself, this is mainly for studying the tau decay. */
  /*   if ((x.second==1) && (nFromFSLParticles>0)) continue; */
  /*   // if this track is the decay product of the fsl, it should already been added to the trackClusters */
  /*   if ((x.first==1) && (nFromFSLParticles>0) && (tracksFromFSLSecondary.find(x.second) != tracksFromFSLSecondary.end())) continue; */
  /*   // if this is the decay product of the final state pizero, it should already been added to the trackClusters */
  /*   if ((tracksFromFSPizeroSecondary.find(x.second) != tracksFromFSPizeroSecondary.end())) continue; */
  /*   // if this is the decay product of the tau decay pizero, it should already been added to the trackClusters */
  /*   if ((tracksFromFSLDecayPizeroSecondary.find(x.second) != tracksFromFSLDecayPizeroSecondary.end())) continue; */
  /*   // add the track to the corresponding cluster if its parent is in the cluster. */
  /*   // one track can have only one parent, break the loop once its parent is found. */
  /*   for (int iPrim= 0; iPrim< nPrimaryParticle; ++iPrim) { */
  /*     if (trackClusters[iPrim].find(x.first) != trackClusters[iPrim].end()) { */
  /*       trackClusters[iPrim].insert(x.second); */
  /*       break; */
  /*     } */
  /*   } */
  /* } */

  /* if (GeometricalParameters::Get()->GetAddFLArE()) */
  /* { */
  /*   G4cout << "Adding FLArE pixel map..." << G4endl; */
  /*   const Double_t res_tpc[3] = {1, 5, 5}; // mm */
  /*   if (neutrino.NuPDG()!=0) { */
  /*     pm3D = new PixelMap3D(evtID, nPrimaryParticle, neutrino.NuPDG(), res_tpc); */
  /*   } else { */
  /*     pm3D = new PixelMap3D(evtID, nPrimaryParticle, primaries[0].PDG(), res_tpc); */
  /*   } */
  /*   // boundary in global coord. */
  /*   pm3D->SetPMBoundary(GeometricalParameters::Get()->GetFLArEPosition()/mm - */
  /*                         GeometricalParameters::Get()->GetTPCSizeXYZ()/mm/2, */
  /*                         GeometricalParameters::Get()->GetFLArEPosition()/mm + */
  /*                         GeometricalParameters::Get()->GetTPCSizeXYZ()/mm/2); */
  /*   pm3D->InitializePM(); */

  /*   /// FillTrueEdep must run after FillPrimaryTruthTree, */
  /*   /// otherwise tracksFromFSL and tracksFromFSLSecondary are invalid */
  /*   /// Pixel map is also filled here */
  /*   for (auto sdname : SDNamelist) { */
  /*     FillTrueEdep(sdname.first, sdname.second); */
  /*   } */

  /*   if (m_save2DEvd) pm3D->Write2DPMToFile(thefile); */

  /*   pm3D->Process3DPM(fH5file, neutrino, m_save3DEvd); */
  /*   sparseFractionMem = pm3D->GetSparseFractionMem(); */
  /*   sparseFractionBins = pm3D->GetSparseFractionBins(); */

  /*

    /* // FillPseudoRecoVar must run after FillTrueEdep, otherwise some of the variables won't be filled */
  /* FillPseudoRecoVar(); */

  /* delete pm3D; */
  /* } */
  // G4cout<<"passed3"<<G4endl;

  
  /* evt->Fill(); */

  /* G4cout << "Total number of recorded hits : " << nHits << std::endl; */
  /* if(m_saveTrack) G4cout << "Total number of recorded track: " << count_tracks << std::endl; */

  /* for (int iPrim= 0; iPrim< nPrimaryParticle; ++iPrim) { */
  /*   trackClusters[iPrim].clear(); */
  /* } */
  /* trackClusters.clear(); */
  /* trackClusters.shrink_to_fit(); */
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::FillPrimariesTree(const G4Event *event)
{
  nPrimaryVertex = event->GetNumberOfPrimaryVertex();
  G4cout << "\nNumber of primary vertices  : " << nPrimaryVertex << G4endl;
  
  /// loop over the vertices, and then over primary particles,
  /// neutrino truth info from event generator.
  for (G4int ivtx = 0; ivtx < event->GetNumberOfPrimaryVertex(); ++ivtx)
  {
    G4cout << "=== Vertex " << ivtx+1 << " of " << nPrimaryVertex << " -> " 
           << event->GetPrimaryVertex(ivtx)->GetNumberOfParticle() << " primaries ===" << G4endl;
    for (G4int ipp = 0; ipp < event->GetPrimaryVertex(ivtx)->GetNumberOfParticle(); ++ipp)
    {
      G4PrimaryParticle *primary_particle = event->GetPrimaryVertex(ivtx)->GetPrimary(ipp);
      if (primary_particle)
      {
        PrimaryParticleInformation *primary_particle_info = dynamic_cast<PrimaryParticleInformation *>(primary_particle->GetUserInformation());
        primary_particle_info->Print();

        primVtxID = ivtx;
        primTrackID = primary_particle_info->GetPartID(); // confirm matches track's id later?

        auto particleId = ActsFatras::Barcode();
        particleId.setVertexPrimary(ivtx);
        particleId.setGeneration(0);
        particleId.setSubParticle(0);
        particleId.setParticle(primTrackID - 1);

        primParticleID = particleId.value();
        primPDG = primary_particle_info->GetPDG();
        primVx = primary_particle_info->GetVertexMC().x();
        primVy = primary_particle_info->GetVertexMC().y();
        primVz = primary_particle_info->GetVertexMC().z();
        primVt = 0;
        primPx = primary_particle_info->GetMomentumMC().x();
        primPy = primary_particle_info->GetMomentumMC().y();
        primPz = primary_particle_info->GetMomentumMC().z();
        primM = primary_particle_info->GetMass()/MeV;
        primQ = primary_particle_info->GetCharge();

        G4double energy = GetTotalEnergy(primPx, primPy, primPz, primM);
        TLorentzVector p4(primPx,primPy,primPz,energy);
        primEta = p4.Eta();
        primPhi = p4.Phi();
        primPt = p4.Pt();
        primP = p4.P();
        primE = energy;
        primKE = energy - primM;

        // store a copy as a FPFParticle for further processing
        primaryIDs.push_back(primTrackID); //store to avoid duplicates
        primaries.push_back(FPFParticle(primPDG, 0, 
		                        primTrackID, primaryIDs.size()-1, 1,
		                        primM,
                            primVx, primVy, primVz, primVt,
                            primPx, primPy, primPz,energy));
        fPrim->Fill();
      }
    }
  }

  G4cout << "\nNumber of primaries  : " << primaryIDs.size() << G4endl;
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::FillTrajectoriesTree(const G4Event* event)
{
  int count_tracks = 0;

  G4cout << "==== Saving track information to tree ====" << G4endl; 
  auto trajectoryContainer = event->GetTrajectoryContainer(); 
  if (!trajectoryContainer)
  {
    G4cout << "No tracks found: did you enable their storage with '/tracking/storeTrajectory 1'?" << G4endl;
    return;
  }

  for (size_t i = 0; i < trajectoryContainer->entries(); ++i) 
  { 
    auto trajectory = static_cast<G4Trajectory*>((*trajectoryContainer)[i]); 
    trackTID = trajectory->GetTrackID();
    trackPID = trajectory->GetParentID();
    trackPDG = trajectory->GetPDGEncoding(); 
    trackKinE = trajectory->GetInitialKineticEnergy(); 
    trackNPoints = trajectory->GetPointEntries(); 
    count_tracks++; 
    for (size_t j = 0; j < trackNPoints; ++j) 
    { 
      G4ThreeVector pos = trajectory->GetPoint(j)->GetPosition(); 
      trackPointX.push_back( pos.x() );
      trackPointY.push_back( pos.y() );
      trackPointZ.push_back( pos.z() );
    }
    fTrk->Fill();
    trackPointX.clear(); 
    trackPointY.clear();
    trackPointZ.clear();
  }
  G4cout << "Total number of recorded track: " << count_tracks << std::endl;
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::FillFLArEOutput()
{
  G4cout << "===== Filling FLArE output trees =====" << G4endl;

  // loop over the detected FLArE sensitive volumes
  // for now, everything is lar_box
  int nHits = 0;
  for(const int sdId : fFlareSDs )
  {
    //  Get and cast hit collection with LArBoxHits
    std::string sdName = fSDNamelist.at(sdId);
    auto hitCollection = dynamic_cast<LArBoxHitsCollection *>(fHCofEvent->GetHC(sdId));
    if (!hitCollection)
    {
      G4cout << "No hits recorded by " << sdName << G4endl;
      continue;
    }
  
    for (auto hit : *hitCollection->GetVector())
    {
      nHits++;
      flareTrackID = hit->GetTID();
      flareParentID = hit->GetPID();
      flarePDG = hit->GetPDG();
      flareCopyNum = hit->GetCopyNum();
      flareT = hit->GetTime();
      flareX = hit->GetPreStepPosition().x();
      flareY = hit->GetPreStepPosition().y();
      flareZ = hit->GetPreStepPosition().z();
      flarePx = hit->GetInitMomentum().x();
      flarePy = hit->GetInitMomentum().y();
      flarePz = hit->GetInitMomentum().z();
      flareDeltaPx = hit->GetDeltaMom().x();
      flareDeltaPy = hit->GetDeltaMom().y();
      flareDeltaPz = hit->GetDeltaMom().z();
      flareEdep = hit->GetEdep();

      auto particleId = ActsFatras::Barcode();
      particleId.setVertexPrimary(1); // fix this value
      particleId.setGeneration(flareParentID);
      particleId.setSubParticle(0);
      particleId.setParticle(flareTrackID);
      flareParticleID = particleId.value();

      if (sdName == "FLArEBoxSD/lar_box") fFLArEHits->Fill();
      else if (sdName == "FLArEHadCalXSD/lar_box" || 
               sdName == "FLArEMuonFinderXSD/lar_box" ||
               sdName == "FLArEBabyMINDHorBarSD/lar_box" )
      {
        // specify this is an ZX hit
        flareIsZX = true;
        fFLArEHCALHits->Fill();
      }
      else if (sdName == "FLArEHadCalYSD/lar_box" || 
               sdName == "FLArEMuonFinderYSD/lar_box" ||
               sdName == "FLArEBabyMINDVerBarSD/lar_box" ||
               sdName == "FLArEHadAbsorbSD/lar_box" ||
               sdName == "FLArEMuonFinderAbsorbSD/lar_box" )
      {
        flareIsZX = false;
        fFLArEHCALHits->Fill();
      }
    }
  }

  G4cout << "Total FLArE + FLArE HCAL recorded hits: " << nHits << G4endl;
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::FillFASER2Output()
{
  G4cout << "==== Filling FASER2 output trees ====" << G4endl;

  // loop over the detected FASER2 sensitive volumes
  int nHits = 0;
  for(const int sdId : fFaser2SDs )
  {
    //  Get and cast hit collection with LArBoxHits
    std::string sdName = fSDNamelist.at(sdId);
    auto hitCollection = dynamic_cast<FASER2TrackerHitsCollection *>(fHCofEvent->GetHC(sdId));
    if (!hitCollection)
    {
      G4cout << "No hits recorded by " << sdName << G4endl;
      continue;
    }
  
    std::map<G4int, G4int> sub_part_map{};
    for (auto hit : *hitCollection->GetVector())
    {
      if (hit->GetCharge() == 0)
        continue; // skip neutral particles, they don't hit

      /*
      * A note on the ActsHitsGeometryID variable
       This variable in Acts keeps track of an Acts::GeometryIdentifier. This is essentially a long unsigned int, the bits of which are used to
       look up the the volume/layer/boundary/sensitive indices of a piece of geometry. In principle it should be possible to assign this variable
       here in GEANT4 but I don't understand the Acts code well enough to do it without adding Acts as a dependancy to this codebase.
       As a result I set `geometry_id` to zero and give the the resposibility of assigning this variable to the user during the reading of the `hits` tree.
      */

      nHits++;
      ActsHitsEventID = evtID;
      ActsHitsGeometryID = 0;

      int hitID = hit->GetTrackID();
      int nPrimaries = ActsParticlesParticleId.size();

      auto particleId = ActsFatras::Barcode();
      particleId.setVertexPrimary(1);
      particleId.setVertexSecondary(0);
      particleId.setParticle(hit->GetTrackID() - 1); // The track ID is the primary particle index plus one
      particleId.setGeneration(hit->GetParentID());

      sub_part_map.try_emplace(hit->GetTrackID() - 1, sub_part_map.size());

      // This is a fudge - assumes that that the secondary particles are always sub-particles of the primary particle
      particleId.setSubParticle(hit->GetParentID() == 0 ? 0 : sub_part_map[hit->GetTrackID() - 1]);
      ActsHitsParticleID = particleId.value();

      ActsHitsX = hit->GetX();
      ActsHitsY = hit->GetY();
      ActsHitsZ = hit->GetZ();
      ActsHitsT = hit->GetT();
      ActsHitsPx = hit->GetPx();
      ActsHitsPy = hit->GetPy();
      ActsHitsPz = hit->GetPz();
      ActsHitsE = hit->GetEnergy();
      ActsHitsDeltaPx = hit->GetDeltaPx();
      ActsHitsDeltaPy = hit->GetDeltaPy();
      ActsHitsDeltaPz = hit->GetDeltaPz();
      ActsHitsDeltaE = hit->GetDeltaE();
      ActsHitsIndex = hit->GetCopyNumSensor(); // index of layer: 0, 1, 2, ...

      // These variables I'm not 100% sure about. I reverse engineered them by matching them to how they're set when writing the hits from the particle gun in Acts
      // In principle with the right headers from Acts we could construct the geometry ID value here
      ActsHitsVolumeID = 1;
      ActsHitsBoundaryID = 0;
      ActsHitsLayerID = (hit->GetCopyNumSensor() + 1) * 2; // Acts specfic layer ID, goes 2, 4, 6, ...
      ActsHitsApproachID = 0;
      ActsHitsSensitiveID = 1;
      fActsHitsTree->Fill();

      // Now fill the Acts particles tree
      bool isDuplicate = false;
      for (const auto &id : ActsParticlesParticleId)
      {
        if (id == particleId.value())
        {
          isDuplicate = true;
        }
      }
      if (isDuplicate) continue; // Skip this particle if it's already been added

      ActsParticlesParticleId.push_back(particleId.value());
      ActsParticlesParticleType.push_back(hit->GetPDGID());
      ActsParticlesProcess.push_back(0);
      ActsParticlesVx.push_back(hit->GetTrackVertex().x());
      ActsParticlesVy.push_back(hit->GetTrackVertex().y());
      ActsParticlesVz.push_back(hit->GetTrackVertex().z());
      ActsParticlesVt.push_back(0);
      ActsParticlesPx.push_back(hit->GetTrackP4().px());
      ActsParticlesPy.push_back(hit->GetTrackP4().py());
      ActsParticlesPz.push_back(hit->GetTrackP4().pz());
      ActsParticlesM.push_back(hit->GetTrackP4().m());
      ActsParticlesQ.push_back(hit->GetCharge());

      ActsParticlesEta.push_back(hit->GetTrackP4().eta());
      ActsParticlesPhi.push_back(hit->GetTrackP4().phi());
      ActsParticlesPt.push_back(pow(pow(hit->GetTrackP4().px(), 2) + pow(hit->GetTrackP4().py(), 2), 0.5));
      ActsParticlesP.push_back(pow(pow(hit->GetTrackP4().px(), 2) + pow(hit->GetTrackP4().py(), 2) + pow(hit->GetTrackP4().pz(), 2), 0.5));
      ActsParticlesVertexPrimary.push_back(hit->GetIsPrimaryTrack());     //? These variables need to be filled, but are unused by Acts
      ActsParticlesVertexSecondary.push_back(hit->GetIsSecondaryTrack()); //? These variables need to be filled, but are unused by Acts
      ActsParticlesParticle.push_back(1);                                 //? These variables need to be filled, but are unused by Acts
      ActsParticlesGeneration.push_back(0);                               //? These variables need to be filled, but are unused by Acts
      ActsParticlesSubParticle.push_back(0);                              //? These variables need to be filled, but are unused by Acts
      ActsParticlesELoss.push_back(0);                                    //? These variables need to be filled, but are unused by Acts
      ActsParticlesPathInX0.push_back(0);                                 //? These variables need to be filled, but are unused by Acts
      ActsParticlesPathInL0.push_back(0);                                 //? These variables need to be filled, but are unused by Acts
      ActsParticlesNumberOfHits.push_back(0);                             //? These variables need to be filled, but are unused by Acts
      ActsParticlesOutcome.push_back(0);                                  //? These variables need to be filled, but are unused by Acts
    } // end of loop over hits
    fActsParticlesTree->Fill();
  }

  G4cout << "Total FASER2 recorded hits: " << nHits << G4endl;
}

  // Get and cast hit collection with LArBoxHits
  // LArBoxHitsCollection* hitCollection = dynamic_cast<LArBoxHitsCollection*>(hcofEvent->GetHC(sdId));
  // if (hitCollection) {
  //   for (auto hit: *hitCollection->GetVector()) {
  //     nHits++;

  //     double pre_x  = hit->GetPreStepPosition().x();
  //     double pre_y  = hit->GetPreStepPosition().y();
  //     double pre_z  = hit->GetPreStepPosition().z();
  //     double post_x = hit->GetPostStepPosition().x();
  //     double post_y = hit->GetPostStepPosition().y();
  //     double post_z = hit->GetPostStepPosition().z();

  //     // energy deposition in different volumes of the detector
  //     if (sdName == "lArBoxSD/lar_box")
  //       edepInLAr += hit->GetEdep();
  //     else if (sdName == "HadCalXSD/lar_box")
  //       edepInHadCalX += hit->GetEdep();
  //     else if (sdName == "HadCalYSD/lar_box")
  //       edepInHadCalY += hit->GetEdep();
  //     else if (sdName == "MuonFinderXSD/lar_box")
  //       edepInMuonFinderX += hit->GetEdep();
  //     else if (sdName == "MuonFinderYSD/lar_box")
  //       edepInMuonFinderY += hit->GetEdep();
  //     else if (sdName == "HadAbsorbSD/lar_box")
  //       edepInHadAborb += hit->GetEdep();
  //     else if (sdName == "MuonFinderAbsorbSD/lar_box")
  //       edepInMuonFinderAbsorb += hit->GetEdep();

  //     //allTracksPTPair.insert(std::make_pair(hit->GetPID(), hit->GetTID()));

  //     // stable final state particles in GENIE, primary particles in Geant4
  //     if (hit->GetCreatorProcess()=="PrimaryParticle") { // i.e. PID==0
  //       if ( std::find(primaryIDs.begin(), primaryIDs.end(), hit->GetTID()) == primaryIDs.end() ) {
  //         // the following line excludes final state lepton tau from the primary particle list
  //         //if (abs(nuPDG)==16 && abs(nuFSLPDG)==15 && abs(hit->GetParticle()==15)) continue;
  //         countPrimaryParticle++;
  //   primaryIDs.push_back(hit->GetTID());
  //         primaries.push_back(FPFParticle(hit->GetParticle(),
  //               hit->GetPID(), hit->GetTID(), countPrimaryParticle-1, 1, hit->GetParticleMass(),
  //               hit->GetTrackVertex().x(), hit->GetTrackVertex().y(), hit->GetTrackVertex().z(), 0,
  //               hit->GetInitMomentum().x(), hit->GetInitMomentum().y(), hit->GetInitMomentum().z(),
  //               GetTotalEnergy(hit->GetInitMomentum().x(), hit->GetInitMomentum().y(),
  //                 hit->GetInitMomentum().z(), hit->GetParticleMass())));
  //       }
  //     }
  //     // in case of the fsl decay, the decay products are counted as primary particles
  //     // * tau- decay (dominant)
  //     // * mu- decay
  //     //if (hit->GetPID()==1 && hit->GetCreatorProcess()=="Decay") {
  //     if (hit->GetIsTrackFromPrimaryLepton()) {
  //       tracksFromFSLSecondary.insert(hit->GetTID());
  //       if (std::find(primaryIDs.begin(), primaryIDs.end(), hit->GetTID()) == primaryIDs.end()) {
  //         countPrimaryParticle++;
  //   primaryIDs.push_back(hit->GetTID());
  //         primaries.push_back(FPFParticle(hit->GetParticle(),
  //               hit->GetPID(), hit->GetTID(), countPrimaryParticle-1, 2, hit->GetParticleMass(),
  //               hit->GetTrackVertex().x(), hit->GetTrackVertex().y(), hit->GetTrackVertex().z(), 0,
  //               hit->GetInitMomentum().x(), hit->GetInitMomentum().y(), hit->GetInitMomentum().z(),
  //               GetTotalEnergy(hit->GetInitMomentum().x(), hit->GetInitMomentum().y(),
  //                 hit->GetInitMomentum().z(), hit->GetParticleMass())));
  //       }
  //     }
  //     // in case of pizero in the list of primary track
  //     // its decay products are also counted as primary particles, mostly 2 gammas
  //     if (hit->GetIsTrackFromPrimaryPizero()) {
  //       tracksFromFSPizeroSecondary.insert(hit->GetTID());
  //       if (std::find(primaryIDs.begin(), primaryIDs.end(), hit->GetTID()) == primaryIDs.end()) {
  //         countPrimaryParticle++;
  //   primaryIDs.push_back(hit->GetTID());
  //         primaries.push_back(FPFParticle(hit->GetParticle(),
  //               hit->GetPID(), hit->GetTID(), countPrimaryParticle-1, 3, hit->GetParticleMass(),
  //               hit->GetTrackVertex().x(), hit->GetTrackVertex().y(), hit->GetTrackVertex().z(), 0,
  //               hit->GetInitMomentum().x(), hit->GetInitMomentum().y(), hit->GetInitMomentum().z(),
  //               GetTotalEnergy(hit->GetInitMomentum().x(), hit->GetInitMomentum().y(),
  //                 hit->GetInitMomentum().z(), hit->GetParticleMass())));
  //       }
  //     }
  //     // in case of tau decay pizero
  //     // decay products of this pizero are also counted as primary particles, mostly 2 gammas
  //     if (hit->GetIsTrackFromFSLPizero()) {
  //       tracksFromFSLDecayPizeroSecondary.insert(hit->GetTID());
  //       if (std::find(primaryIDs.begin(), primaryIDs.end(), hit->GetTID()) == primaryIDs.end()) {
  //         countPrimaryParticle++;
  //   primaryIDs.push_back(hit->GetTID());
  //         primaries.push_back(FPFParticle(hit->GetParticle(),
  //               hit->GetPID(), hit->GetTID(), countPrimaryParticle-1, 4, hit->GetParticleMass(),
  //               hit->GetTrackVertex().x(), hit->GetTrackVertex().y(), hit->GetTrackVertex().z(), 0,
  //               hit->GetInitMomentum().x(), hit->GetInitMomentum().y(), hit->GetInitMomentum().z(),
  //               GetTotalEnergy(hit->GetInitMomentum().x(), hit->GetInitMomentum().y(),
  //                 hit->GetInitMomentum().z(), hit->GetParticleMass())));
  //       }
  //     }
  //   } // end of hit loop
  // }

float_t AnalysisManager::GetTotalEnergy(float_t px, float_t py, float_t pz, float_t m)
{
  return TMath::Sqrt(px * px + py * py + pz * pz + m * m);
}
/*
void AnalysisManager::FillPseudoRecoVar()
{
  //  AngleToBeamDir, dEdx, dEdxInLAr ProngType
  std::cout << std::fixed << std::setw(10) << "PDG"
            << std::setw(12) << "Angle"
            << std::setw(13) << "TrackLength"
            << std::setw(13) << "ShowerLength"
            << std::setw(18) << "ShowerWidthInLAr"
            << std::setw(12) << "EInLAr"
            << std::setw(12) << "EInHadCal"
            << std::setw(12) << "dEdxInLAr"
            << std::setw(10) << "ProngType"
            << std::setw(12) << "Pz" << std::endl;

  for (int iPrim = 0; iPrim < nPrimaryParticle; ++iPrim)
  {
    if (ProngEInDetector[iPrim] > 0)
    {
      ShowerWidth[iPrim] = ShowerWidth[iPrim] / ProngEInDetector[iPrim];
    }
    if (ProngEInLAr[iPrim] > 0)
    {
      ShowerWidthInLAr[iPrim] = ShowerWidthInLAr[iPrim] / ProngEInLAr[iPrim];
    }

    double ShowerP = primaries[iPrim].P();
    double costheta = primaries[iPrim].Pz() / ShowerP;
    ProngAngleToBeamDir[iPrim] = TMath::ACos(costheta);

    ProngAvgdEdx[iPrim] = (ProngEInLAr[iPrim] +
                           ProngEInHadCal[iPrim] +
                           ProngEInMuonFinder[iPrim]) /
                          ShowerLength[iPrim];
    ProngAvgdEdxInLAr[iPrim] = ProngEInLAr[iPrim] / ShowerLengthInLAr[iPrim];

    std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(3);
    std::cout << std::setw(10) << primaries[iPrim].PDG();
    std::cout << std::setw(12) << ProngAngleToBeamDir[iPrim];
    std::cout << std::setw(13) << primaryTrackLength[iPrim];
    std::cout << std::setw(13) << ShowerLength[iPrim];
    std::cout << std::setw(18) << ShowerWidthInLAr[iPrim];
    std::cout << std::setw(12) << ProngEInLAr[iPrim];
    std::cout << std::setw(12) << ProngEInHadCal[iPrim];
    std::cout << std::setw(12) << ProngAvgdEdxInLAr[iPrim];
    std::cout << std::setw(10) << primaries[iPrim].ProngType();
    std::cout << std::setw(12) << primaries[iPrim].Pz() << std::endl;
  }

  slid::ShowerLID *shwlid = new slid::ShowerLID(pm3D->Get3DPixelMap(),
                                                neutrino.NuVx(), neutrino.NuVy(), neutrino.NuVz(), 0., 0., 1.);
  Double_t *ptr_dedx = shwlid->GetTotalDedxLongitudinal();
  std::copy(ptr_dedx, ptr_dedx + 3000, TotalDedxLongitudinal);

  for (int iPrim = 0; iPrim < nPrimaryParticle; ++iPrim)
  {
    directionfitter::LinearFit *linFit = new directionfitter::LinearFit(
        pm3D->Get2DPixelMapZX(iPrim + 1),
        pm3D->Get2DPixelMapZY(iPrim + 1),
        pm3D->Get2DVtxPixelMapZX(iPrim + 1),
        pm3D->Get2DVtxPixelMapZY(iPrim + 1),
        neutrino.NuVx(), neutrino.NuVy(), neutrino.NuVz(),
        primaries[iPrim].Vx(), primaries[iPrim].Vy(), primaries[iPrim].Vz());
    dir_pol_x[iPrim] = linFit->GetDir().X();
    dir_pol_y[iPrim] = linFit->GetDir().Y();
    dir_pol_z[iPrim] = linFit->GetDir().Z();
    dir_coc_x[iPrim] = linFit->GetCOCDir().X();
    dir_coc_y[iPrim] = linFit->GetCOCDir().Y();
    dir_coc_z[iPrim] = linFit->GetCOCDir().Z();
    delete linFit;
  }

  delete shwlid;
}*/
