#include <vector>
#include <functional>
#include <iostream>
#include <string>
#include <unordered_map>
#include <iomanip>
#include <random>

#include "AnalysisManager.hh"
#include "geometry/GeometricalParameters.hh"
#include "LArBoxSD.hh"
#include "LArBoxHit.hh"
#include "PrimaryParticleInformation.hh"
#include "reco/PCAAnalysis3D.hh"
#include "reco/Cluster3D.hh"
#include "reco/LinearFit.hh"
#include "reco/ShowerLID.hh"
#include "reco/Barcode.hpp"
#include "reco/CircleFit.hh"
#include "FPFParticle.hh"
#include "FPFNeutrino.hh"

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

using namespace hep_hpc::hdf5;

AnalysisManager* AnalysisManager::instance = 0;

AnalysisManager* AnalysisManager::GetInstance() {
  if (!instance) {
    std::cout<<"AnalysisManager: Re-initialization"<<std::endl;
    instance = new AnalysisManager();
  }
  return instance;
}

AnalysisManager::AnalysisManager() {
  evt = 0;
  trk = 0;
  messenger = new AnalysisManagerMessenger(this);
  m_filename = "test.root";
  m_addDiffusion = "false";
  m_saveHit = false;
  m_saveFlare = false;
  m_saveTrack = false;
  m_save3DEvd = false;
  m_save2DEvd = false;
  m_circularFit = false;

  fH5Filename = "test.h5";
}

AnalysisManager::~AnalysisManager() {;}

void AnalysisManager::bookEvtTree() {
  evt = new TTree("evt", "evtTreeInfo");
  evt->Branch("evtID"                     , &evtID                     , "evtID/I");
  evt->Branch("FPFNeutrino"                , &neutrino, 96000, 0);
  evt->Branch("FPFParticle"                , &primaries, 96000, 0);
  evt->Branch("TotalDedxLongitudinal"     , TotalDedxLongitudinal      , "TotalDedxLongitudinal[3000]/D");
  evt->Branch("TrueTotalDedxLongitudinal" , TrueTotalDedxLongitudinal  , "TrueTotalDedxLongitudinal[3000]/D");
  evt->Branch("nPrimaryParticle"          , &nPrimaryParticle          , "nPrimaryParticle/I");
  evt->Branch("primaryParentPDG"          , primaryParentPDG           , "primaryParentPDG[nPrimaryParticle]/I");
  evt->Branch("primaryTrackLength"        , primaryTrackLength         , "primaryTrackLength[nPrimaryParticle]/D");
  evt->Branch("primaryTrackLengthInTPC"   , primaryTrackLengthInTPC    , "primaryTrackLengthInTPC[nPrimaryParticle]/D");
  evt->Branch("ProngEInDetector"          , ProngEInDetector           , "ProngEInDetector[nPrimaryParticle]/D");
  evt->Branch("ProngEInLAr"               , ProngEInLAr                , "ProngEInLAr[nPrimaryParticle]/D");
  evt->Branch("ProngEInHadCal"            , ProngEInHadCal             , "ProngEInHadCal[nPrimaryParticle]/D");
  evt->Branch("ProngEInMuonFinder"        , ProngEInMuonFinder         , "ProngEInMuonFinder[nPrimaryParticle]/D");
  evt->Branch("ProngEInMuonFinderLayer1X" , ProngEInMuonFinderLayer1X  , "ProngEInMuonFinderLayer1X[nPrimaryParticle]/D");
  evt->Branch("ProngEInMuonFinderLayer1Y" , ProngEInMuonFinderLayer1Y  , "ProngEInMuonFinderLayer1Y[nPrimaryParticle]/D");
  evt->Branch("ProngEInMuonFinderLayer2X" , ProngEInMuonFinderLayer2X  , "ProngEInMuonFinderLayer2X[nPrimaryParticle]/D");
  evt->Branch("ProngEInMuonFinderLayer2Y" , ProngEInMuonFinderLayer2Y  , "ProngEInMuonFinderLayer2Y[nPrimaryParticle]/D");
  evt->Branch("ProngAngleToBeamDir"       , ProngAngleToBeamDir        , "ProngAngleToBeamDir[nPrimaryParticle]/D");
  evt->Branch("ShowerLength"              , ShowerLength               , "ShowerLength[nPrimaryParticle]/D");
  evt->Branch("ShowerLengthInLAr"         , ShowerLengthInLAr          , "ShowerLengthInLAr[nPrimaryParticle]/D");
  evt->Branch("ShowerWidth"               , ShowerWidth                , "ShowerWidth[nPrimaryParticle]/D");
  evt->Branch("ShowerWidthInLAr"          , ShowerWidthInLAr           , "ShowerWidthInLAr[nPrimaryParticle]/D");
  evt->Branch("ProngAvgdEdx"              , ProngAvgdEdx               , "ProngAvgdEdx[nPrimaryParticle]/D");
  evt->Branch("ProngAvgdEdxInLAr"         , ProngAvgdEdxInLAr          , "ProngAvgdEdxInLAr[nPrimaryParticle]/D");
  evt->Branch("ProngdEdxAlongTrack"       , ProngdEdxAlongTrack        , "ProngdEdxAlongTrack[nPrimaryParticle][100]/D");
  evt->Branch("ProngdEdxTrackLength"      , ProngdEdxTrackLength       , "ProngdEdxTrackLength[nPrimaryParticle][100]/I");
  evt->Branch("dir_pol_x"                 , dir_pol_x                  , "dir_pol_x[nPrimaryParticle]/D");
  evt->Branch("dir_pol_y"                 , dir_pol_y                  , "dir_pol_y[nPrimaryParticle]/D");
  evt->Branch("dir_pol_z"                 , dir_pol_z                  , "dir_pol_z[nPrimaryParticle]/D");
  evt->Branch("dir_coc_x"                 , dir_coc_x                  , "dir_coc_x[nPrimaryParticle]/D");
  evt->Branch("dir_coc_y"                 , dir_coc_y                  , "dir_coc_y[nPrimaryParticle]/D");
  evt->Branch("dir_coc_z"                 , dir_coc_z                  , "dir_coc_z[nPrimaryParticle]/D");

  evt->Branch("nHits"                     , &nHits                     , "nHits/I");
  evt->Branch("sparseFractionMem"         , &sparseFractionMem         , "sparseFractionMem/D");
  evt->Branch("sparseFractionBins"        , &sparseFractionBins        , "sparseFractionBins/D");
  if (m_saveHit) {
    evt->Branch("HitTID"                  , HitTID                     , "HitTID[nHits]/I");
    evt->Branch("HitPID"                  , HitPID                     , "HitPID[nHits]/I");
    evt->Branch("HitPDG"                  , HitPDG                     , "HitPDG[nHits]/I");
    evt->Branch("HitTrackStatus"          , HitTrackStatus             , "HitTrackStatus[nHits]/I");
    evt->Branch("HitPrePositionX"         , HitPrePositionX            , "HitPrePositionX[nHits]/D");
    evt->Branch("HitPrePositionY"         , HitPrePositionY            , "HitPrePositionY[nHits]/D");
    evt->Branch("HitPrePositionZ"         , HitPrePositionZ            , "HitPrePositionZ[nHits]/D");
    evt->Branch("HitPosPositionX"         , HitPosPositionX            , "HitPosPositionX[nHits]/D");
    evt->Branch("HitPosPositionY"         , HitPosPositionY            , "HitPosPositionY[nHits]/D");
    evt->Branch("HitPosPositionZ"         , HitPosPositionZ            , "HitPosPositionZ[nHits]/D");
    evt->Branch("HitEdep"                 , HitEdep                    , "HitEdep[nHits]/D");
  }

  evt->Branch("edepInLAr"                 , &edepInLAr                 , "edepInLAr/D");
  evt->Branch("edepInHadCalX"             , &edepInHadCalX             , "edepInHadCalX/D");
  evt->Branch("edepInHadCalY"             , &edepInHadCalY             , "edepInHadCalY/D");
  evt->Branch("edepInMuonFinderX"         , &edepInMuonFinderX         , "edepInMuonFinderX/D");
  evt->Branch("edepInMuonFinderY"         , &edepInMuonFinderY         , "edepInMuonFinderY/D");
  evt->Branch("edepInHadAborb"            , &edepInHadAborb            , "edepInHadAborb/D");
  evt->Branch("edepInMuonFinderAbsorb"    , &edepInMuonFinderAbsorb    , "edepInMuonFinderAbsorb/D");
  evt->Branch("missCountedEnergy"         , &missCountedEnergy         , "missCountedEnergy/D");

  evt->Branch("nFromFSLParticles"         , &nFromFSLParticles         , "nFromFSLParticles/I");

  if(m_circularFit){
    evt->Branch("circNhits"               , &circNhits                 , "circNhits/I");
    evt->Branch("trkNhits"                , &trkNhits                  , "trkNhits/I");
    evt->Branch("circXc"                  , &xc                        , "circXc/D");
    evt->Branch("circZc"                  , &zc                        , "circZc/D");
    evt->Branch("circRc"                  , &rc                        , "circRc/D");
    evt->Branch("circp0"                  , &p0                        , "circp0/D");
    evt->Branch("circp1"                  , &p1                        , "circp1/D");
    evt->Branch("circcosDip"              , &cosDip                    , "circcosDip/D");
    evt->Branch("Nmagnets"                , &Nmagnets                  , "Nmagnets/I");
    evt->Branch("trkp0"                   , &trkp0                     , "trkp0/D");
    evt->Branch("trkp1"                   , &trkp1                     , "trkp1/D");
    evt->Branch("trkcosDip"               , &trkcosDip                 , "trkcosDip/D");
    evt->Branch("magnetZs"                , &magzpos);
    evt->Branch("trkXc"                   , &trkxc);
    evt->Branch("trkZc"                   , &trkzc);
    evt->Branch("trkRc"                   , &trkrc);
    evt->Branch("trkqIn"                  , &trkqIn);
    evt->Branch("trkmIn"                  , &trkmIn);
    evt->Branch("trkqOut"                 , &trkqOut);
    evt->Branch("trkmOut"                 , &trkmOut);
    evt->Branch("circHitXFSL"             , &hitXFSL);
    evt->Branch("circHitYFSL"             , &hitYFSL);
    evt->Branch("circHitZFSL"             , &hitZFSL);
    evt->Branch("circHitPFSL"             , &hitPFSL);
    evt->Branch("trkHitXFSL"              , &trkXFSL);
    evt->Branch("trkHitYFSL"              , &trkYFSL);
    evt->Branch("trkHitZFSL"              , &trkZFSL);
    evt->Branch("trkHitPFSL"              , &trkPFSL);
  }

}

void AnalysisManager::bookTrkTree() {

  trk = new TTree("trk", "trkTreeInfo");
  trk->Branch("evtID", &evtID, "evtID/I");
  trk->Branch("trackTID", &trackTID, "trackTID/I");
  trk->Branch("trackPID", &trackPID, "trackPID/I");
  trk->Branch("trackPDG", &trackPDG, "trackPDG/I");
  trk->Branch("trackKinE", &trackKinE, "trackKinE/D");
  trk->Branch("trackNPoints", &trackNPoints, "trackNPoints/I");
  trk->Branch("trackPointX", &trackPointX);
  trk->Branch("trackPointY", &trackPointY);
  trk->Branch("trackPointZ", &trackPointZ);

}


void AnalysisManager::bookPrimTree() {


  prim = new TTree("primaries", "primPTreeInfo");
  prim->Branch("primEvtID"					, &primEvtID 			, "primEvtID/I");
  prim->Branch("primVtxID" 					, &primVtxID 			, "primVtxID/I");
  prim->Branch("primPDG" 					, &primPDG				, "primPDG/I");
  prim->Branch("primParticleID" 			, &primParticleID		, "primParticleID/I");
  prim->Branch("primTrackID" 				, &primTrackID 			, "primTrackID/I");

  prim->Branch("primM" 						, &primM				, "primM/F");
  prim->Branch("primQ" 						, &primQ				, "primQ/F");
  prim->Branch("primEta" 					, &primEta				, "primEta/F");
  prim->Branch("primPhi" 					, &primPhi				, "primPhi/F");
  prim->Branch("primPt" 					, &primPt				, "primPt/F");
  prim->Branch("primP" 						, &primP				, "primP/F");

  prim->Branch("primVx" 					, &primVx				, "primVx/F");//position
  prim->Branch("primVy" 					, &primVy				, "primVy/F");
  prim->Branch("primVz" 					, &primVz				, "primVz/F");
  prim->Branch("primVt" 					, &primVt				, "primVt/F");

  prim->Branch("primPx" 					, &primPx				, "primPx/F");//momentum
  prim->Branch("primPy" 					, &primPy				, "primPy/F");
  prim->Branch("primPz" 					, &primPz				, "primPz/F");

  prim->Branch("primE"						, &primE				, "primE/F"); //initial total energy
  prim->Branch("primKE"						, &primKE				, "primKE/F"); //initial kinetic energy

}

void AnalysisManager::bookFLArEHitTree() {
  flarHit = new TTree("FLArE_hits", "flarHitTreeInfo");

  flarHit->Branch("flareTrackID"			, &flareTrackID 		, "flareTrackID/I");
  flarHit->Branch("flareParticleID" 		, &flareParticleID		, "flareParticleID/I");
  flarHit->Branch("flareParentID" 			, &flareParentID		, "flareParentID/I");
  flarHit->Branch("flarePDG" 				, &flarePDG				, "flarePDG/I");
  flarHit->Branch("flareCopyNum" 			, &flareCopyNum			, "flareCopyNum/I");

  flarHit->Branch("flareT" 					, &flareT				, "flareT/I");

  flarHit->Branch("flareX" 					, &flareX				, "flareX/D");//Pre-position
  flarHit->Branch("flareY" 					, &flareY				, "flareY/D");
  flarHit->Branch("flareZ" 					, &flareZ				, "flareZ/D");

  flarHit->Branch("flarePx" 				, &flarePx				, "flarePx/D");//Post-position
  flarHit->Branch("flarePy" 				, &flarePy				, "flarePy/D");
  flarHit->Branch("flarePz" 				, &flarePz				, "flarePz/D");

  flarHit->Branch("flareDeltaPx" 			, &flareDeltaPx			, "flareDeltaPx/D");
  flarHit->Branch("flareDeltaPy" 			, &flareDeltaPy			, "flareDeltaPy/D");
  flarHit->Branch("flareDeltaPz" 			, &flareDeltaPz			, "flareDeltaPz/D");

  flarHit->Branch("flareEdep" 				, &flareEdep			, "flareEdep/D");

}


void AnalysisManager::bookFLArEHitTree() {
  flarHit = new TTree("FLArE_hits", "flarHitTreeInfo");

  flarHit->Branch("flareTrackID"			, &flareTrackID 		, "flareTrackID/I");
  flarHit->Branch("flareParticleID" 		, &flareParticleID		, "flareParticleID/I");
  flarHit->Branch("flareParentID" 			, &flareParentID		, "flareParentID/I");
  flarHit->Branch("flarePDG" 				, &flarePDG				, "flarePDG/I");
  flarHit->Branch("flareCopyNum" 			, &flareCopyNum			, "flareCopyNum/I");

  flarHit->Branch("flareT" 					, &flareT				, "flareT/I");

  flarHit->Branch("flareX" 					, &flareX				, "flareX/D");//Pre-position
  flarHit->Branch("flareY" 					, &flareY				, "flareY/D");
  flarHit->Branch("flareZ" 					, &flareZ				, "flareZ/D");

  flarHit->Branch("flarePx" 				, &flarePx				, "flarePx/D");//Post-position
  flarHit->Branch("flarePy" 				, &flarePy				, "flarePy/D");
  flarHit->Branch("flarePz" 				, &flarePz				, "flarePz/D");

  flarHit->Branch("flareDeltaPx" 			, &flareDeltaPx			, "flareDeltaPx/D");
  flarHit->Branch("flareDeltaPy" 			, &flareDeltaPy			, "flareDeltaPy/D");
  flarHit->Branch("flareDeltaPz" 			, &flareDeltaPz			, "flareDeltaPz/D");

  flarHit->Branch("flareEdep" 				, &flareEdep			, "flareEdep/D");

}






void AnalysisManager::BeginOfRun() {

  G4cout<<"TTree is booked and the run has been started"<<G4endl;
  if (thefile) {
    delete thefile;
  }
  thefile = new TFile(m_filename.c_str(), "RECREATE");
  bookEvtTree();
  bookPrimTree();
  if(m_saveTrack) bookTrkTree();

  fH5Filename = m_filename;
  if(fH5Filename.find(".root") != std::string::npos) {
    const size_t pos = fH5Filename.find(".root");
    fH5Filename.resize(pos);
  }
  fH5Filename += ".h5";

  fH5file = hep_hpc::hdf5::File(fH5Filename, H5F_ACC_TRUNC);

  SDNamelist = GeometricalParameters::Get()->GetSDNamelist();
  NumberOfSDs = SDNamelist.size();
  G4cout<<"Number of SDs : "<<NumberOfSDs<<G4endl;
  for (auto sdname : SDNamelist){
    G4cout<<sdname.first<<" "<<sdname.second<<G4endl;
	//make flareHit tree if it's in this list of names:
	if(sdname.second == "lArBoxSD/lar_box"){
		m_saveFlare = true;
		flareDir = thefile->mkdir("flare");
		bookFLArEHitTree();
	}
  }
}

void AnalysisManager::EndOfRun() {
  thefile->cd();
  evt->Write();
  prim->Write();
  if (m_saveFlare){
	thefile->cd(flareDir->GetName());
  	flarHit->Write();
	thefile->cd();
  }
  if(m_saveTrack) trk->Write();
  thefile->Close();
  fH5file.close();
}

void AnalysisManager::BeginOfEvent() {
  neutrino = FPFNeutrino();
  primaries.clear();
  primaryIDs.clear();
  nHits                        = 0;
  sparseFractionMem            = -1;
  sparseFractionBins           = -1;
  edepInLAr                    = 0;
  edepInHadCalX                = 0;
  edepInHadCalY                = 0;
  edepInMuonFinderX            = 0;
  edepInMuonFinderY            = 0;
  edepInHadAborb               = 0;
  edepInMuonFinderAbsorb       = 0;
  missCountedEnergy            = 0;
  nPrimaryParticle             = 0;
  nTestNPrimaryTrack           = 0;
  nFromFSLParticles            = 0;
  nFromFSPizeroParticles       = 0;
  nFromFSLDecayPizeroParticles = 0;
  countPrimaryParticle         = 0;
  for (G4int j= 0; j< 3000; ++j) {
    TotalDedxLongitudinal[j] = 0;
    TrueTotalDedxLongitudinal[j] = 0;
  }
  for (G4int i= 0; i< 1000; ++i) {
    primaryParentPDG[i]          = 0;
    primaryTrackLength[i]        = 0;
    primaryTrackLengthInTPC[i]   = 0;
    ProngEInDetector[i]          = 0;
    ProngEInLAr[i]               = 0;
    ProngEInHadCal[i]            = 0;
    ProngEInMuonFinder[i]        = 0;
    ProngEInMuonFinderLayer1X[i] = 0;
    ProngEInMuonFinderLayer1Y[i] = 0;
    ProngEInMuonFinderLayer2X[i] = 0;
    ProngEInMuonFinderLayer2Y[i] = 0;
    ProngAngleToBeamDir[i]       = -1;
    ShowerLength[i]              = -1;
    ShowerLengthInLAr[i]         = -1;
    ShowerWidth[i]               = 0;
    ShowerWidthInLAr[i]          = 0;
    ProngAvgdEdx[i]              = -1;
    ProngAvgdEdxInLAr[i]         = -1;
    for (G4int j= 0; j< 100; ++j) {
      ProngdEdxAlongTrack[i][j] = 0;
      ProngdEdxTrackLength[i][j] = -1;
    }
    dir_pol_x[i]                 = -999;
    dir_pol_y[i]                 = -999;
    dir_pol_z[i]                 = -999;
    dir_coc_x[i]                 = -999;
    dir_coc_y[i]                 = -999;
    dir_coc_z[i]                 = -999;
    fromFSLParticlePDG[i]        = 0;
  }
  if (m_saveHit) {
    for (G4int i= 0; i< 40000000; ++i) {
      HitTID[i] = -1;
      HitPID[i] = -1;
      HitPDG[i] = 0;
      HitTrackStatus[i]  = -1;
      HitPrePositionX[i] = -999;
      HitPrePositionY[i] = -999;
      HitPrePositionZ[i] = -999;
      HitPosPositionX[i] = -999;
      HitPosPositionY[i] = -999;
      HitPosPositionZ[i] = -999;
      HitEdep[i]         = 0;
    }
  }

  // vectors that need to be cleared for a new event
  allTracksPTPair.clear();
  trackClusters.clear();
  tracksFromFSL.clear();
  tracksFromFSLSecondary.clear();
  tracksFromFSPizeroSecondary.clear();
  tracksFromFSLDecayPizeroSecondary.clear();
  fPrimIdxFSL = -1;

  magzpos.clear();
  trkxc.clear();
  trkzc.clear();
  trkrc.clear();
  trkqIn.clear();
  trkmIn.clear();
  trkqOut.clear();
  trkmOut.clear();

  hitXFSL.clear();
  hitYFSL.clear();
  hitZFSL.clear();
  hitPFSL.clear();
  trkXFSL.clear();
  trkYFSL.clear();
  trkZFSL.clear();
  trkPFSL.clear();

  trackTID = -1;
  trackPID = -1;
  trackPDG = -1;
  trackKinE = -1;
  trackNPoints = -1;
  trackPointX.clear();
  trackPointY.clear();
  trackPointZ.clear();



  primEvtID = 0;
  primVtxID = 0;
  primParticleID = 0;
  primTrackID = 0;
  primPDG = 0;

  primM = 0;
  primQ = 0;
  primEta = 0;
  primPhi = 0;
  primPt = 0;
  primP = 0;

  primVx = 0;
  primVy = 0;
  primVz = 0;
  primVt = 0;

  primPx = 0;
  primPy = 0;
  primPz = 0;

  primE = 0;
  primKE = 0;
  //add parentID later?

  if (m_saveFlare) {
  	flareTrackID = 0;
  	flareParentID = 0;
  	flarePDG = 0;
  	flareCopyNum = 0;
	flareParticleID = 0;

	flareT = 0;

  	flareX = 0;
  	flareY = 0;
  	flareZ = 0;

  	flarePx = 0;
  	flarePy = 0;
  	flarePz = 0;

  	flareDeltaPx = 0;
  	flareDeltaPy = 0;
  	flareDeltaPz = 0;

  	flareEdep = 0;
  }

}

void AnalysisManager::EndOfEvent(const G4Event* event) {
  /// evtID
  evtID = event->GetEventID();
  primEvtID = evtID;

  /// loop over the vertices, and then over primary particles,
  /// neutrino truth info from event generator.
  for (G4int ivtx = 0; ivtx < event->GetNumberOfPrimaryVertex(); ++ivtx) {

    primVtxID = ivtx;
    for (G4int ipp = 0; ipp < event->GetPrimaryVertex(ivtx)->GetNumberOfParticle(); ++ipp) {

      G4PrimaryParticle* primary_particle =
        event->GetPrimaryVertex(ivtx)->GetPrimary(ipp);
      if (primary_particle) {
        PrimaryParticleInformation* primary_particle_info =
          dynamic_cast<PrimaryParticleInformation*>(primary_particle->GetUserInformation());
        primary_particle_info->Print();

		primTrackID = primary_particle_info->GetPartID(); //confirm matches track's id later?

		auto particleId = ActsFatras::Barcode();
		particleId.setVertexPrimary(ivtx);
		particleId.setGeneration(0);
		particleId.setSubParticle(0);
		particleId.setParticle(primTrackID-1);
		primParticleID = particleId.value();


		primPDG = primary_particle_info->GetPDG();

		primVx = primary_particle_info->GetVertexMC().x();
		primVy = primary_particle_info->GetVertexMC().y();
		primVz = primary_particle_info->GetVertexMC().z();
		primVt = 0;

		float_t p_x = primary_particle_info->GetMomentumMC().x()/1000;
		float_t p_y = primary_particle_info->GetMomentumMC().y()/1000;
		float_t p_z = primary_particle_info->GetMomentumMC().z()/1000;

		primPx = p_x; //div by 1000 for MeV->GeV
		primPy = p_y;
		primPz = p_z;
		primM = primary_particle_info->GetMass()/1000;
		primQ = primary_particle_info->GetCharge();

		TLorentzVector p4;
		G4double energy = GetTotalEnergy(p_x,p_y,p_z,primary_particle_info->GetMass());
		p4.SetPx(p_x);
		p4.SetPy(p_y);
		p4.SetPz(p_z);
		p4.SetE(energy);

		primEta = p4.Eta();
		primPhi = p4.Phi();
		primPt = p4.Pt();
		primP = p4.P();

		primE = energy;
		primKE = energy-primM;

	  	prim->Fill();


      }
    }
  }
  nPrimaryVertex   = event->GetNumberOfPrimaryVertex();
  std::cout<<"\nnumber of primary vertices  : "<<nPrimaryVertex<<std::endl;
  std::cout<<"=============>"<<neutrino.NuPDG()<<" "<<neutrino.NuFSLPDG()<<std::endl;

  // Get the hit collections
  // If there is no hit collection, there is nothing to be done
  hcofEvent = event->GetHCofThisEvent();
  if (!hcofEvent) return;

  // loop over all sensitive detectors to get all the hit collections
  // fill the truth tree for primary particles
  for (auto sdname : SDNamelist) {
    FillPrimaryTruthTree(sdname.first, sdname.second);
  }
}

void AnalysisManager::FillPrimaryTruthTree(G4int sdId, std::string sdName)
{

  // Get and cast hit collection with LArBoxHits
  LArBoxHitsCollection* hitCollection = dynamic_cast<LArBoxHitsCollection*>(hcofEvent->GetHC(sdId));
  if (hitCollection ) {
    for (auto hit: *hitCollection->GetVector()) {
      nHits++;

      double pre_x  = hit->GetPreStepPosition().x(); //using preX-Z as recorded positions
      double pre_y  = hit->GetPreStepPosition().y();
      double pre_z  = hit->GetPreStepPosition().z();

	  UInt_t tid = hit->GetTID();
	  UInt_t pid = hit->GetPID();
	  UInt_t PDG = hit->GetPDG();
	  UInt_t copyNum = hit->GetCopyNum();

	  UInt_t time = hit->GetTime();

	  double Px = hit->GetInitMomentum().x();
  	  double Py = hit->GetInitMomentum().y();
	  double Pz = hit->GetInitMomentum().z();

	  double dPx = hit->GetInitMomentum().x();
  	  double dPy = hit->GetInitMomentum().y();
	  double dPz = hit->GetInitMomentum().z();


	  double edep =  hit->GetEdep();

	  auto particleId = ActsFatras::Barcode();
	  particleId.setVertexPrimary(1);//fix this value
	  particleId.setGeneration(parentID);
	  particleId.setSubParticle(0);
	  particleId.setParticle(trackID);
	  partID = particleId.value();


      // energy deposition in different volumes of the detector
      if (sdName == "lArBoxSD/lar_box"){
	  	flareTrackID = tid
	  	flareParentID = pid;
	  	flarePDG = PDG;
	  	flareCopyNum = copyNum;
		flareParticleID = partID;
		flareT = time;
	  	flareX = pre_x;
	  	flareY = pre_y;
	  	flareZ = pre_z;
	  	flarePx = Px;
	  	flarePy = Py;
	  	flarePz = Pz;
	  	flareDeltaPx = dPx;
	  	flareDeltaPy = dPy;
	  	flareDeltaPz = dPz;
	  	flareEdep = edep;//should this be += or =?


	  	flarHit->Fill();
	  }

      if (sdName == "HadCalXSD/lar_box"){
	  	hadXTrackID = tid
	  	hadXParentID = pid;
	  	hadXPDG = PDG;
	  	hadXCopyNum = copyNum;
		hadXParticleID = partID;
		hadXT = time;
	  	hadXx = pre_x;
	  	hadXy = pre_y;
	  	hadXz = pre_z;
	  	hadXPx = Px;
	  	hadXPy = Py;
	  	hadXPz = Pz;
	  	hadXDeltaPx = dPx;
	  	hadXDeltaPy = dPy;
	  	hadXDeltaPz = dPz;
	  	hadXEdep = edep;

	  	hadXHit->Fill();
	  }


    } // end of hit loop
  }
}

void AnalysisManager::FillTrueEdep(G4int sdId, std::string sdName)
{
}

double AnalysisManager::GetTotalEnergy(double px, double py, double pz, double m) {
  return TMath::Sqrt(px*px+py*py+pz*pz+m*m);
}

void AnalysisManager::FillPseudoRecoVar() {
  //  AngleToBeamDir, dEdx, dEdxInLAr ProngType
  std::cout<<std::fixed<<std::setw(10)<<"PDG"
    <<std::setw(12)<<"Angle"
    <<std::setw(13)<<"TrackLength"
    <<std::setw(13)<<"ShowerLength"
    <<std::setw(18)<<"ShowerWidthInLAr"
    <<std::setw(12)<<"EInLAr"
    <<std::setw(12)<<"EInHadCal"
    <<std::setw(12)<<"dEdxInLAr"
    <<std::setw(10)<<"ProngType"
    <<std::setw(12)<<"Pz"<<std::endl;

  for (int iPrim= 0; iPrim< nPrimaryParticle; ++iPrim) {
    if (ProngEInDetector[iPrim]>0) {
      ShowerWidth[iPrim] = ShowerWidth[iPrim]/ProngEInDetector[iPrim];
    }
    if (ProngEInLAr[iPrim]>0) {
      ShowerWidthInLAr[iPrim] = ShowerWidthInLAr[iPrim]/ProngEInLAr[iPrim];
    }

    double ShowerP = primaries[iPrim].P();
    double costheta = primaries[iPrim].Pz()/ShowerP;
    ProngAngleToBeamDir[iPrim] = TMath::ACos(costheta);

    ProngAvgdEdx[iPrim] = (ProngEInLAr[iPrim] +
                           ProngEInHadCal[iPrim] +
                           ProngEInMuonFinder[iPrim])/ShowerLength[iPrim];
    ProngAvgdEdxInLAr[iPrim] = ProngEInLAr[iPrim]/ShowerLengthInLAr[iPrim];

    std::cout<<std::setiosflags(std::ios::fixed)<<std::setprecision(3);
    std::cout<<std::setw(10)<<primaries[iPrim].PDG();
    std::cout<<std::setw(12)<<ProngAngleToBeamDir[iPrim];
    std::cout<<std::setw(13)<<primaryTrackLength[iPrim];
    std::cout<<std::setw(13)<<ShowerLength[iPrim];
    std::cout<<std::setw(18)<<ShowerWidthInLAr[iPrim];
    std::cout<<std::setw(12)<<ProngEInLAr[iPrim] ;
    std::cout<<std::setw(12)<<ProngEInHadCal[iPrim];
    std::cout<<std::setw(12)<<ProngAvgdEdxInLAr[iPrim];
    std::cout<<std::setw(10)<<primaries[iPrim].ProngType();
    std::cout<<std::setw(12)<<primaries[iPrim].Pz()<<std::endl;
  }

  slid::ShowerLID* shwlid = new slid::ShowerLID(pm3D->Get3DPixelMap(),
      neutrino.NuVx(), neutrino.NuVy(), neutrino.NuVz(), 0., 0., 1.);
  Double_t* ptr_dedx = shwlid->GetTotalDedxLongitudinal();
  std::copy(ptr_dedx, ptr_dedx+3000, TotalDedxLongitudinal);

  for(int iPrim= 0; iPrim< nPrimaryParticle; ++iPrim) {
    directionfitter::LinearFit* linFit = new directionfitter::LinearFit(
        pm3D->Get2DPixelMapZX(iPrim+1),
        pm3D->Get2DPixelMapZY(iPrim+1),
        pm3D->Get2DVtxPixelMapZX(iPrim+1),
        pm3D->Get2DVtxPixelMapZY(iPrim+1),
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
}
