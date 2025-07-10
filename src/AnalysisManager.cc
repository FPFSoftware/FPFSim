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
#include "FASER2TrackerSD.hh"
#include "LArBoxHit.hh"
#include "PrimaryParticleInformation.hh"
#include "reco/PCAAnalysis3D.hh"
#include "reco/Cluster3D.hh"
#include "reco/LinearFit.hh"
#include "reco/ShowerLID.hh"
#include "reco/Barcode.hh"
#include "reco/CircleFit.hh"
#include "reco/Barcode.hh"
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
  m_saveActs = true;
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

  // Defaults to true if FASER2 is enabled. But give the option to disable output if required
  m_saveActs = m_saveActs && GeometricalParameters::Get()->GetAddFASER2();
  if (m_saveActs) {
    //* Acts Hits Tree [i == unsigned int; F == float; l == Long unsigned 64 int]
    acts_hits_tree = new TTree("hits", "ActsHitsTree");
    acts_hits_tree->Branch("event_id"               , &ActsHitsEventID,     "event_id/i");
    acts_hits_tree->Branch("geometry_id"            , &ActsHitsGeometryID,  "geometryid/l");
    acts_hits_tree->Branch("particle_id"            , &ActsHitsParticleID,  "particle_id/l");
    acts_hits_tree->Branch("tx"                     , &ActsHitsX,           "tx/F");
    acts_hits_tree->Branch("ty"                     , &ActsHitsY,           "ty/F");
    acts_hits_tree->Branch("tz"                     , &ActsHitsZ,           "tz/F");
    acts_hits_tree->Branch("tt"                     , &ActsHitsT,           "tt/F");
    acts_hits_tree->Branch("tpx"                    , &ActsHitsPx,          "tpx/F");
    acts_hits_tree->Branch("tpy"                    , &ActsHitsPy,          "tpy/F");
    acts_hits_tree->Branch("tpz"                    , &ActsHitsPz,          "tpz/F");
    acts_hits_tree->Branch("te"                     , &ActsHitsE,           "tpe/F");
    acts_hits_tree->Branch("deltapx"                , &ActsHitsDeltaPx,     "deltapx/F");
    acts_hits_tree->Branch("deltapy"                , &ActsHitsDeltaPy,     "deltapy/F");
    acts_hits_tree->Branch("deltapz"                , &ActsHitsDeltaPz,     "deltapz/F");
    acts_hits_tree->Branch("deltae"                 , &ActsHitsDeltaE,      "deltae/F");
    acts_hits_tree->Branch("index"                  , &ActsHitsIndex,       "index/I");
    acts_hits_tree->Branch("volume_id"              , &ActsHitsVolumeID,    "volume_id/i");
    acts_hits_tree->Branch("boundary_id"            , &ActsHitsBoundaryID,  "boundary_id/i");
    acts_hits_tree->Branch("layer_id"               , &ActsHitsLayerID,     "layer_id/i");
    acts_hits_tree->Branch("approach_id"            , &ActsHitsApproachID,  "approach_id/i");
    acts_hits_tree->Branch("sensitive_id"           , &ActsHitsSensitiveID, "sensitive_id/i");

    //* Acts truth particle tree
    acts_particles_tree = new TTree("particles", "ActsParticlesTree");
    acts_particles_tree->Branch("event_id"               , &ActsHitsEventID,     "event_id/i");
    acts_particles_tree->Branch("particle_id"            , &ActsParticlesParticleId);
    acts_particles_tree->Branch("particle_type"          , &ActsParticlesParticleType);
    acts_particles_tree->Branch("process"                , &ActsParticlesProcess);
    acts_particles_tree->Branch("vx"                     , &ActsParticlesVx);
    acts_particles_tree->Branch("vy"                     , &ActsParticlesVy);
    acts_particles_tree->Branch("vz"                     , &ActsParticlesVz);
    acts_particles_tree->Branch("vt"                     , &ActsParticlesVt);
    acts_particles_tree->Branch("px"                     , &ActsParticlesPx);
    acts_particles_tree->Branch("py"                     , &ActsParticlesPy);
    acts_particles_tree->Branch("pz"                     , &ActsParticlesPz);
    acts_particles_tree->Branch("m"                      , &ActsParticlesM);
    acts_particles_tree->Branch("q"                      , &ActsParticlesQ);
    acts_particles_tree->Branch("eta"                    , &ActsParticlesEta);
    acts_particles_tree->Branch("phi"                    , &ActsParticlesPhi);
    acts_particles_tree->Branch("pt"                     , &ActsParticlesPt);
    acts_particles_tree->Branch("p"                      , &ActsParticlesP);
    acts_particles_tree->Branch("vertex_primary"         , &ActsParticlesVertexPrimary);
    acts_particles_tree->Branch("vertex_secondary"       , &ActsParticlesVertexSecondary);
    acts_particles_tree->Branch("particle"               , &ActsParticlesParticle);
    acts_particles_tree->Branch("generation"             , &ActsParticlesGeneration);
    acts_particles_tree->Branch("sub_particle"           , &ActsParticlesSubParticle);
    acts_particles_tree->Branch("e_loss"                 , &ActsParticlesELoss);
    acts_particles_tree->Branch("total_x0"               , &ActsParticlesPathInX0);
    acts_particles_tree->Branch("total_l0"               , &ActsParticlesPathInL0);
    acts_particles_tree->Branch("number_of_hits"         , &ActsParticlesNumberOfHits);
    acts_particles_tree->Branch("outcome"                , &ActsParticlesOutcome);
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
  flarHit = new TTree("flare_hits", "flarHitTreeInfo");

  flarHit->Branch("flareTrackID"			, &flareTrackID 		, "flareTrackID/I");
  flarHit->Branch("flareParticleID" 		, &flareParticleID		, "flareParticleID/I");
  flarHit->Branch("flareParentID" 			, &flareParentID		, "flareParentID/I");
  flarHit->Branch("flarePDG" 				, &flarePDG				, "flarePDG/I");
  flarHit->Branch("flareCopyNum" 			, &flareCopyNum			, "flareCopyNum/I");

  flarHit->Branch("flareT" 					, &flareT				, "flareT/I");

  flarHit->Branch("flareX" 					, &flareX				, "flareX/D");//Pre-position
  flarHit->Branch("flareY" 					, &flareY				, "flareY/D");
  flarHit->Branch("flareZ" 					, &flareZ				, "flareZ/D");

  flarHit->Branch("flarePx" 				, &flarePx				, "flarePx/D");//momentum
  flarHit->Branch("flarePy" 				, &flarePy				, "flarePy/D");
  flarHit->Branch("flarePz" 				, &flarePz				, "flarePz/D");

  flarHit->Branch("flareDeltaPx" 			, &flareDeltaPx			, "flareDeltaPx/D");
  flarHit->Branch("flareDeltaPy" 			, &flareDeltaPy			, "flareDeltaPy/D");
  flarHit->Branch("flareDeltaPz" 			, &flareDeltaPz			, "flareDeltaPz/D");

  flarHit->Branch("flareEdep" 				, &flareEdep			, "flareEdep/D");

}


void AnalysisManager::bookFlarHadCalHitTree() {
  flarHadHit = new TTree("hadXY_hits", "hadXYHitTreeInfo");

  //x
  flarHadHit->Branch("hadXTrackID"			, &hadXTrackID 			, "hadXTrackID/I");
  flarHadHit->Branch("hadXParticleID" 		, &hadXParticleID		, "hadXParticleID/I");
  flarHadHit->Branch("hadXParentID" 		, &hadXParentID			, "hadXParentID/I");
  flarHadHit->Branch("hadXPDG" 				, &hadXPDG				, "hadXPDG/I");
  flarHadHit->Branch("hadXCopyNum" 			, &hadXCopyNum			, "hadXCopyNum/I");

  flarHadHit->Branch("hadXT" 				, &hadXT				, "hadXT/I");

  flarHadHit->Branch("hadXx" 				, &hadXx				, "hadXx/D");//Pre-position
  flarHadHit->Branch("hadXy" 				, &hadXy				, "hadXy/D");
  flarHadHit->Branch("hadXz" 				, &hadXz				, "hadXz/D");

  flarHadHit->Branch("hadXPx" 				, &hadXPx				, "hadXPx/D");//momentum
  flarHadHit->Branch("hadXPy" 				, &hadXPy				, "hadXPy/D");
  flarHadHit->Branch("hadXPz" 				, &hadXPz				, "hadXPz/D");

  flarHadHit->Branch("hadXDeltaPx" 			, &hadXDeltaPx			, "hadXDeltaPx/D");
  flarHadHit->Branch("hadXDeltaPy" 			, &hadXDeltaPy			, "hadXDeltaPy/D");
  flarHadHit->Branch("hadXDeltaPz" 			, &hadXDeltaPz			, "hadXDeltaPz/D");

  flarHadHit->Branch("hadXEdep" 			, &hadXEdep				, "hadXEdep/D");

  //y
  flarHadHit->Branch("hadYTrackID"			, &hadYTrackID 			, "hadYTrackID/I");
  flarHadHit->Branch("hadYParticleID" 		, &hadYParticleID		, "hadYParticleID/I");
  flarHadHit->Branch("hadYParentID" 		, &hadYParentID			, "hadYParentID/I");
  flarHadHit->Branch("hadYPDG" 				, &hadYPDG				, "hadYPDG/I");
  flarHadHit->Branch("hadYCopyNum" 			, &hadYCopyNum			, "hadYCopyNum/I");

  flarHadHit->Branch("hadYT" 				, &hadYT				, "hadYT/I");

  flarHadHit->Branch("hadYx" 				, &hadYx				, "hadYx/D");//Pre-position
  flarHadHit->Branch("hadYy" 				, &hadYy				, "hadYy/D");
  flarHadHit->Branch("hadYz" 				, &hadYz				, "hadYz/D");

  flarHadHit->Branch("hadYPx" 				, &hadYPx				, "hadYPx/D");//momentum
  flarHadHit->Branch("hadYPy" 				, &hadYPy				, "hadYPy/D");
  flarHadHit->Branch("hadYPz" 				, &hadYPz				, "hadYPz/D");

  flarHadHit->Branch("hadYDeltaPx" 			, &hadYDeltaPx			, "hadYDeltaPx/D");
  flarHadHit->Branch("hadYDeltaPy" 			, &hadYDeltaPy			, "hadYDeltaPy/D");
  flarHadHit->Branch("hadYDeltaPz" 			, &hadYDeltaPz			, "hadYDeltaPz/D");

  flarHadHit->Branch("hadYEdep" 			, &hadYEdep				, "hadYEdep/D");


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
		bookFlarHadCalHitTree();
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
	flarHadHit->Write();
	thefile->cd();
  }
  if(m_saveTrack) trk->Write();
  if (m_saveActs) acts_hits_tree->Write();
  if (m_saveActs) acts_particles_tree->Write();
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

	hadXTrackID = 0;
	hadXParentID = 0;
	hadXPDG = 0;
	hadXCopyNum = 0;
	hadXParticleID = 0;
	hadXT = 0;
	hadXx = 0;
	hadXy = 0;
	hadXz = 0;
	hadXPx = 0;
	hadXPy = 0;
	hadXPz = 0;
	hadXDeltaPx = 0;
	hadXDeltaPy = 0;
	hadXDeltaPz = 0;
	hadXEdep = 0;

	hadYTrackID = 0;
	hadYParentID = 0;
	hadYPDG = 0;
	hadYCopyNum = 0;
	hadYParticleID = 0;
	hadYT = 0;
	hadYx = 0;
	hadYy = 0;
	hadYz = 0;
	hadYPx = 0;
	hadYPy = 0;
	hadYPz = 0;
	hadYDeltaPx = 0;
	hadYDeltaPy = 0;
	hadYDeltaPz = 0;
	hadYEdep = 0;


  }

  ActsHitsEventID = 0;
  ActsHitsGeometryID = 0;
  ActsHitsParticleID = 0;
  ActsHitsX = 0;
  ActsHitsY = 0;
  ActsHitsZ = 0;
  ActsHitsT = 0;
  ActsHitsPx = 0;
  ActsHitsPy = 0;
  ActsHitsPz = 0;
  ActsHitsE = 0;
  ActsHitsDeltaPx = 0;
  ActsHitsDeltaPy = 0;
  ActsHitsDeltaPz = 0;
  ActsHitsDeltaE = 0;
  ActsHitsIndex = 0;
  ActsHitsVolumeID = 0;
  ActsHitsBoundaryID = 0;
  ActsHitsLayerID = 0;
  ActsHitsApproachID = 0;
  ActsHitsSensitiveID = 0;

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

		//G4cout<<"passed1"<<G4endl;

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
  // update number of primary particles after FillPrimaryTruthTree
  // including decay products from primary tau and pizero
  nPrimaryParticle = countPrimaryParticle; // TODO: they're equivalent, should only store one of them
  nFromFSLParticles = tracksFromFSLSecondary.size();
  nFromFSPizeroParticles = tracksFromFSPizeroSecondary.size();
  nFromFSLDecayPizeroParticles = tracksFromFSLDecayPizeroSecondary.size();
//G4cout<<"passed9"<<G4endl;

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
 //G4cout<<"passed10"<<G4endl;

  // tracksFromFSL includes all the tracks orginating from the fsl
  // tracksFromFSLSecondary only inclues the tracks directly decayed from the fsl
  /* std::cout<<"Recorded tracks       : "<<allTracksPTPair.size()<<std::endl; */
  /* std::cout<<"Tracks from FSL       : "<<tracksFromFSL.size()<<std::endl; */
  /* std::cout<<"Tracks from FSL (2nd) : "<<tracksFromFSLSecondary.size()<<std::endl; */
  /* std::cout<<"number of primary particles : "<<nPrimaryParticle */
  /*   <<" , in which number of particles from fsl : "<<nFromFSLParticles<<std::endl; */
  //std::cout<<"Test nTestNPrimaryTrack : "<<nTestNPrimaryTrack<<std::endl;
  // return;
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
//G4cout<<"passed2"<<G4endl;

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

  /*   if (m_circularFit){ */

  /*     circNhits = hitXFSL.size(); */
  /*     trkNhits = trkXFSL.size(); */

  /*     // apply circular fitting to FLArE hadCath/muonFinder */
  /*     if ( circNhits > 0 ){ */
  /*       circularfitter::CircleFit* circFit = new circularfitter::CircleFit(hitXFSL,hitZFSL); */
  /*       circularfitter::LineFit* lineFit = new circularfitter::LineFit(hitZFSL,hitYFSL); */
  /*       xc = circFit->GetXc(); */
  /*       zc = circFit->GetZc(); */
  /*       rc = circFit->GetR(); */
  /*       p0 = lineFit->GetP0(); */
  /*       p1 = lineFit->GetP1(); */
  /*       cosDip = lineFit->GetCosDip(); */
  /*     } */

  /*     // apply circular fitting for FASER2 spectrometer magnet */
  /*     if( trkNhits > 0 ){ */

  /*       Nmagnets = (GeometricalParameters::Get()->GetFASER2MagnetOption() == GeometricalParameters::magnetOption::SAMURAI) ? 1 : */
  /*                   GeometricalParameters::Get()->GetNFASER2Magnets(); */
  /*       G4cout << "Number of FASER2 magnets: " << Nmagnets << G4endl; */

  /*       circularfitter::CircleExtractor* circExtract = new circularfitter::CircleExtractor(trkXFSL,trkYFSL,trkZFSL); */
  /*       magzpos = circExtract->GetMagnetZs(); */
  /*       trkxc = circExtract->GetXc(); */
  /*       trkzc = circExtract->GetZc(); */
  /*       trkrc = circExtract->GetR(); */

  /*       std::vector<circularfitter::line> pre  = circExtract->GetPreLine(); */
  /*       std::vector<circularfitter::line> post = circExtract->GetPostLine(); */
  /*       for(int i=0; i<Nmagnets; i++){ */
          //G4cout << "circle: " << trkzc[i] << " " << trkxc[i] << " " << trkrc[i] << G4endl;
          /* trkqIn.push_back( pre[i].q ); */
          /* trkmIn.push_back( pre[i].m ); */
          /* trkqOut.push_back( post[i].q ); */
          /* trkmOut.push_back( post[i].m ); */
        /* } */

        /* circularfitter::LineFit* lineFit2 = new circularfitter::LineFit(trkZFSL,trkYFSL); */
        /* trkp0 = lineFit2->GetP0(); */
        /* trkp1 = lineFit2->GetP1(); */
        /* trkcosDip = lineFit2->GetCosDip(); */
      /* } */
    /* } */

    /* // FillPseudoRecoVar must run after FillTrueEdep, otherwise some of the variables won't be filled */
    /* FillPseudoRecoVar(); */

    /* delete pm3D; */
  /* } */
//G4cout<<"passed3"<<G4endl;

  /* int count_tracks = 0; */
  /* if(m_saveTrack){ */
  /*   G4cout << "---> Saving track information to tree..." << G4endl; */
  /*   auto trajectoryContainer = event->GetTrajectoryContainer(); */
  /*   if (trajectoryContainer) { */
  /*     for (size_t i = 0; i < trajectoryContainer->entries(); ++i) { */
  /*       auto trajectory = static_cast<G4Trajectory*>((*trajectoryContainer)[i]); */
  /*       trackTID = trajectory->GetTrackID(); */
  /*       trackPID = trajectory->GetParentID(); */
  /*       trackPDG = trajectory->GetPDGEncoding(); */
	/* trackKinE = trajectory->GetInitialKineticEnergy(); */
  /*       trackNPoints = trajectory->GetPointEntries(); */
  /*       count_tracks++; */
	/* for (size_t j = 0; j < trackNPoints; ++j) { */
  /*         G4ThreeVector pos = trajectory->GetPoint(j)->GetPosition(); */
	  /* trackPointX.push_back( pos.x() ); */
	  /* trackPointY.push_back( pos.y() ); */
	  /* trackPointZ.push_back( pos.z() ); */
  /*       } */
  /*       trk->Fill(); */
	/* trackPointX.clear(); */
	/* trackPointY.clear(); */
	/* trackPointZ.clear(); */
  /*     } */
  /*   } else G4cout << "No tracks found: did you enable their storage with '/tracking/storeTrajectory 1'?" << G4endl; */
  /*   G4cout << "---> Done!" << G4endl; */
  /* } */

  /* evt->Fill(); */

  /* G4cout << "Total number of recorded hits : " << nHits << std::endl; */
  /* if(m_saveTrack) G4cout << "Total number of recorded track: " << count_tracks << std::endl; */

  /* for (int iPrim= 0; iPrim< nPrimaryParticle; ++iPrim) { */
  /*   trackClusters[iPrim].clear(); */
  /* } */
  /* trackClusters.clear(); */
  /* trackClusters.shrink_to_fit(); */

}

void AnalysisManager::FillPrimaryTruthTree(G4int sdId, std::string sdName)
{
  //G4cout<<"passed4"<<G4endl;
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

	  double time = hit->GetTime();

	  double Px = hit->GetInitMomentum().x();
  	  double Py = hit->GetInitMomentum().y();
	  double Pz = hit->GetInitMomentum().z();

	  double dPx = hit->GetDeltaMom().x();
  	  double dPy = hit->GetDeltaMom().y();
	  double dPz = hit->GetDeltaMom().z();


	  double edep =  hit->GetEdep();

	  auto particleId = ActsFatras::Barcode();
	  particleId.setVertexPrimary(1);//fix this value
	  particleId.setGeneration(pid);
	  particleId.setSubParticle(0);
	  particleId.setParticle(tid);
	  double partID = particleId.value();


      // energy deposition in different volumes of the detector
      if (sdName == "lArBoxSD/lar_box"){
	  	flareTrackID = tid;
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

      else if (sdName == "HadCalXSD/lar_box"){
		G4cout<<"hadX hit"<<G4endl;

	  	hadXTrackID = tid;
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

	  	flarHadHit->Fill();
		//G4cout<<"passed5.5"<<G4endl;
	  }

		else if (sdName == "HadCalXSD/lar_box"){
		G4cout<<"hadX hit"<<G4endl;

	  	hadYTrackID = tid;
	  	hadYParentID = pid;
	  	hadYPDG = PDG;
	  	hadYCopyNum = copyNum;
		hadYParticleID = partID;
		hadYT = time;
	  	hadYx = pre_x;
	  	hadYy = pre_y;
	  	hadYz = pre_z;
	  	hadYPx = Px;
	  	hadYPy = Py;
	  	hadYPz = Pz;
	  	hadYDeltaPx = dPx;
	  	hadYDeltaPy = dPy;
	  	hadYDeltaPz = dPz;
	  	hadYEdep = edep;

	  	flarHadHit->Fill();
		//G4cout<<"passed5.4"<<G4endl;

	  }

    }
  }

  if (m_saveActs)
  {
    auto hitCollection = dynamic_cast<FASER2TrackerHitsCollection*>(hcofEvent->GetHC(sdId));
    if (!hitCollection) return;
    std::map<G4int, G4int> sub_part_map{};
    for (auto hit: *hitCollection->GetVector())
    {
      if (hit->GetCharge() == 0) continue; // skip neutral particles, they don't hit

      /*
      * A note on the ActsHitsGeometryID variable
       This variable in Acts keeps track of an Acts::GeometryIdentifier. This is essentially a long unsigned int, the bits of which are used to
       look up the the volume/layer/boundary/sensitive indices of a piece of geometry. In principle it should be possible to assign this variable
       here in GEANT4 but I don't understand the Acts code well enough to do it without adding Acts as a dependancy to this codebase.
       As a result I set `geometry_id` to zero and give the the resposibility of assigning this variable to the user during the reading of the `hits` tree.
      */

      ActsHitsEventID = evtID;
      ActsHitsGeometryID = 0;

      int hitID = hit->GetTrackID();
      int nPrimaries = ActsParticlesParticleId.size();

      auto particleId = ActsFatras::Barcode();
      particleId.setVertexPrimary(1);
      particleId.setVertexSecondary(0);
      particleId.setParticle(hit->GetTrackID() - 1); // The track ID is the primary particle index plus one
      particleId.setGeneration(hit->GetParentID());

      sub_part_map.try_emplace(hit->GetTrackID()-1, sub_part_map.size());
	  //G4cout<<"passed6"<<G4endl;

      // This is a fudge - assumes that that the secondary particles are always sub-particles of the primary particle
      particleId.setSubParticle(hit->GetParentID() == 0 ? 0 : sub_part_map[hit->GetTrackID()-1]);
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
      ActsHitsLayerID = (hit->GetCopyNumSensor()+1)*2; // Acts specfic layer ID, goes 2, 4, 6, ...
      ActsHitsApproachID = 0;
      ActsHitsSensitiveID = 1;
      acts_hits_tree->Fill();

      // Now fill the Acts particles tree
      bool isDuplicate = false;
      for (const auto& id : ActsParticlesParticleId) {
        if (id == particleId.value()) {
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
      ActsParticlesVertexPrimary.push_back(hit->GetIsPrimaryTrack()); //? These variables need to be filled, but are unused by Acts
      ActsParticlesVertexSecondary.push_back(hit->GetIsSecondaryTrack()); //? These variables need to be filled, but are unused by Acts
      ActsParticlesParticle.push_back(1); //? These variables need to be filled, but are unused by Acts
      ActsParticlesGeneration.push_back(0); //? These variables need to be filled, but are unused by Acts
      ActsParticlesSubParticle.push_back(0); //? These variables need to be filled, but are unused by Acts
      ActsParticlesELoss.push_back(0); //? These variables need to be filled, but are unused by Acts
      ActsParticlesPathInX0.push_back(0); //? These variables need to be filled, but are unused by Acts
      ActsParticlesPathInL0.push_back(0); //? These variables need to be filled, but are unused by Acts
      ActsParticlesNumberOfHits.push_back(0); //? These variables need to be filled, but are unused by Acts
      ActsParticlesOutcome.push_back(0); //? These variables need to be filled, but are unused by Acts
    } // end of loop over hits
    acts_particles_tree->Fill();
	//G4cout<<"passed7"<<G4endl;

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

  //     if (m_saveHit) {
  //       if (nHits<=40000000) {
  //         HitTID[nHits-1] = hit->GetTID();
  //         HitPID[nHits-1] = hit->GetPID();
  //         HitPDG[nHits-1] = hit->GetParticle();
  //         HitTrackStatus[nHits-1]  = hit->GetTrackStatus();
  //         HitPrePositionX[nHits-1] = pre_x;
  //         HitPrePositionY[nHits-1] = pre_y;
  //         HitPrePositionZ[nHits-1] = pre_z;
  //         HitPosPositionX[nHits-1] = post_x;
  //         HitPosPositionY[nHits-1] = post_y;
  //         HitPosPositionZ[nHits-1] = post_z;
  //         HitEdep[nHits-1] = hit->GetEdep();
  //       }
  //     }

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

  //     // save FSL (only muons!) hits for circle fitting
  //     // check that particle is a primary muon
  //     if( TMath::Abs(hit->GetParticle())==13 && hit->GetPID() == 0 ){
  //       double px = hit->GetInitMomentum().x();
  //       double pz = hit->GetInitMomentum().z();
  //       double p_perp = TMath::Sqrt(px*px+pz*pz);
  //       // save hits in FLArE hadCather+muonFinder
  //       if( (sdName == "HadCalXSD/lar_box") || (sdName == "HadCalYSD/lar_box") ||
  //           (sdName == "MuonFinderXSD/lar_box") || (sdName == "MuonFinderYSD/lar_box") ||
  //           (sdName == "HadAbsorbSD/lar_box") || (sdName == "MuonFinderAbsorbSD/lar_box") ||
  //           (sdName == "BabyMINDHorBarSD/lar_box") || (sdName == "BabyMINDVerBarSD/lar_box")){
  //         hitXFSL.push_back(post_x);
  //         hitYFSL.push_back(post_y);
  //         hitZFSL.push_back(post_z);
  //         hitPFSL.push_back(p_perp);
  //       }
  //       else if ((sdName == "TrkHorScinSD/lar_box") || (sdName == "TrkVerScinSD/lar_box")){
  //         trkXFSL.push_back(post_x);
  //         trkYFSL.push_back(post_y);
  //         trkZFSL.push_back(post_z);
  //         trkPFSL.push_back(p_perp);
  //       }
  //     }

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
}

void AnalysisManager::FillTrueEdep(G4int sdId, std::string sdName)
{
}

double AnalysisManager::GetTotalEnergy(double px, double py, double pz, double m) {
  return TMath::Sqrt(px*px+py*py+pz*pz+m*m);
}

void AnalysisManager::FillPseudoRecoVar() {
  //G4cout<<"passed8"<<G4endl;

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
