#include <vector>
#include <functional>
#include <iostream>
#include <string>
#include <sstream>
#include <map>
#include <iomanip>
#include <random>

#include <G4Event.hh>
#include <G4SDManager.hh>
#include <G4SystemOfUnits.hh>
#include <Randomize.hh>
#include <G4Poisson.hh>
#include <G4LorentzVector.hh>
#include <G4Exception.hh>

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
#include "EventInformation.hh"
#include "generators/GeneratorVertexMetadata.hh"
#include "reco/PCAAnalysis3D.hh"
#include "reco/Cluster3D.hh"
#include "reco/LinearFit.hh"
#include "reco/ShowerLID.hh"
#include "reco/Barcode.hh"
#include "FPFParticle.hh"
#include "FPFTrajectory.hh"

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
  fPar = nullptr;
  fFLArEHits = nullptr;
  fFLArEHCALHits = nullptr;
  fFLArEPseudoReco = nullptr;
  fActsHitsTree = nullptr;

  fSaveTrack = false;
  fParKinECut = 50.0*keV;
  fSave3DEvd = false;
  fSave2DEvd = false;
  fSavePseudoReco = false;
  fEnableFLArE = true;
  fSaveActs = true;
  fAddDiffusion = "false";
}

AnalysisManager::~AnalysisManager() {}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::bookEvtTree()
{
  fEvt = new TTree("event", "event info");
  fEvt->Branch("event_id", &evtID, "event_id/I");
  fEvt->Branch("vertex_id", &vertexID, "vertex_id/I");
  fEvt->Branch("weight", &weight, "weight/D");
  fEvt->Branch("generatorType", &genType);
  fEvt->Branch("processName", &processName);
  fEvt->Branch("vtxX", &vtxX, "vtxX/D");
  fEvt->Branch("vtxY", &vtxY, "vtxY/D");
  fEvt->Branch("vtxZ", &vtxZ, "vtxZ/D");
  fEvt->Branch("vtxT", &vtxT, "vtxT/D");
  fEvt->Branch("initPDG", &initPDG, "initPDG/I");
  fEvt->Branch("initPx", &initPx, "initPx/D");
  fEvt->Branch("initPy", &initPy, "initPy/D");
  fEvt->Branch("initPz", &initPz, "initPz/D");
  fEvt->Branch("initE", &initE, "initE/D");
  fEvt->Branch("initM", &initM, "initM/D");
  fEvt->Branch("initQ", &initQ, "initQ/D");
  fEvt->Branch("tgtPDG", &tgtPDG, "tgtPDG/I");
  fEvt->Branch("tgtA", &tgtA, "tgtA/I");
  fEvt->Branch("tgtZ", &tgtZ, "tgtZ/I");
  fEvt->Branch("hitnucPDG", &hitnucPDG, "hitnucPDG/I");
  fEvt->Branch("fslPDG", &fslPDG, "fslPDG/I");
  fEvt->Branch("intType", &intType, "intType/I");
  fEvt->Branch("scatteringType", &scatteringType, "scatteringType/I");
  fEvt->Branch("xs", &xs, "xs/D");
  fEvt->Branch("Q2", &Q2, "Q2/D");
  fEvt->Branch("xBj", &xBj, "xBj/D");
  fEvt->Branch("y", &y, "y/D");
  fEvt->Branch("W", &W, "W/D");
  fEvt->Branch("nPrimaries", &nPrimaries, "nPrimaries/I");
  fEvt->Branch("primTID", &primTID);
  fEvt->Branch("primPDG", &primPDG);
  fEvt->Branch("primPx", &primPx);
  fEvt->Branch("primPy", &primPy);
  fEvt->Branch("primPz", &primPz);
  fEvt->Branch("primE", &primE);
}

void AnalysisManager::bookParTree()
{
  fPar = new TTree("particles", "particle info");
  fPar->Branch("event_id", &evtID, "event_id/I");
  fPar->Branch("particle_id", &particle_id, "particle_id/I");
  fPar->Branch("particle_pdg", &particle_PDG, "particle_pdg/I");
  fPar->Branch("track_id", &particle_TID, "track_id/I");
  fPar->Branch("parent_id", &particle_PID, "parent_id/I");
  fPar->Branch("ancestor_id", &particle_ancestor, "ancestor_id/I");
  fPar->Branch("process", &particle_process);
  fPar->Branch("vx", &particle_vx, "vx/F");
  fPar->Branch("vy", &particle_vy, "vy/F");
  fPar->Branch("vz", &particle_vz, "vz/F");
  fPar->Branch("vt", &particle_vt, "vt/F");
  fPar->Branch("px", &particle_px, "px/F");
  fPar->Branch("py", &particle_py, "py/F");
  fPar->Branch("pz", &particle_pz, "pz/F");
  fPar->Branch("m", &particle_m, "m/F");
  fPar->Branch("q", &particle_q, "q/F");
  fPar->Branch("eta", &particle_eta, "eta/F");
  fPar->Branch("phi", &particle_phi, "phi/F");
  fPar->Branch("pt", &particle_pt, "pt/F");
  fPar->Branch("p", &particle_p, "p/F");
  fPar->Branch("ke", &particle_ke, "ke/F");
  if (fSaveTrack) // if saving full trajectory, add vector of trajectory points
  {
    fPar->Branch("traj_Npoints", &traj_Npoints, "traj_Npoints/I");
    fPar->Branch("traj_pointX", &traj_pointX);
    fPar->Branch("traj_pointY", &traj_pointY);
    fPar->Branch("traj_pointZ", &traj_pointZ);
  }
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::bookFLArETrees()
{
  // create subdirectory in file
  fFLArEDir = fFile->mkdir("flare", "FLArE output", kTRUE);
  fFile->cd(fFLArEDir->GetName());

  // FLArE hits tree
  fFLArEHits = new TTree("flare_hits", "FLArE hits info");
  fFLArEHits->Branch("flareEvtID", &evtID, "flareEvtID/I");
  fFLArEHits->Branch("flareTrackID", &flareTrackID, "flareTrackID/I");
  fFLArEHits->Branch("flareBarcode", &flareParticleID, "flareParticleID/I");
  fFLArEHits->Branch("flareParentID", &flareParentID, "flareParentID/I");
  fFLArEHits->Branch("flareAncestorID", &flareAncestorID, "flareAncestorID/I");
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
  fFLArEHCALHits->Branch("hadAncestorID", &flareAncestorID, "hadAncestorID/I");
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

  // FLArE pseudo-reco
  if (fSavePseudoReco)
  {
    fFLArEPseudoReco = new TTree("pseudo_reco", "FLArE pseudo-reco info");
    fFLArEPseudoReco->Branch("evtID", &evtID, "evtID/I");
    fFLArEPseudoReco->Branch("edepInLAr", &edepInLAr, "edepInLAr/F");
    fFLArEPseudoReco->Branch("edepInHCAlX", &edepInHCALX, "edepInHCALX/F");
    fFLArEPseudoReco->Branch("edepInHCALY", &edepInHCALY, "edepInHCALY/F");
    fFLArEPseudoReco->Branch("sparseFractionMem", &sparseFractionMem, "sparseFractionMem/F");
    fFLArEPseudoReco->Branch("sparseFractionBins", &sparseFractionBins, "sparseFractionBins/F");
    fFLArEPseudoReco->Branch("TotalDedxLongitudinal", &TotalDedxLongitudinal);
    fFLArEPseudoReco->Branch("TrueTotalDedxLongitudinal", &TrueTotalDedxLongitudinal);
    fFLArEPseudoReco->Branch("nPrimaryParticle", &nprimaries, "nPrimaryParticle/I");
    fFLArEPseudoReco->Branch("primaryPDG", &primaryPDG);
    fFLArEPseudoReco->Branch("primaryTrackLength", &primaryTrackLength);
    fFLArEPseudoReco->Branch("primaryTrackLengthInTPC", &primaryTrackLengthInTPC);
    fFLArEPseudoReco->Branch("ProngEInLAr", &ProngEInLAr);
    fFLArEPseudoReco->Branch("ProngEInHadCal", &ProngEInHadCal);
    fFLArEPseudoReco->Branch("ProngAngleToBeamDir", &ProngAngleToBeamDir);
    fFLArEPseudoReco->Branch("ShowerLength", &ShowerLength);
    fFLArEPseudoReco->Branch("ShowerLengthInLAr", &ShowerLengthInLAr);
    fFLArEPseudoReco->Branch("ShowerWidth", &ShowerWidth);
    fFLArEPseudoReco->Branch("ShowerWidthInLAr", &ShowerWidthInLAr);
    fFLArEPseudoReco->Branch("ProngAvgdEdx", &ProngAvgdEdx);
    fFLArEPseudoReco->Branch("ProngAvgdEdxInLAr", &ProngAvgdEdxInLAr);
    fFLArEPseudoReco->Branch("dir_pol_x", &dir_pol_x);
    fFLArEPseudoReco->Branch("dir_pol_y", &dir_pol_y);
    fFLArEPseudoReco->Branch("dir_pol_z", &dir_pol_z);
    fFLArEPseudoReco->Branch("dir_coc_x", &dir_coc_x);
    fFLArEPseudoReco->Branch("dir_coc_y", &dir_coc_y);
    fFLArEPseudoReco->Branch("dir_coc_z", &dir_coc_z);
  }

  fFile->cd();
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::bookFASER2Trees()
{
  // create subdirectory in file
  fFASER2Dir = fFile->mkdir("faser2", "FASER2 output", kTRUE);
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

  fFile->cd();
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::BeginOfRun()
{
  G4cout << "Run has been started, preparing output" << G4endl;

  if (fFile)
    delete fFile;

  // Preparing output file
  fFile = new TFile(fFilename.c_str(), "RECREATE");

  // Booking common output trees
  bookEvtTree();
  bookParTree();

  fSDNamelist = GeometricalParameters::Get()->GetSDNamelist();
  G4cout << "Number of SDs : " << fSDNamelist.size() << G4endl;

  // Loop through active sentive detectors
  // book trees as necessary
  for (auto sdname : fSDNamelist)
  {
    G4cout << sdname.first << " " << sdname.second << G4endl;

    // if FLArE is enabled, book trees
    if (fEnableFLArE && sdname.second.find("FLArE") != std::string::npos)
    {
      fFlareSDs.push_back(sdname.first);
      bookFLArETrees();
    }

    // if FASER2 is enabled, book ACTS trees
    // but give the option to disable output if required
    else if (fSaveActs && sdname.second.find("FASER2") != std::string::npos)
    {
      fFaser2SDs.push_back(sdname.first);
      bookFASER2Trees();
    }
  }

  // if FLArE is enabled & its geometry was found
  // prepare .h5 output file
  if (fFlareSDs.size() > 0)
  {
    fH5Filename = fFilename;
    if (fH5Filename.find(".root") != std::string::npos)
    {
      const size_t pos = fH5Filename.find(".root");
      fH5Filename.resize(pos);
    }
    fH5Filename += ".h5";
    fH5file = hep_hpc::hdf5::File(fH5Filename, H5F_ACC_TRUNC);
  }
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::EndOfRun()
{
  // save common trees at the top of the output file
  fFile->cd();
  fEvt->Write();
  fPar->Write();

  // save detector-specif trees in their directories
  if (fFlareSDs.size() > 0)
  {
    fFile->cd(fFLArEDir->GetName());
    fFLArEHits->Write();
    fFLArEHCALHits->Write();
    if (fSavePseudoReco)
      fFLArEPseudoReco->Write();
    fFile->cd(); // go back to top
    fH5file.close();
  }
  if (fFaser2SDs.size() > 0)
  {
    fFile->cd(fFASER2Dir->GetName());
    fActsHitsTree->Write();
    fFile->cd(); // go back to top
  }

  fFile->Close();
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::BeginOfEvent()
{
  // reset vectors that need to be cleared for a new event
  // only reset arrays or vectors, tipically no need for other defaults

  // primaries in event tree
  // (actually need reset at every vertex)
  primaries.clear();
  primaryIDs.clear();
  primTID.clear();
  primPDG.clear();
  primPx.clear();
  primPy.clear();
  primPz.clear();
  primE.clear();

  // track ID to primary ancestor association
  trackToPrimaryAncestor.clear();

  // track ID to particle ID association
  trackIDtoParticleID.clear();

  // trajectory points (if enabled)
  // (actually need reset at every track)
  traj_pointX.clear();
  traj_pointY.clear();
  traj_pointZ.clear();

  // these are used to accumulate
  // so need to be reset
  edepInLAr = 0;
  edepInHCALX = 0;
  edepInHCALY = 0;

  TotalDedxLongitudinal.clear();
  TrueTotalDedxLongitudinal.clear();
  TotalDedxLongitudinal.resize(3000, 0.0);
  ;
  TrueTotalDedxLongitudinal.resize(3000, 0.0);

  primaryPDG.clear();
  primaryTrackLength.clear();
  primaryTrackLengthInTPC.clear();
  ProngEInLAr.clear();
  ProngEInHadCal.clear();
  ProngAngleToBeamDir.clear();
  ShowerLength.clear();
  ShowerLengthInLAr.clear();
  ShowerWidth.clear();
  ShowerWidthInLAr.clear();
  ProngAvgdEdx.clear();
  ProngAvgdEdxInLAr.clear();
  dir_pol_x.clear();
  dir_pol_y.clear();
  dir_pol_z.clear();
  dir_coc_x.clear();
  dir_coc_y.clear();
  dir_coc_z.clear();
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::EndOfEvent(const G4Event *event)
{

  /// evtID
  evtID = event->GetEventID();

  //-----------------------------------------------------------
  // FILL EVENT TREE
  FillEventTree(event);

  //-----------------------------------------------------------
  // FILL PARTICLES/TRAJECTORIES TREE
  FillParticlesTree(event);

  //-----------------------------------------------------------
  // FILL DETECTOR HITS

  // Get the hit collections
  // If there is no hit collection, there is nothing to be done
  fHCofEvent = event->GetHCofThisEvent();
  if (!fHCofEvent)
  {
    G4cout << "No hits recorded in any sensitive volume --> nothing to save!" << G4endl;
    return;
  }

  if (fFlareSDs.size() > 0)
    FillFLArEOutput();

  if (fFaser2SDs.size() > 0)
    FillFASER2Output();
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::FillEventTree(const G4Event *event)
{
  EventInformation *eventInfo = static_cast<EventInformation *>(event->GetUserInformation());
  auto metadata = eventInfo->GetEventMetadata();

  int nVertexes = metadata.size();
  int nG4Vertexes = event->GetNumberOfPrimaryVertex();
  if (nVertexes != nG4Vertexes)
  {
    std::ostringstream oss;
    oss << "Mismatch between vertexes in generator metadata vs G4Event (" << nVertexes << " vs " << nG4Vertexes << ")";
    G4Exception("AnalysisManager", "LogicError", JustWarning, oss.str().c_str());
  }

  /// loop over the vertices, and then over primary particles,
  /// metadata info comes event generator
  /// primaries come from G4event vertex
  for (int ivtx = 0; ivtx < nVertexes; ivtx++)
  {
    eventInfo->Print(ivtx);

    vertexID = ivtx;
    weight = metadata[ivtx].weight;
    genType = metadata[ivtx].generatorType;
    processName = metadata[ivtx].processName;
    vtxX = metadata[ivtx].x4.x();
    vtxY = metadata[ivtx].x4.y();
    vtxZ = metadata[ivtx].x4.z();
    vtxT = metadata[ivtx].x4.t();
    initPDG = metadata[ivtx].pdg;
    initPx = metadata[ivtx].p4.x();
    initPy = metadata[ivtx].p4.y();
    initPz = metadata[ivtx].p4.z();
    initE = metadata[ivtx].p4.e();
    initM = metadata[ivtx].mass;
    initQ = metadata[ivtx].charge;
    fslPDG = metadata[ivtx].fsl_pdg;
    tgtPDG = metadata[ivtx].tgt_pdg;
    tgtZ = metadata[ivtx].tgt_Z;
    tgtA = metadata[ivtx].tgt_A;
    hitnucPDG = metadata[ivtx].hitnuc_pdg;
    intType = metadata[ivtx].intType;
    scatteringType = metadata[ivtx].scatteringType;
    xs = metadata[ivtx].xs;
    Q2 = metadata[ivtx].Q2;
    xBj = metadata[ivtx].xBj;
    y = metadata[ivtx].y;
    W = metadata[ivtx].W;

    nPrimaries = event->GetPrimaryVertex(ivtx)->GetNumberOfParticle();
    primTID.clear();
    primPDG.clear();
    primPx.clear();
    primPy.clear();
    primPz.clear();
    primE.clear();

    G4cout << "\n-> number of primaries: " << nPrimaries << G4endl;

    for (int ipp = 0; ipp < event->GetPrimaryVertex(ivtx)->GetNumberOfParticle(); ++ipp)
    {
      G4PrimaryParticle *primary_particle = event->GetPrimaryVertex(ivtx)->GetPrimary(ipp);
      if (primary_particle)
      {
        int tid = ipp + 1; // confirm matches track id?
        int pdg = primary_particle->GetPDGcode();
        double px = primary_particle->GetMomentum().x();
        double py = primary_particle->GetMomentum().y();
        double pz = primary_particle->GetMomentum().z();
        double m = primary_particle->GetMass() / MeV;
        double e = GetTotalEnergy(px, py, pz, m);

        primTID.push_back(tid);
        primPDG.push_back(pdg);
        primPx.push_back(px);
        primPy.push_back(py);
        primPz.push_back(pz);
        primE.push_back(e);

        // store a copy as a FPFParticle for further processing
        // this accumulates across entire event (no vertex reset)
        primaryIDs.push_back(tid); // store to avoid duplicates
        primaries.push_back(FPFParticle(pdg, 0,
                                        tid, primaryIDs.size() - 1, 1,
                                        m, vtxX, vtxY, vtxZ, vtxT,
                                        px, py, pz, e));

        G4cout << "TID: " << tid << ", PDG: " << pdg << ", p=(" << px << ", " << py << ", " << pz << ") MeV" << G4endl;
      }
    }
    fEvt->Fill();
  }
  G4cout << "\nTotal number of vertexes  : " << nVertexes << G4endl;
  G4cout << "Total number of primaries  : " << primaryIDs.size() << G4endl;
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::FillParticlesTree(const G4Event *event)
{
  int count_tracks = 0;
  G4cout << "==== Saving particle information to tree ====" << G4endl;
  auto trajectoryContainer = event->GetTrajectoryContainer();
  if (!trajectoryContainer)
  {
    G4cout << "No trajectories found in the event!" << G4endl;
    return;
  }

  std::map<G4int, G4int> sub_part_map{};
  for (size_t i = 0; i < trajectoryContainer->entries(); ++i)
  {
    auto trajectory = static_cast<FPFTrajectory *>((*trajectoryContainer)[i]);

    particle_TID = trajectory->GetTrackID();
    particle_PID = trajectory->GetParentID();
    particle_PDG = trajectory->GetPDGEncoding();
    particle_ancestor = GetTrackPrimaryAncestor(particle_TID);
    particle_process = trajectory->GetProcessName();

    auto particleId = ActsFatras::Barcode();
    particleId.setVertexPrimary(1);
    particleId.setVertexSecondary(0);
    particleId.setParticle(particle_TID - 1); // The track ID is the primary particle index plus one
    particleId.setGeneration(particle_PID);
    sub_part_map.try_emplace(particle_TID- 1, sub_part_map.size());
    // This is a fudge - assumes that the secondary particles are always sub-particles of the primary particle
    particleId.setSubParticle(particle_PID == 0 ? 0 : sub_part_map[particle_TID - 1]);
    particle_id = particleId.value();

    // store in map for use in other trees!
    trackIDtoParticleID[particle_TID] = particle_id;

    // first point is vertex (always stored even if traj is disabled)
    particle_vx = trajectory->GetPoint(0)->GetPosition().x();
    particle_vy = trajectory->GetPoint(0)->GetPosition().y();
    particle_vz = trajectory->GetPoint(0)->GetPosition().z();
    particle_vt = 0; 
    particle_px = trajectory->GetInitialP4().px();
    particle_py = trajectory->GetInitialP4().py();
    particle_pz = trajectory->GetInitialP4().pz();
    particle_m = trajectory->GetInitialP4().m();
    particle_q = trajectory->GetCharge();
    particle_eta = trajectory->GetInitialP4().eta();
    particle_phi = trajectory->GetInitialP4().phi();
    particle_pt = trajectory->GetInitialP4().vect().perp();
    particle_p = trajectory->GetInitialP4().vect().mag();
    
    // do not save if below cut threshold
    particle_ke = trajectory->GetInitialKineticEnergy();
    if(particle_ke < fParKinECut) continue;
  
    traj_Npoints = trajectory->GetPointEntries();
    traj_pointX.clear();
    traj_pointY.clear();
    traj_pointZ.clear();

    for (size_t j = 0; j < traj_Npoints; ++j)
    {
      G4ThreeVector pos = trajectory->GetPoint(j)->GetPosition();
      traj_pointX.push_back(pos.x());
      traj_pointY.push_back(pos.y());
      traj_pointZ.push_back(pos.z());
    }
    count_tracks++;
    fPar->Fill();
  }

  G4cout << "Total number of recorded particles: " << count_tracks << std::endl;
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::FillFLArEOutput()
{
  G4cout << "===== Filling FLArE output trees =====" << G4endl;

  // prepare the pixel map
  const double_t res_tpc[3] = {1, 5, 5}; // mm

  pm3D = new PixelMap3D(evtID, primaries.size(), initPDG, res_tpc);
  // boundaries in global coordinates
  pm3D->SetPMBoundary(GeometricalParameters::Get()->GetFLArEPosition() / mm -
                          GeometricalParameters::Get()->GetTPCSizeXYZ() / mm / 2,
                      GeometricalParameters::Get()->GetFLArEPosition() / mm +
                          GeometricalParameters::Get()->GetTPCSizeXYZ() / mm / 2);
  pm3D->InitializePM();

  // accumulate values per primary particle
  ProngEInLAr.resize(primaries.size(), 0.);
  ProngEInHadCal.resize(primaries.size(), 0.);
  ShowerLength.resize(primaries.size(), 0.);
  ShowerLengthInLAr.resize(primaries.size(), 0.);
  ShowerWidth.resize(primaries.size(), 0.);
  ShowerWidthInLAr.resize(primaries.size(), 0.);
  primaryTrackLength.resize(primaries.size(), 0.);
  primaryTrackLengthInTPC.resize(primaries.size(), 0.);

  // loop over the detected FLArE sensitive volumes
  // for now, everything is lar_box
  int nLArHits = 0;
  int nHCALHits = 0;
  for (const int sdId : fFlareSDs)
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
      flareTrackID = hit->GetTID();
      flareParentID = hit->GetPID();
      flareAncestorID = GetTrackPrimaryAncestor(flareTrackID);
      flarePDG = hit->GetPDG();
      flareCopyNum = hit->GetCopyNum();
      flareT = hit->GetTime();
      flareX = hit->GetPostStepPosition().x();
      flareY = hit->GetPostStepPosition().y();
      flareZ = hit->GetPostStepPosition().z();
      flarePx = hit->GetInitMomentum().x();
      flarePy = hit->GetInitMomentum().y();
      flarePz = hit->GetInitMomentum().z();
      flareDeltaPx = hit->GetDeltaMom().x();
      flareDeltaPy = hit->GetDeltaMom().y();
      flareDeltaPz = hit->GetDeltaMom().z();
      flareEdep = hit->GetEdep();
      flareParticleID = trackIDtoParticleID.at(flareTrackID);

      double hit_position_xyz[3] = {flareX, flareY, flareZ};
      double vtx_xyz[3] = {primaries[0].Vx(), primaries[0].Vy(), primaries[0].Vz()};

      // which primary ancestor does this hit belong to?
      G4int whichPrim = flareAncestorID;
      G4int whichIndex = whichPrim - 1; // need to start from zero

      // pseudo reco: track/shower length and width
      primaryTrackLength[whichIndex] += hit->GetStepLength();
      double ShowerP = primaries[whichIndex].P();
      double dsquare_hit_vtx = TMath::Power((flareX - primaries[whichIndex].Vx()), 2) +
                               TMath::Power((flareY - primaries[whichIndex].Vy()), 2) +
                               TMath::Power((flareZ - primaries[whichIndex].Vz()), 2);
      double product_hit_p = (flareX - primaries[whichIndex].Vx()) * primaries[whichIndex].Px() +
                             (flareY - primaries[whichIndex].Vy()) * primaries[whichIndex].Py() +
                             (flareZ - primaries[whichIndex].Vz()) * primaries[whichIndex].Pz();
      double len_hit = TMath::Abs(product_hit_p) / ShowerP;
      double width_hit = TMath::Sqrt((dsquare_hit_vtx - product_hit_p * product_hit_p / ShowerP / ShowerP));
      ShowerLength[whichIndex] = (ShowerLength[whichIndex] > len_hit) ? ShowerLength[whichIndex] : len_hit;
      double weighted_width_hit = width_hit * flareEdep;
      if (!std::isnan(weighted_width_hit))
        ShowerWidth[whichIndex] += weighted_width_hit;

      if (sdName == "FLArEBoxSD/lar_box")
      {
        nLArHits++;
        fFLArEHits->Fill();

        if (fAddDiffusion == "toy")
          pm3D->FillEntryWithToyElectronTransportation(hit_position_xyz, vtx_xyz, flareEdep, whichIndex);
        else if (fAddDiffusion == "single")
          pm3D->FillEntryWithToySingleElectronTransportation(hit_position_xyz, vtx_xyz, flareEdep, whichIndex);
        else
          pm3D->FillEntry(hit_position_xyz, vtx_xyz, flareEdep, whichIndex);

        // accumulate Edep in LAr
        edepInLAr += flareEdep;
        ProngEInLAr[whichIndex] += flareEdep;
        primaryTrackLengthInTPC[whichIndex] += hit->GetStepLength();
        ShowerLengthInLAr[whichIndex] = (ShowerLengthInLAr[whichIndex] > len_hit) ? ShowerLengthInLAr[whichIndex] : len_hit;
        if (!std::isnan(weighted_width_hit))
          ShowerWidthInLAr[whichIndex] += weighted_width_hit;

        // accumulate dE/dx in LAr
        float_t longitudinal_distance_to_vtx = ((flareX - vtx_xyz[0]) * primaries[0].Px() +
                                                (flareY - vtx_xyz[1]) * primaries[0].Py() +
                                                (flareZ - vtx_xyz[2]) * primaries[0].Pz()) /
                                               primaries[0].P();
        if (longitudinal_distance_to_vtx >= 0 && longitudinal_distance_to_vtx < 3000) // within 3000 mm
          TrueTotalDedxLongitudinal[Int_t(longitudinal_distance_to_vtx)] += hit->GetEdep();
      }
      else if (sdName == "FLArEHadCalXSD/lar_box" ||
               sdName == "FLArEMuonFinderXSD/lar_box" ||
               sdName == "FLArEBabyMINDHorBarSD/lar_box")
      {
        // specify this is an ZX hit
        nHCALHits++;
        flareIsZX = true;
        fFLArEHCALHits->Fill();

        // accumulate Edep in HCAL
        edepInHCALX += flareEdep;
        ProngEInHadCal[whichIndex] += flareEdep;
      }
      else if (sdName == "FLArEHadCalYSD/lar_box" ||
               sdName == "FLArEMuonFinderYSD/lar_box" ||
               sdName == "FLArEBabyMINDVerBarSD/lar_box" ||
               sdName == "FLArEHadAbsorbSD/lar_box" ||
               sdName == "FLArEMuonFinderAbsorbSD/lar_box")
      {
        nHCALHits++;
        flareIsZX = false;
        fFLArEHCALHits->Fill();

        // accumulate Edep in HCAL
        edepInHCALY += flareEdep;
        ProngEInHadCal[whichIndex] += flareEdep;
      }
    }
  }

  if (fSave2DEvd)
    pm3D->Write2DPMToFile(fFile, fFLArEDir);

  pm3D->Process3DPM(fH5file, initPDG, fslPDG, intType, scatteringType, initE, fSave3DEvd);

  sparseFractionMem = pm3D->GetSparseFractionMem();
  sparseFractionBins = pm3D->GetSparseFractionBins();

  if (fSavePseudoReco)
    FillFLArEPseudoReco();

  G4cout << "Total FLArE recorded hits: " << nLArHits << G4endl;
  G4cout << "Total FLArE HCAL recorded hits: " << nHCALHits << G4endl;

  delete pm3D;
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::FillFASER2Output()
{
  G4cout << "==== Filling FASER2 output trees ====" << G4endl;

  // loop over the detected FASER2 sensitive volumes
  int nHits = 0;
  for (const int sdId : fFaser2SDs)
  {
    //  Get and cast hit collection with LArBoxHits
    std::string sdName = fSDNamelist.at(sdId);
    auto hitCollection = dynamic_cast<FASER2TrackerHitsCollection *>(fHCofEvent->GetHC(sdId));
    if (!hitCollection)
    {
      G4cout << "No hits recorded by " << sdName << G4endl;
      continue;
    }

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
      ActsHitsParticleID = trackIDtoParticleID.at(hitID);

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

    } // end of loop over hits
  }
  G4cout << "Total FASER2 recorded hits: " << nHits << G4endl;
}

float_t AnalysisManager::GetTotalEnergy(float_t px, float_t py, float_t pz, float_t m)
{
  return TMath::Sqrt(px * px + py * py + pz * pz + m * m);
}

void AnalysisManager::FillFLArEPseudoReco()
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

  nprimaries = primaries.size();
  primaryPDG.resize(nprimaries, 0.);
  ProngAngleToBeamDir.resize(nprimaries, 0.);
  ProngAvgdEdx.resize(nprimaries, 0.);
  ProngAvgdEdxInLAr.resize(nprimaries, 0.);
  dir_pol_x.resize(nprimaries, 0.);
  dir_pol_y.resize(nprimaries, 0.);
  dir_pol_z.resize(nprimaries, 0.);
  dir_coc_x.resize(nprimaries, 0.);
  dir_coc_y.resize(nprimaries, 0.);
  dir_coc_z.resize(nprimaries, 0.);

  for (auto iPrim : primaryIDs)
  {
    // trackIDs go from 1 to N, you need index: 0 to N-1
    G4int iiPrim = iPrim - 1;

    primaryPDG[iiPrim] = primaries[iiPrim].PDG();

    float_t totProngE = ProngEInLAr[iiPrim] + ProngEInHadCal[iiPrim];
    if (totProngE > 0)
    {
      ShowerWidth[iiPrim] = ShowerWidth[iiPrim] / totProngE;
    }
    if (ProngEInLAr[iiPrim] > 0)
    {
      ShowerWidthInLAr[iiPrim] = ShowerWidthInLAr[iiPrim] / ProngEInLAr[iiPrim];
    }

    double ShowerP = primaries[iiPrim].P();
    double costheta = primaries[iiPrim].Pz() / ShowerP;
    ProngAngleToBeamDir[iiPrim] = TMath::ACos(costheta);

    ProngAvgdEdx[iiPrim] = (ProngEInLAr[iiPrim] +
                            ProngEInHadCal[iiPrim]) /
                           ShowerLength[iiPrim];
    ProngAvgdEdxInLAr[iiPrim] = ProngEInLAr[iiPrim] / ShowerLengthInLAr[iiPrim];

    std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(3);
    std::cout << std::setw(10) << primaries[iiPrim].PDG();
    std::cout << std::setw(12) << ProngAngleToBeamDir[iiPrim];
    std::cout << std::setw(13) << primaryTrackLength[iiPrim];
    std::cout << std::setw(13) << ShowerLength[iiPrim];
    std::cout << std::setw(18) << ShowerWidthInLAr[iiPrim];
    std::cout << std::setw(12) << ProngEInLAr[iiPrim];
    std::cout << std::setw(12) << ProngEInHadCal[iiPrim];
    std::cout << std::setw(12) << ProngAvgdEdxInLAr[iiPrim];
    std::cout << std::setw(10) << primaries[iiPrim].ProngType();
    std::cout << std::setw(12) << primaries[iiPrim].Pz() << std::endl;
  }

  slid::ShowerLID *shwlid = new slid::ShowerLID(pm3D->Get3DPixelMap(),
                                                vtxX, vtxY, vtxZ, 0., 0., 1.);
  Double_t *ptr_dedx = shwlid->GetTotalDedxLongitudinal();
  std::copy(ptr_dedx, ptr_dedx + 3000, TotalDedxLongitudinal.begin());

  for (int iPrim = 0; iPrim < nprimaries; ++iPrim)
  {
    directionfitter::LinearFit *linFit = new directionfitter::LinearFit(
        pm3D->Get2DPixelMapZX(iPrim + 1),
        pm3D->Get2DPixelMapZY(iPrim + 1),
        pm3D->Get2DVtxPixelMapZX(iPrim + 1),
        pm3D->Get2DVtxPixelMapZY(iPrim + 1),
        vtxX, vtxY, vtxZ,
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

  fFLArEPseudoReco->Fill();
}
