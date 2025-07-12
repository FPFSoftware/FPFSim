#ifndef ANALYSISMANAGER_HH
#define ANALYSISMANAGER_HH

#include <set>
#include <vector>
#include <string>

#include "globals.hh"
#include "G4Event.hh"
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"

#include "AnalysisManagerMessenger.hh"
#include "PixelMap3D.hh"
#include "FPFParticle.hh"
#include "FPFNeutrino.hh"

#include "hep_hpc/hdf5/File.hpp"

class AnalysisManager {
  public:

    AnalysisManager();
    ~AnalysisManager();
    static AnalysisManager* GetInstance();

    //------------------------------------------------
    // Functions executed at specific times
    void BeginOfRun(); 
    void EndOfRun();
    void BeginOfEvent();
    void EndOfEvent(const G4Event* event);

    //------------------------------------------------
    // functions for controlling from the configuration file
    void setFileName(std::string val) { fFilename = val; }
    void saveTrack(G4bool val) { fSaveTrack = val; }
    void saveActs(G4bool val) { fSaveActs = val; }
    void addDiffusion(G4String val) { fAddDiffusion = val; } // TODO: needed?
    void save3DEvd(G4bool val) { fSave3DEvd = val; } // TODO: needed?
    void save2DEvd(G4bool val) { fSave2DEvd = val; } // TODO: needed?

    // TODO: needed???
    void SetTrackPTPair(G4int PID, G4int TID) { allTracksPTPair.insert(std::make_pair(PID, TID)); }
    void AddOnePrimaryTrack() { nTestNPrimaryTrack++; }

  private:

    //------------------------------------------------
    // Book ROOT output TTrees
    // common + detector specific
    void bookEvtTree();  
    void bookTrkTree();  
    void bookPrimTree(); 
    void bookFLArETrees();      
    void bookFASER2Trees();

    void FillPrimariesTree(const G4Event* event);
    void FillTrajectoriesTree(const G4Event* event);
    
    void FillFLArEOutput();
    void FillFASER2Output();

    double GetTotalEnergy(double px, double py, double pz, double m);

    // TODO: needed??
    void FillPseudoRecoVar();

    static AnalysisManager* fInstance;
    AnalysisManagerMessenger* fMessenger;

    G4bool fSaveTrack;
    G4bool fSave3DEvd;
    G4bool fSave2DEvd;
    G4bool fSaveActs;
    TString fAddDiffusion;

    std::map<int, std::string> fSDNamelist;
    std::vector<int> fFlareSDs;
    std::vector<int> fFaser2SDs;
    
    G4HCofThisEvent* fHCofEvent;
    
    G4int nPrimaryVertex;
    std::vector<FPFParticle> primaries;
    std::vector<int> primaryIDs;

    //------------------------------------------------
    // output files and trees
    std::string fFilename;
    std::string fH5Filename;
    hep_hpc::hdf5::File fH5file;
    TFile*   fFile;
    TTree*   fEvt;
    TTree*   fTrk;
    TTree*   fPrim;

    TDirectory* fFLArEDir;
    TTree*   fFLArEHits;
	  TTree*	 fFLArEHCALHits;
    TDirectory* fFASER2Dir;
    TTree*   fActsHitsTree;
    TTree*   fActsParticlesTree;

    // TODO: needed???  
    std::set<std::pair<int, int> > allTracksPTPair;
    G4int    nTestNPrimaryTrack;
    /* 
    std::vector<std::set<int> > trackClusters;
    std::set<int> tracksFromFSL;
    std::set<int> tracksFromFSLSecondary;
    std::set<int> tracksFromFSPizeroSecondary;
    std::set<int> tracksFromFSLDecayPizeroSecondary;
    int fPrimIdxFSL;
    */

    //---------------------------------------------------
    // OUTPUT VARIABLES FOR COMMON TREES
    // TODO: review carefully 
    // TODO: need to make evt tree less FLARE-centric
    // TODO: remove pseudo-reco or add it as FLARE-only tree
    // TODO: turn arrays in std::vector!

    G4int    evtID;
    /* FPFNeutrino neutrino;
 
    // Truth information from genie
    G4int    nuIdx;             ///<- neutrino index (for genie neutrino interaction)
    G4int    nuPDG;             ///<- neutrino PDG code (for genie neutrino interaction)
    G4double nuE;               ///<- neutrino energy
    G4double nuX;               ///<- neutrino vertex X
    G4double nuY;               ///<- neutrino vertex Y
    G4double nuZ;               ///<- neutrino vertex Z
    G4int    nuIntType;         ///<- interaction type: CC, NC, et.al.
    G4int    nuScatteringType;  ///<- scattering type: QE, DIS, RES, MEC, et. al.
    G4double nuW;               ///<- invariant hadronic mass
    G4int    nuFSLPDG;          ///<- Final state lepton PDG code (for genie neutrino interaction)
    G4double nuFSLPx;           ///<- Final state lepton Px
    G4double nuFSLPy;           ///<- Final state lepton Py
    G4double nuFSLPz;           ///<- Final state lepton Pz
    G4double nuFSLE;            ///<- Final state lepton total energy (GeV)

    G4int    countPrimaryParticle;
    G4int    nPrimaryParticle;  ///<- number of primary particle
                                ///<- (in case of genie neutrino interaction, number of stable particle in the final state)
                                ///<- (in case of the FSL decay, decay products counted as primary particle)
                                ///<- (in case of the final state pizero, decay products counted as primary particle)
    //// Geant4 truth
    G4int    primaryParentPDG[1000];        ///<- parent PDG of primary particles
    G4double primaryTrackLength[1000];      ///<- track length of primary particles
    G4double primaryTrackLengthInTPC[1000]; ///<- track length of primary particles in TPC region

    // pseudo-reco
    G4double ProngEInDetector[1000];
    G4double ProngEInLAr[1000];
    G4double ProngEInHadCal[1000];
    G4double ProngEInMuonFinder[1000];
    G4double ProngEInMuonFinderLayer1X[1000];
    G4double ProngEInMuonFinderLayer1Y[1000];
    G4double ProngEInMuonFinderLayer2X[1000];
    G4double ProngEInMuonFinderLayer2Y[1000];
    G4double ProngAngleToBeamDir[1000];
    G4double ShowerLength[1000];
    G4double ShowerLengthInLAr[1000];
    G4double ShowerWidth[1000];
    G4double ShowerWidthInLAr[1000];
    G4double ProngAvgdEdx[1000];
    G4double ProngAvgdEdxInLAr[1000];
    G4double ProngdEdxAlongTrack[1000][100];
    G4int    ProngdEdxTrackLength[1000][100];
    G4double TotalDedxLongitudinal[3000];
    G4double TrueTotalDedxLongitudinal[3000];
    // reco
    // direction
    G4double dir_pol_x[1000];
    G4double dir_pol_y[1000];
    G4double dir_pol_z[1000];
    G4double dir_coc_x[1000];
    G4double dir_coc_y[1000];
    G4double dir_coc_z[1000];

    G4double edepInLAr;
    G4double edepInHadCalX;
    G4double edepInHadCalY;
    G4double edepInMuonFinderX;
    G4double edepInMuonFinderY;
    G4double edepInHadAborb;
    G4double edepInMuonFinderAbsorb;
    G4double missCountedEnergy;

    G4int    nFromFSLParticles;
    G4int    nFromFSPizeroParticles;
    G4int    nFromFSLDecayPizeroParticles;
    G4int    fromFSLParticlePDG[2000000];

    PixelMap3D* pm3D;
    G4double sparseFractionMem;
    G4double sparseFractionBins;
    */

    //---------------------------------------------------
    // Output variables for TRAJECTORIES tree
    int trackTID;
    int trackPID;
    int trackPDG;
    double trackKinE;
    int trackNPoints;
    std::vector<double> trackPointX;
    std::vector<double> trackPointY;
    std::vector<double> trackPointZ;

    //---------------------------------------------------
    // Output variables for PRIMARIES tree
    UInt_t primVtxID;
    UInt_t primParticleID;
    UInt_t primTrackID;
    UInt_t primPDG; // why unsigned?
    float_t primM;
    float_t primQ;
    float_t primEta;
    float_t primPhi;
    float_t primPt;
    float_t primP;
    float_t primVx;
    float_t primVy;
    float_t primVz;
    float_t primVt;
    float_t primPx;
    float_t primPy;
    float_t primPz;
    float_t primE;
    float_t primKE;

    //---------------------------------------------------
    // OUTPUT VARIABLES FOR FLArE TREES
    // TODO: merge hit variables? no need to use different names?
    // TODO: somehow need to add back here info from evt tree

    UInt_t flareTrackID;
    UInt_t flareParticleID;
    UInt_t flareParentID;
    UInt_t flarePDG;
    UInt_t flareCopyNum;
    double flareT;
    double flareX;
    double flareY;
    double flareZ;
    double flarePx;
    double flarePy;
    double flarePz;
    double flareDeltaPx;
    double flareDeltaPy;
    double flareDeltaPz;
    double flareEdep;
    bool flareIsZX;

    //---------------------------------------------------
    // OUTPUT VARIABLES FOR FASER2 TREES

    // Acts Hit Information - the types are set to match the types expected by Acts::RootSimHitReader
    UInt_t ActsHitsEventID;
    ULong64_t ActsHitsGeometryID;
    ULong64_t ActsHitsParticleID;
    Float_t ActsHitsX;
    Float_t ActsHitsY;
    Float_t ActsHitsZ;
    Float_t ActsHitsT;
    Float_t ActsHitsPx;
    Float_t ActsHitsPy;
    Float_t ActsHitsPz;
    Float_t ActsHitsE;
    Float_t ActsHitsDeltaPx;
    Float_t ActsHitsDeltaPy;
    Float_t ActsHitsDeltaPz;
    Float_t ActsHitsDeltaE;
    Int_t ActsHitsIndex;
    UInt_t ActsHitsVolumeID;
    UInt_t ActsHitsBoundaryID;
    UInt_t ActsHitsLayerID;
    UInt_t ActsHitsApproachID;
    UInt_t ActsHitsSensitiveID;

    // Acts Particle Information - need the truth info on the particles in order to do the truth tracking
    std::vector<std::uint64_t> ActsParticlesParticleId;
    std::vector<std::int32_t> ActsParticlesParticleType;
    std::vector<std::uint32_t> ActsParticlesProcess;
    std::vector<float> ActsParticlesVx;
    std::vector<float> ActsParticlesVy;
    std::vector<float> ActsParticlesVz;
    std::vector<float> ActsParticlesVt;
    std::vector<float> ActsParticlesPx;
    std::vector<float> ActsParticlesPy;
    std::vector<float> ActsParticlesPz;
    std::vector<float> ActsParticlesM;
    std::vector<float> ActsParticlesQ;
    std::vector<float> ActsParticlesEta;
    std::vector<float> ActsParticlesPhi;
    std::vector<float> ActsParticlesPt;
    std::vector<float> ActsParticlesP;
    std::vector<std::uint32_t> ActsParticlesVertexPrimary;
    std::vector<std::uint32_t> ActsParticlesVertexSecondary;
    std::vector<std::uint32_t> ActsParticlesParticle;

    std::vector<std::uint32_t> ActsParticlesGeneration;
    std::vector<std::uint32_t> ActsParticlesSubParticle;
    std::vector<float> ActsParticlesELoss;
    std::vector<float> ActsParticlesPathInX0;
    std::vector<float> ActsParticlesPathInL0;
    std::vector<std::int32_t> ActsParticlesNumberOfHits;
    std::vector<std::uint32_t> ActsParticlesOutcome;

};

#endif
