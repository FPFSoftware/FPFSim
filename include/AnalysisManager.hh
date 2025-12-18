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
    void parKinECut(G4double val) { fParKinECut = val; }
    void saveActs(G4bool val) { fSaveActs = val; }
    void savePseudoReco(G4bool val) { fSavePseudoReco = val; }
    void addDiffusion(G4String val) { fAddDiffusion = val; } 
    void save3DEvd(G4bool val) { fSave3DEvd = val; } 
    void save2DEvd(G4bool val) { fSave2DEvd = val; }
    void enableFLArE(G4bool val) { fEnableFLArE = val;}

    // build TID to primary ancestor association
    // filled progressively from StackingAction
    void SetTrackPrimaryAncestor(G4int trackID, G4int ancestorID) { trackToPrimaryAncestor[trackID] = ancestorID; }
    G4int GetTrackPrimaryAncestor(G4int trackID) { return trackToPrimaryAncestor.at(trackID); }

    // return whether saving full tracks in trajectories
    G4bool GetSaveTrack() { return fSaveTrack; }

  private:

    //------------------------------------------------
    // Book ROOT output TTrees
    // common + detector specific
    void bookEvtTree();  
    void bookParTree();  
    void bookFLArETrees();      
    void bookFASER2Trees();

    void FillEventTree(const G4Event* event);
    void FillParticlesTree(const G4Event* event);
    
    void FillFLArEOutput();
    void FillFLArEPseudoReco();

    void FillFASER2Output();

    float_t GetTotalEnergy(float_t px, float_t py, float_t pz, float_t m);

    static AnalysisManager* fInstance;
    AnalysisManagerMessenger* fMessenger;

    G4bool fSaveTrack;
    G4double fParKinECut;
    G4bool fSave3DEvd;
    G4bool fSave2DEvd;
    G4bool fSavePseudoReco;
    TString fAddDiffusion;
    G4bool fSaveActs;
    G4bool fEnableFLArE;

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
    TTree*   fPar;
    TTree*   fPrim;

    TDirectory* fFLArEDir;
    TTree*   fFLArEHits;
	  TTree*	 fFLArEHCALHits;
    TTree*   fFLArEPseudoReco; 

    TDirectory* fFASER2Dir;
    TTree*   fActsHitsTree;

    // track to primary ancestor (track id to track id)
    std::map<G4int, G4int> trackToPrimaryAncestor;

    // map trackID to its barcode
    std::map<G4int, ULong64_t> trackIDtoParticleID;

    //---------------------------------------------------
    // OUTPUT VARIABLES FOR COMMON TREES
    //---------------------------------------------------
    // Output variables for EVENT tree

    G4int evtID;
    G4int vertexID;
    double weight;
    std::string genType;
    std::string processName; 
    double vtxX, vtxY, vtxZ, vtxT;
    int initPDG;           
    double initPx, initPy, initPz, initE;
    double initM;     
    double initQ;  
    int fslPDG;           
    int tgtPDG;     
    int tgtA;      
    int tgtZ;      
    int hitnucPDG; 
    int intType;           
    int scatteringType;    
    double xs;
    double Q2;  
    double xBj; 
    double y;   
    double W;  
    int nPrimaries;
    std::vector<int> primTID;
    std::vector<int> primPDG;
    std::vector<float> primPx;
    std::vector<float> primPy;
    std::vector<float> primPz;
    std::vector<float> primE;

    //---------------------------------------------------
    // Output variables for PARTICLES/TRAJECTORIES tree
    
    std::uint64_t particle_id; //barcode
    int particle_TID;
    int particle_PID;
    int particle_PDG;
    std::string particle_process;
    float particle_vx;
    float particle_vy;
    float particle_vz;
    float particle_vt;
    float particle_px;
    float particle_py;
    float particle_pz;
    float particle_m;
    float particle_q;
    float particle_eta;
    float particle_phi;
    float particle_pt;
    float particle_p;
    float particle_ke;
    int traj_Npoints;
    std::vector<float> traj_pointX;
    std::vector<float> traj_pointY;
    std::vector<float> traj_pointZ;

    //---------------------------------------------------
    // OUTPUT VARIABLES FOR FLArE TREES
    // TODO: merge hit variables? no need to use different names?

    PixelMap3D* pm3D;

    UInt_t flareTrackID;
    UInt_t flareParticleID;
    UInt_t flareParentID;
    UInt_t flarePDG;
    UInt_t flareCopyNum;
    float_t flareT;
    float_t flareX;
    float_t flareY;
    float_t flareZ;
    float_t flarePx;
    float_t flarePy;
    float_t flarePz;
    float_t flareDeltaPx;
    float_t flareDeltaPy;
    float_t flareDeltaPz;
    float_t flareEdep;
    bool flareIsZX;

    // flare pseudo-reco
    float_t edepInLAr;
    float_t edepInHCALX;
    float_t edepInHCALY;
    float_t sparseFractionMem;
    float_t sparseFractionBins;
    std::vector<float_t> TotalDedxLongitudinal;
    std::vector<float_t> TrueTotalDedxLongitudinal;
    G4int nprimaries;
    std::vector<float_t> primaryPDG;
    std::vector<float_t> primaryTrackLength;      
    std::vector<float_t> primaryTrackLengthInTPC; 
    std::vector<float_t> ProngEInLAr;
    std::vector<float_t> ProngEInHadCal;
    std::vector<float_t> ProngAngleToBeamDir;
    std::vector<float_t> ShowerLength;
    std::vector<float_t> ShowerLengthInLAr;
    std::vector<float_t> ShowerWidth;
    std::vector<float_t> ShowerWidthInLAr;
    std::vector<float_t> ProngAvgdEdx;
    std::vector<float_t> ProngAvgdEdxInLAr;
    std::vector<float_t> dir_pol_x;
    std::vector<float_t> dir_pol_y;
    std::vector<float_t> dir_pol_z;
    std::vector<float_t> dir_coc_x;
    std::vector<float_t> dir_coc_y;
    std::vector<float_t> dir_coc_z;

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

};

#endif
