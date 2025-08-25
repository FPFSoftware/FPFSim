#ifndef GENIEGenerator_HH
#define GENIEGenerator_HH

#include "generators/GeneratorBase.hh"
#include "G4ParticleDefinition.hh"

#include "TFile.h"
#include "TTree.h"
#include "globals.hh"

class G4Event;

class GENIEGenerator : public GeneratorBase
{
  public:
    GENIEGenerator();
    ~GENIEGenerator();

    // override methods from common base class
    void LoadData() override;
    void GeneratePrimaries(G4Event *anEvent) override;

    // setter methods for messenger
    void SetGSTFilename(G4String val) { fGSTFilename = val; }
    void SetEvtStartIdx(G4int val) { fEvtStartIdx = val; }
    void SetRandomVertex(G4bool val) { fRandomVtx = val; }

  private:
    G4String fGSTFilename;
    G4int fEventCounter;
    G4int fEvtStartIdx;
    G4bool fRandomVtx;
    TFile *fGSTFile;
    TTree *fGSTTree;

    // gst tree branches
    // define the branches we are interested in from the GST tree
    const G4int fkNPmax = 250;
    G4int fNEntries;
    G4double m_wght;
    G4double m_Ev, m_pxv, m_pyv, m_pzv;
    G4double m_El, m_pxl, m_pyl, m_pzl;
    G4int m_nf;
    G4int m_pdgf[250];
    G4double m_Ef[250], m_pxf[250], m_pyf[250], m_pzf[250];
    G4bool m_qel, m_mec, m_res, m_dis, m_coh, m_dfr, m_imd;
    G4bool m_imdanh, m_nuel, m_em, m_cc, m_nc, m_charm;
    G4bool m_singlek, m_amnugamma;
    G4int m_neuPDG, m_fslPDG;
    G4int m_tgt, m_Z, m_A, m_hitnuc;
    G4double m_x, m_y, m_Q2, m_W;

    // specific internal functions
    G4bool FindParticleDefinition(G4int pdg, G4ParticleDefinition* &particleDefinition) const;
    G4int DecodeInteractionType() const;
    G4int DecodeScatteringType() const;
    G4String EncodeProcessName() const;
};

#endif
