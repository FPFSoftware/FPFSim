#include "generators/GeneratorBase.hh"
#include "generators/GENIEGenerator.hh"
#include "generators/GENIEGeneratorMessenger.hh"
#include "generators/GeneratorVertexMetadata.hh"

#include "geometry/GeometricalParameters.hh"

#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Exception.hh"
#include "G4LorentzVector.hh"
#include "Randomize.hh"

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"

GENIEGenerator::GENIEGenerator()
{
  fGeneratorName = "genie";
  fEvtStartIdx = 0;
  fMessenger = new GENIEGeneratorMessenger(this);

  fGSTFile = nullptr;
  fGSTTree = nullptr;
  fRandomVtx = true;
  fEventCounter = 0;
}

GENIEGenerator::~GENIEGenerator()
{
  delete fGSTTree;
  if(fGSTFile) fGSTFile->Close();
  delete fMessenger;
}

void GENIEGenerator::LoadData()
{

  fGSTFile = new TFile(fGSTFilename, "read");
  if (!fGSTFile->IsOpen()) {
    G4String err = "Cannot open GST file : " + fGSTFilename;
    G4Exception("GENIEGenerator", "FileError", FatalErrorInArgument, err.c_str());
  }

  fGSTTree = (TTree*)fGSTFile->Get("gst");
  if (!fGSTTree) {
    G4String err = "No GST event tree in input file : " + fGSTFilename;
    G4Exception("GENIEGenerator", "FileError", FatalErrorInArgument, err.c_str());
  }

  fNEntries = fGSTTree->GetEntries();
  G4cout << "Input GST tree has " << fNEntries << ((fNEntries==1)? " entry." : " entries.") << G4endl;

  fGSTTree->SetBranchAddress("qel",&m_qel); // is QEL?   
  fGSTTree->SetBranchAddress("mec",&m_mec); // is MEC?
  fGSTTree->SetBranchAddress("res",&m_res); // is RES?
  fGSTTree->SetBranchAddress("dis",&m_dis); // is DIS?
  fGSTTree->SetBranchAddress("coh",&m_coh); // is Coherent?
  fGSTTree->SetBranchAddress("dfr",&m_dfr); // id Diffractive?
  fGSTTree->SetBranchAddress("imd",&m_imd); // is IMD?
  fGSTTree->SetBranchAddress("imdanh",&m_imdanh); // is IMD annihilation?
  fGSTTree->SetBranchAddress("singlek",&m_singlek); // is single Kaon?
  fGSTTree->SetBranchAddress("nuel",&m_nuel);  // is ve elastic?
  fGSTTree->SetBranchAddress("em",&m_em); // is EM process?
  fGSTTree->SetBranchAddress("cc",&m_cc); // is Weak CC?
  fGSTTree->SetBranchAddress("nc",&m_nc); // is Weak NC?
  fGSTTree->SetBranchAddress("charm",&m_charm); // produces charm?
  fGSTTree->SetBranchAddress("amnugamma",&m_amnugamma); // is anomaly mediated nu gamma?

  fGSTTree->SetBranchAddress("neu",&m_neuPDG); //neutrino PDG
  fGSTTree->SetBranchAddress("Ev",&m_Ev); // neutrino energy (GeV)
  fGSTTree->SetBranchAddress("pxv",&m_pxv); // neutrino px (GeV)
  fGSTTree->SetBranchAddress("pyv",&m_pyv); // neutrino py (Gev)
  fGSTTree->SetBranchAddress("pzv",&m_pzv); // neutrino pz (GeV)
 
  fGSTTree->SetBranchAddress("fspl",&m_fslPDG); // primary letpton PDG
  fGSTTree->SetBranchAddress("El",&m_El); // primary lepton energy (GeV)
  fGSTTree->SetBranchAddress("pxl",&m_pxl); // primary lepton px (GeV)
  fGSTTree->SetBranchAddress("pyl",&m_pyl); // primary lepton py (Gev)
  fGSTTree->SetBranchAddress("pzl",&m_pzl); // primary lepton pz (GeV)
  
  fGSTTree->SetBranchAddress("nf",&m_nf); // number of final state hadrons
  fGSTTree->SetBranchAddress("pdgf",&m_pdgf); // hadrons PDG
  fGSTTree->SetBranchAddress("Ef",&m_Ef); // hadrons energy (GeV)
  fGSTTree->SetBranchAddress("pxf",&m_pxf); // hadrons px (GeV)
  fGSTTree->SetBranchAddress("pyf",&m_pyf); // hadrons py (Gev)
  fGSTTree->SetBranchAddress("pzf",&m_pzf); // hadrons pz (GeV)
  
  fGSTTree->SetBranchAddress("W",&m_W); // invariant hadronic mass (GeV)
  fGSTTree->SetBranchAddress("Q2",&m_Q2); // momentum transfer (GeV^2)
  fGSTTree->SetBranchAddress("x",&m_x); // Bjorken x
  fGSTTree->SetBranchAddress("y",&m_y); // inelasticity

  fGSTTree->SetBranchAddress("wght",&m_wght); // event weigth
  
  fGSTTree->SetBranchAddress("tgt",&m_tgt); // nuclear target pdg
  fGSTTree->SetBranchAddress("Z",&m_Z); // nuclear target Z
  fGSTTree->SetBranchAddress("A",&m_A); // nuclear target A
  fGSTTree->SetBranchAddress("hitnuc",&m_hitnuc); // hit nucleon pfg

}

G4bool GENIEGenerator::FindParticleDefinition(G4int const pdg, G4ParticleDefinition* &particleDefinition) const
{
  // unknown pgd codes in GENIE --> skip it!
  // ref: https://internal.dunescience.org/doxygen/ConvertMCTruthToG4_8cxx_source.html
  // This has been a known issue with GENIE
  const int genieLo = 2000000001;
  const int genieHi = 2000000202;
  if ( pdg >= genieLo && pdg <= genieHi) {
    G4cout<< "This unknown PDG code [" << pdg << "] was present in the GENIE input, "
             << "but will not be processed by Geant4."
             << G4endl;
    return false; // return bad
  }
  
  if ( pdg == 0) {
    particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");
  } else {
    particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle(pdg);
  }
  
  if ( pdg > 1000000000) { // If the particle is a nucleus
    // If the particle table doesn't have a definition yet, ask the ion
    // table for one. This will create a new ion definition as needed.
    if (!particleDefinition) {
      int Z = (pdg % 10000000) / 10000; // atomic number
      int A = (pdg % 10000) / 10;       // mass number
      particleDefinition = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(Z, A, 0.);
    }
  }

  return true; //return good
}

void GENIEGenerator::GeneratePrimaries(G4Event* anEvent)
{

  // complete line from PrimaryGeneratorAction...
  G4cout << ") : GENIE Generator ===oooOOOooo===" << G4endl;
  
  G4cout << "oooOOOooo Event # " << fEventCounter << " oooOOOooo" << G4endl;
  G4cout << "GeneratePrimaries from file " << fGSTFilename << ", evtID starts from "<< fEvtStartIdx << ", now at " << fEvtStartIdx+fEventCounter << G4endl;

  G4int currentIdx = fEvtStartIdx+fEventCounter;
  anEvent->SetEventID(currentIdx);

  if ( currentIdx >= fNEntries ) {
    G4cerr << "** event index beyond range !! **" << G4endl;
  }
  
  // fetch a single entry from GENIE input file
  fGSTTree->GetEntry(currentIdx); 

  // compute/repackage what is not directly available from the tree
  // position is randomly extracted in the detector fiducial volume
  // or set to the center according to config parameter 
  G4LorentzVector neuP4(m_pxv*GeV,m_pyv*GeV,m_pzv*GeV,m_Ev*GeV);
  G4LorentzVector fslP4(m_pxl*GeV,m_pyl*GeV,m_pzl*GeV,m_El*GeV);
  G4LorentzVector neuX4;

  G4Random::setTheSeed(currentIdx+1);
  if(fRandomVtx){
    neuX4.setX(GeometricalParameters::Get()->GetFLArEPosition().x() +
              (G4UniformRand()-0.5) * GeometricalParameters::Get()->GetFLArEFidVolSize().x());
    neuX4.setY(GeometricalParameters::Get()->GetFLArEPosition().y() +
              (G4UniformRand()-0.5) * GeometricalParameters::Get()->GetFLArEFidVolSize().y());
    neuX4.setZ(GeometricalParameters::Get()->GetFLArEPosition().z() +
              (G4UniformRand()-0.5) * GeometricalParameters::Get()->GetFLArEFidVolSize().z());
    neuX4.setT(0.);
  } else {
    neuX4.setX(0.*m);
    neuX4.setY(0.*m);
    neuX4.setZ(GeometricalParameters::Get()->GetFLArEPosition().z() -
                 GeometricalParameters::Get()->GetFLArEFidVolSize().z()/2);
    neuX4.setT(0.);
  }

  // create primary vertex (neutrino)
  G4PrimaryVertex* vtx = new G4PrimaryVertex(neuX4.x(), neuX4.y(), neuX4.z(), neuX4.t()); 
  // now add all the final state particles into the vertex
  // - final state lepton (if NC, it's the neutrino!)
  G4ParticleDefinition* particleDefinition;
  if ( FindParticleDefinition(m_fslPDG, particleDefinition) ){

    G4PrimaryParticle* plepton = new G4PrimaryParticle(particleDefinition, fslP4.x(), fslP4.y(), fslP4.z(), fslP4.e());
    /* G4cout << "Lepton PDG " << m_fslPDG << " mass " << particleDefinition->GetPDGMass()*MeV << G4endl;
    G4cout << "p4  " << m_fslP4.X() << " " << m_fslP4.Y() << " " <<  m_fslP4.Z() << " " <<  m_fslP4.E() << G4endl;
    G4cout << "kinE " << ( m_fslP4.E()*GeV - particleDefinition->GetPDGMass()*MeV)  << G4endl; */
    vtx->SetPrimary(plepton);

  }

  // - all final state hadrons
  for (int ipar=0; ipar<m_nf; ++ipar) {

    G4LorentzVector p( m_pxf[ipar]*GeV, m_pyf[ipar]*GeV, m_pzf[ipar]*GeV, m_Ef[ipar]*GeV );
    if ( !FindParticleDefinition( m_pdgf[ipar], particleDefinition) ) continue; //skip bad pdgs
    G4PrimaryParticle* prim = new G4PrimaryParticle(particleDefinition, p.x(), p.y(), p.z(), p.t()); 
    /* G4cout << "Particle PDG " << m_pdgf[ipar] << " mass " << particleDefinition->GetPDGMass()*MeV << G4endl;
    G4cout << "p4  " << p.X() << " " << p.Y() << " " << p.Z() << " " << p.E() << G4endl;
    G4cout << "kinE " << (p.E()*GeV - particleDefinition->GetPDGMass()*MeV)  << G4endl; */
    vtx->SetPrimary(prim);

  }

  // package and ship metadata
  GeneratorVertexMetadata metadata;
  metadata.generatorType = fGeneratorName;
  metadata.processName = EncodeProcessName();
  metadata.weight = m_wght;
  metadata.pdg = m_neuPDG; 
  metadata.x4 = neuX4;
  metadata.p4 = neuP4;
  metadata.mass = 0.; // neutrinos always massless
  metadata.charge = 0.; // neutrinos always no charge
  metadata.intType = DecodeInteractionType();
  metadata.scatteringType = DecodeScatteringType();
  metadata.fsl_pdg = m_fslPDG;
  metadata.tgt_pdg = m_tgt;
  metadata.tgt_A = m_A;
  metadata.tgt_Z = m_Z;
  metadata.hitnuc_pdg = m_hitnuc;
  metadata.Q2 = m_Q2;
  metadata.xBj = m_x;
  metadata.y = m_y;
  metadata.W = m_W;
  fVertexMetadata.push_back(metadata);

  anEvent->AddPrimaryVertex(vtx);
  fEventCounter++;

}

G4int GENIEGenerator::DecodeInteractionType() const
{
  // using ref: https://internal.dunescience.org/doxygen/Generator_2src_2Framework_2Interaction_2InteractionType_8h_source.html
  // returning the value of the corresponding enum member from the flags that are in output.
  // there are no flags in the gst tree for some of the interaction types...
 
  enum InteractionType {
    kIntNull   = 0,
    kIntEM,          //
    kIntWeakCC,      //
    kIntWeakNC,      //
    kIntWeakMix,     // CC + NC + interference
    kIntDarkMatter,  //
    kIntNDecay,      //
    kIntNOsc,        //
    kIntNHL,         //
    kIntDarkNC       // 
  };
	
  if (m_cc) return kIntWeakCC;
  else if (m_nc) return kIntWeakNC;
  else if (m_em) return kIntEM;

  return kIntNull;

}

G4int GENIEGenerator::DecodeScatteringType() const
{
  // using ref: https://internal.dunescience.org/doxygen/ScatteringType_8h_source.html
  // returning the value of the corresponding enum member from the flags that are in output.
  // there are no flags in the gst tree for some of the scattering types...
  
  enum ScatteringType {
    kScUnknown = -100,
    kScNull = 0,
    kScQuasiElastic,
    kScSingleKaon,
    kScDeepInelastic,
    kScResonant,
    kScCoherentProduction,
    kScDiffractive,
    kScNuElectronElastic,
    kScInverseMuDecay,
    kScAMNuGamma,
    kScMEC,
    kScCoherentElastic,
    kScInverseBetaDecay,
    kScGlashowResonance,
    kScIMDAnnihilation,
    kScPhotonCoherent,
    kScPhotonResonance,
    kScSinglePion,
    kScDarkMatterElastic = 101,
    kScDarkMatterDeepInelastic,
    kScDarkMatterElectron,
    kScNorm
  };

  if(m_qel) return kScQuasiElastic;
  else if(m_res) return kScResonant;
  else if(m_dis) return kScDeepInelastic;
  else if(m_coh) return kScCoherentProduction;
  else if(m_dfr) return kScDiffractive;
  else if(m_imd) return kScInverseMuDecay;
  else if(m_imdanh) return kScIMDAnnihilation;
  else if(m_mec) return kScMEC;
  else if(m_nuel) return kScNuElectronElastic;
  else if(m_singlek) return kScSingleKaon;
  else if(m_amnugamma) return kScAMNuGamma;

  return kScUnknown;
}

G4String GENIEGenerator::EncodeProcessName() const
{
  G4String process = "";
  if (m_cc) process += "WeakCC";
  else if (m_nc) process += "WeakNC";
  else if (m_em) process += "EM";

  if(m_qel) process += " QE";
  else if(m_res) process += " RES";
  else if(m_dis) process += " DIS";
  else if(m_coh) process += " COH";
  else if(m_imd) process += " IMD";
  else if(m_mec) process += " MEC";
  else if(m_nuel) process += " nuELASTIC";

  return process;
}