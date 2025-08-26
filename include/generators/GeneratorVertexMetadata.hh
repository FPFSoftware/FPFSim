#ifndef GENERATOR_VERTEX_METDATA_HH
#define GENERATOR_VERTEX_METDATA_HH

#include <string>
#include "G4LorentzVector.hh"

// Metadata information on each generated vertex
// Store input interaction/decay process information
struct GeneratorVertexMetadata 
{
  std::string generatorType;  ///< which generator class
  std::string processName;    ///< which process
  double weight = 1.0;        ///< process weight

  int pdg = -1;      ///< initiator pdg
  G4LorentzVector x4 = G4LorentzVector(); ///< initiator 4-position (vertex)
  G4LorentzVector p4 = G4LorentzVector(); ///< initiator 4-momentum
  double mass = -1;  ///< initiator mass
  double charge = 0; ///< initiator charge

  int intType = -1;           ///< interaction type 
  int scatteringType = -1;    ///< scattering type 
  int fsl_pdg = -1;           ///< final state lepton pdg

  int tgt_pdg = -1;     ///< target pdg
  int tgt_A = -1;       ///< nuclear target A
  int tgt_Z = -1;       ///< nuclear target Z
  int hitnuc_pdg = -1;  ///< hit nucleon pdg
  
  double xs = -1.0;  ///< cross-section
  double Q2 = -1.0;  ///< momentum transfer
  double xBj = 0.0;  ///< Bjorken x
  double y = -1.0;   ///< inelasticity
  double W = -1.0;   ///< hadronic invariant mass

};

#endif