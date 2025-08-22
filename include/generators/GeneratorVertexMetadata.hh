#ifndef GENERATOR_VERTEX_METDATA_HH
#define GENERATOR_VERTEX_METDATA_HH

#include <string>
#include "TLorentzVector.h"

// Metadata information on each generated vertex
// Store input interaction/decay process information
struct GeneratorVertexMetadata 
{
  std::string generatorType;  ///< which generator class
  std::string processName;    ///< which process
  double weight;              ///< process weight

  int pdg;           ///< initiator pdg
  TLorentzVector x4; ///< initiator 4-position (vertex)
  TLorentzVector p4; ///< initiator 4-momentum
  double mass;     ///< initiator mass
  double charge;   ///< initiator charge

  int tgt_pdg;     ///< target pdg
  int tgt_A;       ///< nuclear target A
  int tgt_Z;       ///< nuclear target Z
  
  double xs = 0.0;   ///< cross-section
  double Q2 = -1.0;  ///< momentum transfer
  double xBj =  0.0; ///< Bjorken x
  double y = -1.0;   ///< inelasticity
  double W = -1.0;   ///< hadronic invariant mass
};

#endif