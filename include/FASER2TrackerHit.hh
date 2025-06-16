#ifndef FASER2TRACKERHIT_HH
#define FASER2TRACKERHIT_HH

#include "G4THitsCollection.hh"
#include "G4VHit.hh"
#include "G4Types.hh"
#include <vector>
#include "G4ThreeVector.hh"
#include  "G4LorentzVector.hh"

class FASER2TrackerHit : public G4VHit {
public:
  FASER2TrackerHit(){};
  ~FASER2TrackerHit(){};

  ///// Get hit ID calculated as 1000 * sensorID + cellID
  //G4int ID() { return 1000 * fCopyNumSensor + fCopyNumCell; }
  inline void SetPosition(G4double x, G4double y, G4double z) {
    fPosX = x;
    fPosY = y;
    fPosZ = z;
  }

  inline void SetPDGID(G4int pdgid){
    fPDGID = pdgid;
  }

  inline void SetEnergy(G4double E){
    fEnergy = E;
  }

  inline void SetCharge(G4double charge){
    fCharge = charge;
  }

  inline void SetPx(G4double px){
    fPx = px;
  }
  
  inline void SetPy(G4double py){
    fPy = py;
  }

  inline void SetPz(G4double pz){
    fPz = pz;
  }
  
  inline void SetT(G4double t){
    fT = t;
  }

  inline void SetDeltaPx(G4double deltaPx){
    fDeltaPx = deltaPx;
  }
  inline void SetDeltaPy(G4double deltaPy){
    fDeltaPy = deltaPy;
  }
  inline void SetDeltaPz(G4double deltaPz){
    fDeltaPz = deltaPz;
  }
  inline void SetDeltaE(G4double deltaE){
    fDeltaE = deltaE;
  }

  inline void SetMass(G4double mass){
    fMass = mass;
  }
  inline void SetTrackID(G4int trackID){
    fTrackID = trackID;
  }
  inline void SetParentID(G4int parentID){
    fParentID = parentID;
  }

  inline void SetCopyNumSensor(G4int copyNumSensor){
    fCopyNumSensor = copyNumSensor;
  }


  /// Get hit X position
  inline G4double GetX() const { return fPosX; }
  /// Get hit Y position
  inline G4double GetY() const { return fPosY; }
  /// Get hit Z position
  inline G4double GetZ() const { return fPosZ; }
  /// Get hit time
  inline G4double GetT() const { return fT; }
  /// Get hit pdgID
  inline G4double GetPDGID() const { return fPDGID; }
  /// Get hit Energy
  inline G4double GetEnergy() const { return fEnergy; }
  /// Get hit Charge
  inline G4double GetCharge() const { return fCharge; }
  /// Get hit px
  inline G4double GetPx() const { return fPx; }
  /// Get hit py
  inline G4double GetPy() const { return fPy; }
  /// Get hit pz
  inline G4double GetPz() const { return fPz; }
  /// Get hit mass
  inline G4double GetMass() const { return fMass; }
  /// Get hit track ID
  inline G4int GetTrackID() const { return fTrackID; }
  /// Get hit parent ID
  inline G4int GetParentID() const { return fParentID; }
  /// Get change in Px
  inline G4double GetDeltaPx() const { return fDeltaPx; }
  /// Get change in Py
  inline G4double GetDeltaPy() const { return fDeltaPy; }
  /// Get change in Pz
  inline G4double GetDeltaPz() const { return fDeltaPz; }
  /// Get change in E
  inline G4double GetDeltaE() const { return fDeltaE; }
  /// Get copy number of the sensor
  inline G4int GetCopyNumSensor() const { return fCopyNumSensor; }

  inline G4ThreeVector GetTrackVertex() const {return fTrackVertex;};
  inline G4LorentzVector GetTrackP4() const {return fTrackP4;};
  inline G4int GetIsPrimaryTrack() const { return fIsPrimaryTrack; }
  inline G4int GetIsSecondaryTrack() const { return fIsSecondaryTrack; }

  inline void SetTrackVertex(const G4ThreeVector& vertex) { fTrackVertex = vertex; }
  inline void SetTrackP4(const G4LorentzVector& p4) { fTrackP4 = p4; }
  inline void SetIsPrimaryTrack(G4int isPrimary) { fIsPrimaryTrack = isPrimary; }
  inline void SetIsSecondaryTrack(G4int isSecondary) { fIsSecondaryTrack = isSecondary; }


private:
  /// Position along x axis
  G4double fPosX = -1;
  /// Position along y axis
  G4double fPosY = -1;
  /// Position along z axis
  G4double fPosZ = -1;
  /// Get time of hit
  G4double fT = -1;
  /// PDGID
  G4int fPDGID = -999;
  /// Energy
  G4double fEnergy = -999.;
  //Charge
  G4double fCharge = -999.;
  // Momentum
  G4double fPx = -999.;
  G4double fPy = -999.;
  G4double fPz = -999.;
  //Change in momentum
  G4double fDeltaPx = -999.;
  G4double fDeltaPy = -999.;
  G4double fDeltaPz = -999.;
  //Change in energy of primary
  G4double fDeltaE = -999.;
  // Mass
  G4double fMass = -999.;
  // Track ID
  G4int fTrackID = -999;
  // Parent ID
  G4int fParentID = -999;
  // Copy number of the sensor
  G4int fCopyNumSensor = -1;

  G4ThreeVector fTrackVertex{-999., -999., -999.};
  G4LorentzVector fTrackP4{0,0,0,0};
  G4int fIsPrimaryTrack = 0; // 1 if primary track, 0 otherwise
  G4int fIsSecondaryTrack = 0; // 1 if secondary track, 0 otherwise

};

typedef G4THitsCollection<FASER2TrackerHit> FASER2TrackerHitCollection;

#endif /* EXN04DETECTORHIT_HH */