#ifndef LARBOXHIT_HH
#define LARBOXHIT_HH

#include <G4VHit.hh>
#include <G4THitsMap.hh>
#include <G4ThreeVector.hh>

/**
 * Custom hit class
 */
class LArBoxHit : public G4VHit {
  public:
    // Memory allocation and de-allocation
    inline void* operator new(size_t);
    inline void operator delete(void*);

    //setter
    void SetTrackStatus(G4int trackstatus) { fTrackStatus = trackstatus; }
    void SetTrackVertex(G4ThreeVector& TrackVertex) { fTrackVertex = TrackVertex; }
    void SetTrackLength(G4double TrackLength)       { fTrackLength = TrackLength; }
    void SetParticle(G4int particle) { fParticle = particle; }
    void SetParticleMass(G4double mass) { fParticleMass = mass; }
    void SetPID(G4int pid) { fPID = pid; }
    void SetTID(G4int tid) { fTID = tid; }
    void SetStepNo(G4int stepno) { fStepno = stepno; }
    void SetPreStepPosition(G4ThreeVector& PreStepPosition) { fPreStepPosition = PreStepPosition; }
    void SetPostStepPosition(G4ThreeVector& PostStepPosition) { fPostStepPosition = PostStepPosition; }
    void SetInitMomentum(G4ThreeVector& InitMomentum) { fInitMomentum = InitMomentum; }
    void SetInitKinEnergy(G4double InitKinEne) { fInitKinEne = InitKinEne; }
    void SetCreatorProcess(G4String processname) { fCreatorProcess = processname; }
    void SetProcessName(G4String processname) { fProcessName = processname; }
    void SetStepLength(G4double steplength) { fStepLength = steplength; }
    void SetEdep(G4double edep) { fEdep = edep; }
    void SetEdepPosition(G4ThreeVector& edepPos) { fEdepPosition = edepPos; }
    void SetPreVolume(G4String volumeName) { fPreVolumeName = volumeName; }
    void SetCopyNumPreVolume(G4int copyNumber) { fCopyNumPreVol = copyNumber; }
    void SetPostVolume(G4String volumeName) { fPostVolumeName = volumeName; }
    void SetCopyNumPostVolume(G4int copyNumber) { fCopyNumPostVol = copyNumber; }
    void SetStepStatus(G4int stepstatus) { fStepStatus = stepstatus; }
    void SetTrackIsFromPrimaryPizero(G4int i) { fTrackIsFromPrimaryPizero = i; }
    void SetTrackIsFromFSLPizero(G4int i) { fTrackIsFromFSLPizero = i; }
    void SetTrackIsFromPrimaryLepton(G4int i) { fTrackIsFromPrimaryLepton = i; }

    //getter
    G4int GetTrackStatus() const { return fTrackStatus; }
    G4ThreeVector GetTrackVertex() const { return fTrackVertex; }
    G4double GetTrackLength() const { return fTrackLength; }
    G4int GetParticle() const { return fParticle; }
    G4double GetParticleMass() const { return fParticleMass; }
    G4int GetPID() const { return fPID; }
    G4int GetTID() const { return fTID; }
    G4int GetStepNo() const { return fStepno; }
    G4ThreeVector GetPreStepPosition() const { return fPreStepPosition; }
    G4ThreeVector GetPostStepPosition() const { return fPostStepPosition; }
    G4ThreeVector GetInitMomentum()    const { return fInitMomentum; }
    G4double GetInitKinEnergy() const { return fInitKinEne; }
    G4String GetCreatorProcess() const { return fCreatorProcess; }
    G4String GetProcessName() const { return fProcessName; }
    G4double GetStepLength() const { return fStepLength; }
    G4double GetEdep() const { return fEdep; }
    G4ThreeVector GetEdepPosition() const { return fEdepPosition; }
    G4String GetPostVolume() const { return fPostVolumeName; }
    G4int GetCopyNumPostVolume() const { return fCopyNumPostVol; }
    G4String GetPreVolume() const { return fPreVolumeName; }
    G4int GetCopyNumPreVolume() const { return fCopyNumPreVol; }
    G4int GetStepStatus() const { return fStepStatus; }
    G4int GetIsTrackFromPrimaryPizero() const { return fTrackIsFromPrimaryPizero; }
    G4int GetIsTrackFromFSLPizero() const { return fTrackIsFromFSLPizero; }
    G4int GetIsTrackFromPrimaryLepton() const { return fTrackIsFromPrimaryLepton; }

  private:
    G4int fTrackStatus;
    G4ThreeVector fTrackVertex;
    G4double fTrackLength;
    G4int fParticle;
    G4double fParticleMass;
    G4int fPID;
    G4int fTID;
    G4int fStepno;
    G4ThreeVector fPreStepPosition;
    G4ThreeVector fPostStepPosition;
    G4ThreeVector fInitMomentum;
    G4double fInitKinEne;
    G4String fCreatorProcess;
    G4String fProcessName;
    G4double fStepLength;
    G4double fEdep;
    G4ThreeVector fEdepPosition;
    G4String fPreVolumeName;
    G4int fCopyNumPreVol;
    G4String fPostVolumeName;
    G4int fCopyNumPostVol;
    G4int fStepStatus;
    G4bool fHitFromFSL;
    G4int fTrackIsFromPrimaryPizero;
    G4int fTrackIsFromFSLPizero;
    G4int fTrackIsFromPrimaryLepton;
};

using LArBoxHitsCollection = G4THitsCollection<LArBoxHit>;

extern G4ThreadLocal G4Allocator<LArBoxHit> *hitAllocator;

inline void* LArBoxHit::operator new(size_t) {
  if (!hitAllocator) {
    hitAllocator = new G4Allocator<LArBoxHit>;
  }
  return hitAllocator->MallocSingle();
}

inline void LArBoxHit::operator delete(void* aHit) {
  if (!hitAllocator) {
    hitAllocator = new G4Allocator<LArBoxHit>;
  }
  hitAllocator->FreeSingle((LArBoxHit*) aHit);
}

#endif
