#ifndef FPFTRAJECTORY_HH
#define FPFTRAJECTORY_HH

#include "G4VTrajectory.hh"
#include "G4VTrajectoryPoint.hh"
#include "G4Allocator.hh"
#include "G4ios.hh"               
#include "G4ParticleDefinition.hh"  
#include "G4TrajectoryPoint.hh"     
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4LorentzVector.hh"

/// FPFTrajectory is a custom trajectory class that allows to disable
/// storing the full trajectory point collections (memory intensive)
/// while still bookkeeping the full hierarchial particle tree info.

typedef std::vector<G4VTrajectoryPoint *> TrajectoryPointContainer;
class FPFTrajectory : public G4VTrajectory
{
public:

    FPFTrajectory();
    // if storePoints is enabled, full trajectories are saved
    // if not (default), single steps are not saved
    FPFTrajectory(const G4Track *aTrack, G4bool storePoint = false);
    FPFTrajectory(FPFTrajectory &);
    virtual ~FPFTrajectory();

    inline void *operator new(size_t);
    inline void operator delete(void *);
    inline int operator==(const FPFTrajectory &right) const
    {
        return (this == &right);
    }

    inline G4int GetTrackID() const { return fTrackID; }
    inline G4int GetParentID() const { return fParentID; }
    inline G4String GetParticleName() const { return fParticleName; }
    inline G4double GetCharge() const { return fPDGCharge; }
    inline G4int GetPDGEncoding() const { return fPDGEncoding; }
    inline G4String GetProcessName() const { return fProcessName; }
    inline G4LorentzVector GetInitialP4() const {return fInitialP4; }
    G4ThreeVector GetInitialMomentum() const;
    G4double GetInitialKineticEnergy() const; 

    virtual void ShowTrajectory(std::ostream& os=G4cout) const;
    virtual void DrawTrajectory(G4int i_mode = 0) const;

    virtual void AppendStep(const G4Step *aStep);

    virtual int GetPointEntries() const { return fPositionRecord->size(); }
    virtual G4VTrajectoryPoint *GetPoint(G4int i) const { return (*fPositionRecord)[i]; }

    virtual void MergeTrajectory(G4VTrajectory *secondTrajectory);

    G4ParticleDefinition *GetParticleDefinition() const;

    virtual const std::map<G4String, G4AttDef> *GetAttDefs() const;
    virtual std::vector<G4AttValue> *CreateAttValues() const;

private:

    // from standard G4Trajectory implementation
    TrajectoryPointContainer *fPositionRecord;
    G4int fTrackID;
    G4int fParentID;
    G4int fPDGEncoding;
    G4double fPDGCharge;
    G4String fParticleName;
    G4LorentzVector fInitialP4;

    // custom additions
    G4String fProcessName; // creator process
    G4bool fStorePoints;   // whether to save trajectory points

};

#if defined G4TRACKING_ALLOC_EXPORT
extern G4DLLEXPORT G4Allocator<FPFTrajectory> aTrajAllocator;
#else
extern G4DLLIMPORT G4Allocator<FPFTrajectory> aTrajAllocator;
#endif

inline void* FPFTrajectory::operator new(size_t) {
    void* aTrajectory = (void*) aTrajAllocator.MallocSingle();
    return aTrajectory;
}

inline void FPFTrajectory::operator delete(void* aTrajectory) {
    aTrajAllocator.FreeSingle((FPFTrajectory*)aTrajectory);
}

#endif