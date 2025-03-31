#ifndef TRACKINFORMATION_HH
#define TRACKINFORMATION_HH

#include "G4Track.hh"
#include "G4Allocator.hh"
#include "G4VUserTrackInformation.hh"

class TrackInformation : public G4VUserTrackInformation
{
  public:
    TrackInformation();
    TrackInformation(const G4Track* aTrack);
    //TrackInformation(const TrackInformation* aTrackInfo);
    virtual ~TrackInformation();

    inline void *operator new(size_t);
    inline void operator delete(void *aTrackInfo);

    TrackInformation& operator =(const TrackInformation& right);

    virtual void Print() const;

  public:
    inline void SetTrackIsFromPrimaryPizero(G4int i) {fFromPrimaryPizero = i;}
    inline void SetTrackIsFromFSLPizero(G4int i) {fFromFSLPizero = i;}
    inline void SetTrackIsFromPrimaryLepton(G4int i) {fFromPrimaryLepton = i;}
    inline G4int IsTrackFromPrimaryPizero() const {return fFromPrimaryPizero;}
    inline G4int IsTrackFromFSLPizero() const {return fFromFSLPizero;}
    inline G4int IsTrackFromPrimaryLepton() const {return fFromPrimaryLepton;}

  private:
    G4int fFromPrimaryPizero;
    G4int fFromFSLPizero;
    G4int fFromPrimaryLepton;
};

extern G4ThreadLocal
 G4Allocator<TrackInformation> * aTrackInformationAllocator;

inline void* TrackInformation::operator new(size_t)
{
  if(!aTrackInformationAllocator)
    aTrackInformationAllocator = new G4Allocator<TrackInformation>;
  return (void*)aTrackInformationAllocator->MallocSingle();
}

inline void TrackInformation::operator delete(void *aTrackInfo)
{ aTrackInformationAllocator->FreeSingle((TrackInformation*)aTrackInfo);}

#endif
