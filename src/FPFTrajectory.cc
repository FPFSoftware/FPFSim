#include "G4VProcess.hh"
#include "G4ParticleTable.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UnitsTable.hh"
#include "G4LorentzVector.hh"

#include "FPFTrajectory.hh"

G4Allocator<FPFTrajectory> aTrajAllocator;

FPFTrajectory::FPFTrajectory()
    : fPositionRecord(0), fTrackID(0), fParentID(0),
      fPDGEncoding(0), fPDGCharge(0.0), fParticleName(""),
      fProcessName(""), fInitialP4(G4LorentzVector())
{}

FPFTrajectory::FPFTrajectory(const G4Track *aTrack, G4bool storePoints)
{
    fTrackID = aTrack->GetTrackID();
    fParentID = aTrack->GetParentID();
    fInitialP4 = aTrack->GetDynamicParticle()->Get4Momentum();

    G4ParticleDefinition *fpParticleDefinition = aTrack->GetDefinition();
    fPDGEncoding = fpParticleDefinition->GetPDGEncoding();
    fPDGCharge = fpParticleDefinition->GetPDGCharge();
    fParticleName = fpParticleDefinition->GetParticleName();

    // set creator process
    const G4VProcess *proc = aTrack->GetCreatorProcess();
    if (proc)
        fProcessName = proc->GetProcessName();
    else
        fProcessName = "primary";

    // set full trajectory storing
    fStorePoints = storePoints;

    // anyways always initialize the container
    fPositionRecord = new TrajectoryPointContainer();
    // and add the starting point (vertex)
    fPositionRecord->push_back(new G4TrajectoryPoint(aTrack->GetPosition()));
}

FPFTrajectory::FPFTrajectory(FPFTrajectory &right) : G4VTrajectory()
{
    fParticleName = right.fParticleName;
    fPDGCharge = right.fPDGCharge;
    fPDGEncoding = right.fPDGEncoding;
    fTrackID = right.fTrackID;
    fParentID = right.fParentID;
    fInitialP4 = right.fInitialP4;
    fPositionRecord = new TrajectoryPointContainer();
    fStorePoints = right.fStorePoints;
    fProcessName = right.fProcessName;

    for (size_t i = 0; i < right.fPositionRecord->size(); i++)
    {
        G4TrajectoryPoint *rightPoint = (G4TrajectoryPoint *)((*(right.fPositionRecord))[i]);
        fPositionRecord->push_back(new G4TrajectoryPoint(*rightPoint));
    }
}

FPFTrajectory::~FPFTrajectory()
{
    for (size_t i = 0; i < fPositionRecord->size(); i++)
    {
        delete (*fPositionRecord)[i];
    }
    fPositionRecord->clear();
    delete fPositionRecord;
}

void FPFTrajectory::ShowTrajectory(std::ostream &os) const
{
    G4VTrajectory::ShowTrajectory(os);
}

void FPFTrajectory::DrawTrajectory(G4int i_mode) const
{
    //G4VTrajectory::DrawTrajectory(i_mode);
    G4VTrajectory::DrawTrajectory();
}

const std::map<G4String, G4AttDef> *FPFTrajectory::GetAttDefs() const
{
    G4bool isNew;
    std::map<G4String, G4AttDef> *store = G4AttDefStore::GetInstance("FPFTrajectory", isNew);
    if (isNew)
    {
        G4String ID("ID");
        (*store)[ID] = G4AttDef(ID, "Track ID", "Bookkeeping", "", "G4int");

        G4String PID("PID");
        (*store)[PID] = G4AttDef(PID, "Parent ID", "Bookkeeping", "", "G4int");

        G4String PN("PN");
        (*store)[PN] = G4AttDef(PN, "Particle Name", "Physics", "", "G4String");

        G4String Ch("Ch");
        (*store)[Ch] = G4AttDef(Ch, "Charge", "Physics", "", "G4double");

        G4String PDG("PDG");
        (*store)[PDG] = G4AttDef(PDG, "PDG Encoding", "Physics", "", "G4int");

        G4String IMom("IMom");
        (*store)[IMom] = G4AttDef(IMom, "4Momentum of track at start of trajectory", "Physics", "", "G4LorentzVector");

        G4String NTP("NTP");
        (*store)[NTP] = G4AttDef(NTP, "No. of points", "Physics", "", "G4int");
    }
    return store;
}

std::vector<G4AttValue> *FPFTrajectory::CreateAttValues() const
{
    std::string c;
    std::ostringstream s(c);
    std::vector<G4AttValue> *values = new std::vector<G4AttValue>;

    s.seekp(std::ios::beg);
    s << fTrackID << std::ends;
    values->push_back(G4AttValue("ID", c.c_str(), ""));

    s.seekp(std::ios::beg);
    s << fParentID << std::ends;
    values->push_back(G4AttValue("PID", c.c_str(), ""));

    values->push_back(G4AttValue("PN", fParticleName, ""));

    s.seekp(std::ios::beg);
    s << fPDGCharge << std::ends;
    values->push_back(G4AttValue("Ch", c.c_str(), ""));

    s.seekp(std::ios::beg);
    s << fPDGEncoding << std::ends;
    values->push_back(G4AttValue("PDG", c.c_str(), ""));

    s.seekp(std::ios::beg);
    s << G4BestUnit(fInitialP4, "Energy") << std::ends;
    values->push_back(G4AttValue("IMom", c.c_str(), ""));

    s.seekp(std::ios::beg);
    s << GetPointEntries() << std::ends;
    values->push_back(G4AttValue("NTP", c.c_str(), ""));

    return values;
}

void FPFTrajectory::AppendStep(const G4Step *aStep)
{
    // If we are not interested in trajectory points, do nothing
    if (!fStorePoints)
        return;

    // otherwise, keep the standard behavior (store positions)
    fPositionRecord->push_back(new G4TrajectoryPoint(aStep->GetPostStepPoint()->GetPosition()));
}

G4ParticleDefinition *FPFTrajectory::GetParticleDefinition() const
{
    return (G4ParticleTable::GetParticleTable()->FindParticle(fParticleName));
}

void FPFTrajectory::MergeTrajectory(G4VTrajectory *secondTrajectory)
{
    if (!secondTrajectory) 
        return;
    FPFTrajectory *second = (FPFTrajectory *)secondTrajectory;
    G4int ent = second->GetPointEntries();
    // initial point of the second trajectory should not be merged
    for (G4int i = 1; i < ent; ++i)
    {
        fPositionRecord->push_back((*(second->fPositionRecord))[i]);
    }
    delete (*second->fPositionRecord)[0];
    second->fPositionRecord->clear();
}

G4ThreeVector FPFTrajectory::GetInitialMomentum() const 
{
    return G4ThreeVector(fInitialP4.px(), fInitialP4.py(), fInitialP4.pz());
}

G4double FPFTrajectory::GetInitialKineticEnergy() const 
{
    const G4ParticleDefinition* p = GetParticleDefinition();
    double mom = GetInitialMomentum().mag();
    if (!p) return mom;
    double mass = p->GetPDGMass();
    double kin = std::sqrt(mom*mom + mass*mass) - mass;
    if (kin<0.0) return 0.0;
    return kin;
}