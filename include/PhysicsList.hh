#ifndef PHYSICSLIST_HH
#define PHYSICSLIST_HH

#include "G4VModularPhysicsList.hh"

class G4PhysListFactory;
class G4StepLimiter;
class G4UserSpecialCuts;

class PhysicsList: public G4VModularPhysicsList
{
  public:
    PhysicsList();
    virtual ~PhysicsList();

    void AddStepMax();

    virtual void ConstructParticle();
    virtual void ConstructProcess();

  private:
    G4PhysListFactory* factory;
    G4StepLimiter* fStepLimiter;
    G4UserSpecialCuts* fUserSpecialCuts;
};

#endif
