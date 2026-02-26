#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGeneratorAction.hh"

#include "generators/GeneratorBase.hh"
#include "generators/GENIEGenerator.hh"
#include "generators/BackgroundGenerator.hh"
#include "generators/HepMCGenerator.hh"
#include "generators/GPSGenerator.hh"

#include "geometry/GeometricalParameters.hh"

#include "EventInformation.hh"

#include "G4Event.hh"
#include "G4Exception.hh"


PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  // create a messenger for this class
  fGenMessenger = new PrimaryGeneratorMessenger(this);

  // start with default generator
  fGenerator = new GPSGenerator();
  fInitialized = false;

}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fGenerator;
  delete fGenMessenger;
}

void PrimaryGeneratorAction::SetGenerator(G4String name)
{
  G4StrUtil::to_lower(name);

  if( name == "genie" )
    fGenerator = new GENIEGenerator();
  else if( name == "background" )
    fGenerator = new BackgroundGenerator();
  else if( name == "hepmc" )
    fGenerator = new HepMCGenerator();
  else if ( name == "gun" )
    fGenerator = new GPSGenerator();
  else{
    G4String err = "Unknown generator option " + name;
    G4Exception("PrimaryGeneratorAction",
                "UnknownOption",
                FatalErrorInArgument,
                err.c_str());
  }
}


void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // load generator data at first event
  // this function opens files, reads trees, etc (if required)
  if(!fInitialized){
    fGenerator->LoadData();
    fInitialized = true;
  }

  G4cout << G4endl;
  G4cout << "===oooOOOooo=== Event Generator (# " << anEvent->GetEventID();

  // reset event metadata
  fGenerator->ResetEventMetadata();

  // produce an event with current generator
  fGenerator->GeneratePrimaries(anEvent);

  // save vertex metadata information into the event
  anEvent->SetUserInformation(new EventInformation(fGenerator->GetEventMetadata()));

}
