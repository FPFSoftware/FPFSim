#include "generators/BackgroundGeneratorMessenger.hh"
#include "generators/BackgroundGenerator.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BackgroundGeneratorMessenger::BackgroundGeneratorMessenger(BackgroundGenerator* action) 
  : fBkgAction(action) 
{
  fBkgGeneratorDir = new G4UIdirectory("/gen/bkg/");
  fBkgGeneratorDir->SetGuidance("background generator control");

  fBkgInputFileCmd = new G4UIcmdWithAString("/gen/bkg/backgroundInput", this);
  fBkgInputFileCmd->SetGuidance("set input filename of the background generator");
  fBkgInputFileCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);

  fBkgTimeWindowCmd = new G4UIcmdWithADoubleAndUnit("/gen/bkg/backgroundWindow", this);
  fBkgTimeWindowCmd->SetGuidance("set time window for background extraction");
  fBkgTimeWindowCmd->SetUnitCategory("Time");
  fBkgTimeWindowCmd->SetUnitCandidates("s us ms ns");

  fEventOffsetCmd = new G4UIcmdWithAnInteger("/gen/bkg/eventOffset", this);
  fEventOffsetCmd->SetGuidance("set the event numbering offset");
  fEventOffsetCmd->SetDefaultValue((G4int)0);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BackgroundGeneratorMessenger::~BackgroundGeneratorMessenger()
{
  delete fBkgInputFileCmd;
  delete fBkgTimeWindowCmd;
  delete fBkgGeneratorDir;
  delete fEventOffsetCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BackgroundGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValues)
{
  if (command == fBkgInputFileCmd) fBkgAction->SetBkgFilename(newValues);
  else if (command == fBkgTimeWindowCmd) fBkgAction->SetBkgTimeWindow(fBkgTimeWindowCmd->ConvertToDimensionedDouble(newValues));
  else if (command == fEventOffsetCmd) fBkgAction->SetEventOffset(fEventOffsetCmd->GetNewIntValue(newValues));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
