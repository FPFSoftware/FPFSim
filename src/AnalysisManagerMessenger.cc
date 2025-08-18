//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "AnalysisManagerMessenger.hh"

//#include <sstream>

#include "AnalysisManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  AnalysisManagerMessenger::AnalysisManagerMessenger(AnalysisManager* manager)
  :fAnalysisManager (manager)
{
  fOutDir = new G4UIdirectory("/out/");
  fOutDir->SetGuidance("output control");

  fFileCmd = new G4UIcmdWithAString("/out/fileName", this);
  fFileCmd->SetGuidance("set name for the histograms file");
  fFileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fSaveTrackCmd = new G4UIcmdWithABool("/out/saveTrack", this);
  fSaveTrackCmd->SetGuidance("whether save the information of all tracks");
  fSaveTrackCmd->SetParameterName("saveTrack", true);
  fSaveTrackCmd->SetDefaultValue(false);

  fFLArEDir = new G4UIdirectory("/out/flare/");
  fFLArEDir->SetGuidance("flare output control");

  fSave3DEvdCmd = new G4UIcmdWithABool("/out/flare/save3DEvd", this);
  fSave3DEvdCmd->SetGuidance("whether save 3D Pixel Map");
  fSave3DEvdCmd->SetParameterName("save3DEvd", true);
  fSave3DEvdCmd->SetDefaultValue(false);

  fSave2DEvdCmd = new G4UIcmdWithABool("/out/flare/save2DEvd", this);
  fSave2DEvdCmd->SetGuidance("whether save 2D Pixel Map");
  fSave2DEvdCmd->SetParameterName("save2DEvd", true);
  fSave2DEvdCmd->SetDefaultValue(false);

  fAddDiffusionCmd = new G4UIcmdWithAString("/out/flare/addDiffusion", this);
  fAddDiffusionCmd->SetGuidance("add toy diffusion effect");
  fAddDiffusionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fPseudoRecoCmd = new G4UIcmdWithABool("/out/flare/pseudoReco", this);
  fPseudoRecoCmd->SetGuidance("whether saving pseudo reco variables");
  fPseudoRecoCmd->SetParameterName("pseudoReco", true);
  fPseudoRecoCmd->SetDefaultValue(false);

  fFASER2Dir = new G4UIdirectory("/out/faser/");
  fFASER2Dir->SetGuidance("flare output control");

  fSaveActsCmd = new G4UIcmdWithABool("/out/faser/actsHits", this);
  fSaveActsCmd->SetGuidance("save hits in Acts format");
  fSaveActsCmd->SetParameterName("actsHits", true);
  fSaveActsCmd->SetDefaultValue(true);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AnalysisManagerMessenger::~AnalysisManagerMessenger()
{
  delete fFileCmd;
  delete fSaveTrackCmd;
  delete fSave3DEvdCmd;
  delete fSave2DEvdCmd;
  delete fAddDiffusionCmd;
  delete fPseudoRecoCmd;
  delete fSaveActsCmd;
  delete fOutDir;
  delete fFLArEDir;
  delete fFASER2Dir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManagerMessenger::SetNewValue(G4UIcommand* command,G4String newValues)
{
  if (command == fFileCmd) fAnalysisManager->setFileName(newValues);
  if (command == fSaveTrackCmd) fAnalysisManager->saveTrack(fSaveTrackCmd->GetNewBoolValue(newValues));

  if (command == fSave3DEvdCmd) fAnalysisManager->save3DEvd(fSave3DEvdCmd->GetNewBoolValue(newValues));
  if (command == fSave2DEvdCmd) fAnalysisManager->save2DEvd(fSave2DEvdCmd->GetNewBoolValue(newValues));
  if (command == fPseudoRecoCmd) fAnalysisManager->savePseudoReco(fPseudoRecoCmd->GetNewBoolValue(newValues));
  if (command == fAddDiffusionCmd) fAnalysisManager->addDiffusion(newValues);

  if (command == fSaveActsCmd) fAnalysisManager->saveActs(fSaveActsCmd->GetNewBoolValue(newValues));

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
