#include <G4RunManager.hh>
#include <G4UImanager.hh>

#include <G4VisExecutive.hh>

#include <G4UIExecutive.hh>

#include <G4String.hh>
#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "AnalysisManager.hh"
#include "PhysicsList.hh"

using namespace std;

/* Main function that enables to:
 * - run macros
 * - start interactive UI mode (no arguments)
 */
int main(int argc, char** argv) {
  std::cout<<"Application starting..."<<std::endl;
//  G4long myseed = 345354;
//  CLHEP::HepRandom::setTheSeed(myseed);

  // invoke analysis manager before ui manager to invoke analysis manager messenger
  AnalysisManager* analysis = AnalysisManager::GetInstance();

  // Create the run manager (MT or non-MT) and make it a bit verbose.
  auto runManager = new G4RunManager();
  runManager->SetVerboseLevel(1);

  // Set mandatory initialization classes
  runManager->SetUserInitialization(new DetectorConstruction());
  runManager->SetUserInitialization(new PhysicsList());

  // Set user action classes
  runManager->SetUserInitialization(new ActionInitialization());

  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive();
  visManager->SetVerboseLevel(1);   // Default, you can always override this using macro commands
  visManager->Initialize();

  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Parse command line arguments
  if (argc==1) {
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
    UImanager->ApplyCommand("/control/execute macros/vis.mac");
    ui->SessionStart();
    delete ui;
  } else {
    //G4UIExecutive* ui = new G4UIExecutive(argc, argv);
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
    if (argc==3) {
      G4String mode = argv[2];
      if (mode == "vis") {
        G4UIExecutive* ui = new G4UIExecutive(argc, argv);
        ui->SessionStart();
        delete ui;
      } else {
        std::cout<<"Please specify the second argument as vis to visualize the event"<<std::endl;
        return 0;
      }
    }
  }

  delete visManager;
  delete runManager;

  std::cout<<"Application sucessfully ended.\nBye :-)"<<std::endl;

  return 0;
}
