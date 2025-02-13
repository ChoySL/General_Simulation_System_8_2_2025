#include <iostream>
#include <fstream>
#include <string>
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisManager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"

#include "UDInCodeControllingSurface.hh"
#include "UDEnvironmentConstruction.hh"  
#include "UDPhysicsList.hh"  
#include "UDAction.hh"





int main(int argc, char** argv) {
std::vector< G4String > Detectors_name_list;
	G4RunManager* runManager = new G4RunManager();

	runManager->SetUserInitialization(new RPCConstruction());   
	//
	G4cout << "RPC constructed " << G4endl;
	//
	runManager->SetUserInitialization(new PhysicsList());
	//
	runManager->Initialize();
	runManager->SetUserInitialization(new Action());
	//
	UDInCodeControllingSurface ControllingSurface;
	int DN1=ControllingSurface.GetNumSensitiveDetectors(Detectors_name_list);
	G4cout << "Detected detector after Initialize: " << DN1 << G4endl;
	G4UIExecutive* ui = 0;
	if (argc == 1) ui = new G4UIExecutive(argc, argv);    // argc == 1 : No parameters from command line

	G4VisManager* visManager = new G4VisExecutive("Quiet");
	visManager->Initialize();


	G4UImanager* UImanager = G4UImanager::GetUIpointer();
	if (ui) {
		UImanager->ApplyCommand("/control/execute vis.mac");
		int en;
		ControllingSurface.GetNumberOfParticles(en);
		G4String ex="/run/beamOn ";
		ex=ex+en;
		UImanager->ApplyCommand(ex);
		ui->SessionStart();
	}
	else {
		G4String command = "/control/execute ";
		G4String fileName = argv[1];
		UImanager->ApplyCommand(command + fileName);
	}
	return 0;
}
