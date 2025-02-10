#include "UDRunAction.hh"

void UDRunAction::BeginOfRunAction(const G4Run *run)
{
	G4UImanager *UImanager = G4UImanager::GetUIpointer();

	G4AnalysisManager *manager = G4AnalysisManager::Instance();
	G4int runID = run->GetRunID();

	std::stringstream strRunID;
	strRunID << runID;
	G4String fileName;
	ControllingSurface.Getfilename("Ntuple_temp", fileName);
	manager->OpenFile(fileName);
}

void UDRunAction::EndOfRunAction(const G4Run *run)
{
	std::vector<G4String> outtypelist;
	G4String fileName;
	G4AnalysisManager *manager = G4AnalysisManager::Instance();
	manager->Write();
	manager->CloseFile();
	ControllingSurface.Getfilename("Ntuple_Output", fileName);
	ControllingSurface.Getdeleted(outtypelist, fileName);
	ControllingSurface.CreateHistogram(fileName);
	int revt = run->GetNumberOfEventToBeProcessed();
	printf("Recorded Events in this run: %d\n", revt);
}