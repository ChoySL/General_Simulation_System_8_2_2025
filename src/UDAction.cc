#include "UDAction.hh"


Action::Action() {}
Action::~Action() {}

int GetNumSensitiveDetectors2(std::vector<G4String>& arr2) {
    G4LogicalVolumeStore* lvStore = G4LogicalVolumeStore::GetInstance();
    int numSensitiveDetectors = 0;

    for (size_t i = 0; i < lvStore->size(); i++) {
        G4LogicalVolume* logicVolume = (*lvStore)[i];
        
        if (logicVolume->GetSensitiveDetector()) {
            numSensitiveDetectors++;
			G4String name22=logicVolume->GetName();
			arr2.push_back(name22);
        }
    }
return numSensitiveDetectors;
}
void Action::Build() const {
std::vector< G4String > Detector_names;
	muonGenerator *generator = new muonGenerator();
	SetUserAction(generator);
	int DN4=GetNumSensitiveDetectors2(Detector_names);
	CreateNtuple *createNtuple = new CreateNtuple(DN4,Detector_names); 
	SetUserAction(createNtuple);
    UDRunAction* runAction = new UDRunAction();
    SetUserAction(runAction);
    UDEventAction* eventAction = new UDEventAction();
    SetUserAction(eventAction);
    UDSteppingAction* SteppingAction = new UDSteppingAction();
    SetUserAction(SteppingAction);
}