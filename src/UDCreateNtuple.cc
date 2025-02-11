#include "UDcreateNtuple.hh"

int GetNumSensitiveDetectors20(std::vector<G4String>& arr2) {
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



CreateNtuple::CreateNtuple(int x ,std::vector< G4String >& Detector_names) {
G4AnalysisManager *manager = G4AnalysisManager::Instance();
for(int i=0;i<=x+1;i++){
	if(i==0){
	manager->CreateNtuple("generator", "generator");
	manager->CreateNtupleIColumn("ParticleNumber");
	manager->CreateNtupleDColumn("muonMomentum");
	manager->CreateNtupleDColumn("hitTime");
	manager->CreateNtupleDColumn("hitPositionX_truth");//3, mm
	manager->CreateNtupleDColumn("hitPositionY_truth");//4
	manager->CreateNtupleDColumn("hitPositionZ_truth");//5
	manager->CreateNtupleDColumn("Inclination [0,Pi]");//6
	manager->CreateNtupleDColumn("Azimuth [0,2Pi)");//7
	manager->FinishNtuple(0);

	}else if (i>0&&i<=x){

	
	//detectors
	int id=3000+i;
	G4String NtupleName = "Event received from_"+Detector_names[i-1];
	G4String NtupleNameo = to_string(i-1);
	G4String NtupleID = Detector_names[i-1];
	manager->CreateNtuple(NtupleID , NtupleName);
	manager->CreateNtupleIColumn("eventNumber"); //0
	manager->CreateNtupleDColumn("muonMomentum");  //1, GeV
	manager->CreateNtupleDColumn("hitTime");  //2, ns
	manager->CreateNtupleDColumn("hitPositionX_truth");  //3, mm
	manager->CreateNtupleDColumn("hitPositionY_truth");//4
	manager->CreateNtupleDColumn("hitPositionZ_truth");//5
	manager->CreateNtupleDColumn("hitPixelX");//6
	manager->CreateNtupleDColumn("hitPixelY");//7
	manager->CreateNtupleDColumn("hitPixelZ");//8
	manager->CreateNtupleDColumn("Inclination [0,Pi]");//6
	manager->CreateNtupleDColumn("Azimuth [0,2Pi)");//7
	manager->CreateNtupleDColumn("Hitted strip");//12
	manager->FinishNtuple(i-1);

}else if(i==x+1){
	manager->CreateNtuple("Collective Event" , "events");
	manager->CreateNtupleIColumn("eventNumber"); //0
	manager->CreateNtupleDColumn("muonMomentum");  //1, GeV
	manager->CreateNtupleDColumn("hitTime");  //2, ns
	manager->CreateNtupleDColumn("hitPositionX_truth");  //3, mm
	manager->CreateNtupleDColumn("hitPositionY_truth");//4
	manager->CreateNtupleDColumn("hitPositionZ_truth");//5
	manager->CreateNtupleDColumn("hitPixelX");//6
	manager->CreateNtupleDColumn("hitPixelY");//7
	manager->CreateNtupleDColumn("hitPixelZ");//8
	manager->CreateNtupleDColumn("hitDirectionX");//9
	manager->CreateNtupleDColumn("hitDirectionY");//10
	manager->CreateNtupleDColumn("hitDirectionZ");//11
	manager->FinishNtuple(i);


}
}

std::vector<G4String> particle_list;
ControllingSurface.Targeted_Particle_list(x,particle_list);
int num_particles = particle_list.size();
for(int k=1;k<=num_particles;k++){
	G4String sepname=particle_list[(k-1)]+"_Collective Event modified";

	manager->CreateNtuple(sepname , particle_list[(k-1)]);
	manager->CreateNtupleIColumn("eventNumber");//0
	manager->CreateNtupleDColumn("muonMomentum_diff");//1
	manager->CreateNtupleDColumn("hitTime");  //2, ns
	manager->CreateNtupleDColumn("hitPositionX_truth");//3, mm
	manager->CreateNtupleDColumn("hitPositionY_truth");//4
	manager->CreateNtupleDColumn("hitPositionZ_truth");//5
	manager->CreateNtupleDColumn("Inclination [0,Pi]");//6
	manager->CreateNtupleDColumn("Azimuth [-Pi,Pi)");//7
	manager->CreateNtupleDColumn("Inclination_diff [0,Pi]");//8
	manager->CreateNtupleDColumn("Step Length");
	manager->FinishNtuple(x+1+k);
}
	manager->CreateNtuple("Collective Event Modified" , "events");
	manager->CreateNtupleIColumn("eventNumber");//0
	manager->CreateNtupleDColumn("muonMomentum_diff");//1
	manager->CreateNtupleDColumn("hitTime");  //2, ns
	manager->CreateNtupleDColumn("hitPositionX_truth");//3, mm
	manager->CreateNtupleDColumn("hitPositionY_truth");//4
	manager->CreateNtupleDColumn("hitPositionZ_truth");//5
	manager->CreateNtupleDColumn("Inclination [0,Pi]");//6
	manager->CreateNtupleDColumn("Azimuth [-Pi,Pi)");//7
	manager->CreateNtupleDColumn("Inclination_diff [0,Pi]");//8
	manager->CreateNtupleDColumn("Step Length");
	manager->FinishNtuple(x+2+num_particles);

	manager->CreateNtuple("Primany Particle Stepping" , "events");
	manager->CreateNtupleIColumn("eventNumber");//0
	manager->CreateNtupleIColumn("NumberofStep");//1
	manager->CreateNtupleDColumn("muonMomentum_diff");//2
	manager->CreateNtupleDColumn("hitTime");  //3, ns
	manager->CreateNtupleDColumn("hitPositionX_truth");//4, mm
	manager->CreateNtupleDColumn("hitPositionY_truth");//5
	manager->CreateNtupleDColumn("hitPositionZ_truth");//6
	manager->CreateNtupleDColumn("Inclination [0,Pi]");//7
	manager->CreateNtupleDColumn("Azimuth [-Pi,Pi)");//8
	manager->CreateNtupleDColumn("Inclination_diff [0,Pi]");//9
	manager->CreateNtupleDColumn("Azimuth_diff [-Pi,Pi)");//10
	manager->CreateNtupleIColumn("Particle Type Number");//11
	manager->CreateNtupleDColumn("EnergyDeposit");//12
	manager->CreateNtupleDColumn("StepLength");//13
	manager->CreateNtupleDColumn("DeltaEnergy");//14
	manager->FinishNtuple(x+3+num_particles);


}

CreateNtuple::~CreateNtuple() {}

