#include "UDDataCollection.hh"
#define M_PI 3.14159265358979323846

RPCDetector::RPCDetector(G4String name) : G4VSensitiveDetector(name) {}

RPCDetector::~RPCDetector() {}

G4bool RPCDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROIhist)
{
	std::vector<G4String> Detector_names;
	int net = ControllingSurface.GetNumSensitiveDetectors(Detector_names);

	G4Track *track = aStep->GetTrack();

	G4ParticleDefinition *particleDef = track->GetDefinition();
	G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
	G4StepPoint *postStepPoint = aStep->GetPostStepPoint();

	if (preStepPoint->GetStepStatus() != fGeomBoundary)
		return true;

	G4int eventNum = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
	G4double stepLength = aStep->GetStepLength();
	G4ThreeVector posMuon = preStepPoint->GetPosition();
	G4ThreeVector momMuon = preStepPoint->GetMomentum();
	G4ThreeVector momMuonpost = postStepPoint->GetMomentum();
	const G4VTouchable *touchable = aStep->GetPreStepPoint()->GetTouchable();
	G4int copyNo = touchable->GetCopyNumber();
	G4VPhysicalVolume *physVol = touchable->GetVolume();
	G4ThreeVector posDetector = physVol->GetTranslation();
	G4double hitTime = preStepPoint->GetGlobalTime();
	G4LogicalVolume *asd = physVol->GetLogicalVolume();
	G4double muonMomentum = momMuon.mag() / 1000;
	G4double muonMomentumpost = momMuonpost.mag() / 1000;
	G4String dname = asd->GetName();
	G4ThreeVector momUD = preStepPoint->GetMomentumDirection();
	G4AnalysisManager *manager = G4AnalysisManager::Instance();

	for (int jumper = 0; jumper < net; jumper++)
	{
		if (dname == Detector_names[jumper])
		{
			manager->GetNtuple(Detector_names[jumper]);
			manager->FillNtupleIColumn((jumper + 1), 0, eventNum);
			manager->FillNtupleDColumn((jumper + 1), 1, muonMomentum);
			manager->FillNtupleDColumn((jumper + 1), 2, hitTime);
			manager->FillNtupleDColumn((jumper + 1), 3, posMuon[0] / 1000);
			manager->FillNtupleDColumn((jumper + 1), 4, posMuon[1] / 1000);
			manager->FillNtupleDColumn((jumper + 1), 5, posMuon[2] / 1000);
			manager->FillNtupleDColumn((jumper + 1), 6, posDetector[0]);
			manager->FillNtupleDColumn((jumper + 1), 7, posDetector[1]);
			manager->FillNtupleDColumn((jumper + 1), 8, posDetector[2]);
			double inclination = std::acos(momUD[2]);
			double azimuth = std::atan2(momUD[1], momUD[0]);
			manager->FillNtupleDColumn((jumper + 1), 9, inclination);
			manager->FillNtupleDColumn((jumper + 1), 10, azimuth);
			manager->FillNtupleDColumn((jumper + 1), 11, copyNo);
			manager->AddNtupleRow((jumper + 1));
		}
	}

	int Collective = net + 1;
	manager->FillNtupleIColumn(Collective, 0, eventNum);
	manager->FillNtupleDColumn(Collective, 1, muonMomentum);
	manager->FillNtupleDColumn(Collective, 2, hitTime);
	manager->FillNtupleDColumn(Collective, 3, posMuon[0]);
	manager->FillNtupleDColumn(Collective, 4, posMuon[1]);
	manager->FillNtupleDColumn(Collective, 5, posMuon[2]);
	manager->FillNtupleDColumn(Collective, 6, posDetector[0]);
	manager->FillNtupleDColumn(Collective, 7, posDetector[1]);
	manager->FillNtupleDColumn(Collective, 8, posDetector[2]);
	double inclination1 = std::acos(momUD[2]);
	double azimuth1 = std::atan2(momUD[1], momUD[0]);

	manager->FillNtupleDColumn(Collective, 9, inclination1);
	manager->FillNtupleDColumn(Collective, 10, azimuth1);
	manager->AddNtupleRow(Collective);

	if ((posMuon[1] / 1000) < 0)
	{ // changing needed*** Recevicing
		int integerPart = static_cast<int>(floor(muonMomentum)) + 1;

		int tag;
		ControllingSurface.Particle_Dealer(net, particleDef->GetParticleName(), tag);
		double incidentI;
		ControllingSurface.Getincident(incidentI);
		manager->FillNtupleIColumn(tag, 0, eventNum);
		manager->FillNtupleDColumn(tag, 1, (muonMomentum-muonMomentumpost));
		manager->FillNtupleDColumn(tag, 2, hitTime);
		manager->FillNtupleDColumn(tag, 3, posMuon[0]);
		manager->FillNtupleDColumn(tag, 4, posMuon[1]);
		manager->FillNtupleDColumn(tag, 5, posMuon[2]);
		double inclination3 = std::acos(momUD[2]);
		double azimuth3 = std::atan2(momUD[1], momUD[0]);
		manager->FillNtupleDColumn(tag, 6, inclination3);
		manager->FillNtupleDColumn(tag, 7, azimuth3);
		manager->FillNtupleDColumn(tag, 8, (incidentI - inclination3));
		manager->FillNtupleDColumn(tag, 9, stepLength);
		manager->AddNtupleRow(tag);

		// collective event modified
		std::vector<G4String> particle_list;
		ControllingSurface.Targeted_Particle_list(net, particle_list);
		int num_particles = particle_list.size();
		int fnet = Collective + num_particles + 1;
		manager->FillNtupleIColumn(fnet, 0, eventNum);
		manager->FillNtupleDColumn(fnet, 1, (muonMomentum-muonMomentumpost)); //???
		manager->FillNtupleDColumn(fnet, 2, hitTime);
		manager->FillNtupleDColumn(fnet, 3, posMuon[0]);
		manager->FillNtupleDColumn(fnet, 4, posMuon[1]);
		manager->FillNtupleDColumn(fnet, 5, posMuon[2]);
		manager->FillNtupleDColumn(fnet, 6, inclination3);
		manager->FillNtupleDColumn(fnet, 7, azimuth3);
		manager->FillNtupleDColumn(tag, 8, (incidentI - inclination3));
		manager->FillNtupleDColumn(tag, 9, stepLength);
		ControllingSurface.GetHit(particleDef->GetParticleName());
		manager->AddNtupleRow(fnet);
	}

	return true;
}
