#include "UDparticlegenerator.hh"

double out0, out1, out2, out3, out4, out5, out6;
int entries;
int counter = 0;
G4double momGenerate = 0;
std::vector<int> arr;
std::vector<double> arre;
std::vector<double> arrpz;
std::vector<double> arrpy;
std::vector<double> arrpx;
std::vector<double> arrdx;
std::vector<double> arrdy;
std::vector<double> arrdz;

muonGenerator::muonGenerator()
{
	G4String input;
	ControllingSurface.Getfilename("generator", input);
	TFile file2(input.c_str());

	if (!file2.IsOpen() || file2.IsZombie())
	{
		std::cerr << "Error opening file." << std::endl;
		return;
	}
	TNtuple *ntuple = dynamic_cast<TNtuple *>(file2.Get("ntuple"));

	if (!ntuple)
	{
		std::cerr << "Error: NTuple not found in file." << std::endl;
		file2.Close();
		return;
	}
	printf("UDParticleGenerator reading data from: ");
	printf(input.c_str());
	printf("\n");
	entries = ntuple->GetEntries();
	fParticleGun = new G4ParticleGun(1); // number of particles per event (several events = 1 run)
	G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName = "mu-";
	G4ParticleDefinition *particle = particleTable->FindParticle(particleName);
	fParticleGun->SetParticleDefinition(particle);
	Float_t e, px, py, pz, dx, dy, dz;

	G4cout << "Root entries: " << entries << G4endl; // data reading here
	ntuple->SetBranchAddress("e", &e);
	ntuple->SetBranchAddress("px", &px);
	ntuple->SetBranchAddress("py", &py);
	ntuple->SetBranchAddress("pz", &pz);
	ntuple->SetBranchAddress("dx", &dx);
	ntuple->SetBranchAddress("dy", &dy);
	ntuple->SetBranchAddress("dz", &dz);
	for (int k = 0; k < entries; k++)
	{
		ntuple->GetEntry(k);
		out0 = e;
		arre.push_back(out0);
		out1 = px;
		arrpx.push_back(out1);
		out2 = py;
		arrpy.push_back(out2);
		out3 = pz;
		arrpz.push_back(out3);
		out4 = dx;
		arrdx.push_back(out4);
		out5 = dy;
		arrdy.push_back(out5);
		out6 = dz;
		arrdz.push_back(out6);
	}

	printf("Please Run \n/run/beamOn %d\n", entries);
	printf("In the GEANT4 terminal\n");

	file2.Close();
}

muonGenerator::~muonGenerator()
{
	delete fParticleGun;
}

void muonGenerator::GeneratePrimaries(G4Event *anEvent)
{

	const RPCConstruction *rpcConstruction = static_cast<const RPCConstruction *>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
	const double PI = 3.14159265358979323846;

	G4double xRdmGenerate = G4UniformRand() * rpcConstruction->getxRPCFull() - rpcConstruction->getxRPCFull() / 2;
	G4double yHighestRPC = rpcConstruction->getyHighestRPC();
	G4double zRdmGenerate = G4UniformRand() * rpcConstruction->getzRPCFull() - rpcConstruction->getzRPCFull() / 2;

	G4double thetaRdmGenerate = G4UniformRand() * PI / 4;
	G4double phiRdmGenerate = G4UniformRand() * 2 * PI;
	G4double xMomRdmGenerate = std::sin(thetaRdmGenerate) * std::cos(phiRdmGenerate);
	G4double yMomRdmGenerate = fabs(std::cos(thetaRdmGenerate));
	G4double zMomRdmGenerate = std::sin(thetaRdmGenerate) * std::sin(phiRdmGenerate);

	G4double xPos = xRdmGenerate + 10 * m * xMomRdmGenerate;
	G4double yPos = yHighestRPC + 10 * m * yMomRdmGenerate;
	G4double zPos = zRdmGenerate + 10 * m * zMomRdmGenerate;

	G4AnalysisManager *manager = G4AnalysisManager::Instance();
	G4double hitTime = 0.0;
	int EventNumber = anEvent->GetEventID();
	G4ThreeVector pos((arrpx[EventNumber]) * m, (arrpy[EventNumber]) * m, (arrpz[EventNumber]) * m); // starting pos
	fParticleGun->SetParticlePosition(pos);
	G4double direction_x = arrdx[EventNumber];
	G4double direction_y = arrdy[EventNumber];
	G4double direction_z = arrdz[EventNumber];
	G4ThreeVector mom(direction_x, direction_y, direction_z); // firing direction
	fParticleGun->SetParticleMomentumDirection(mom);
	fParticleGun->SetParticleMomentum((arre[EventNumber]) * GeV);
	fParticleGun->GeneratePrimaryVertex(anEvent);

	manager->FillNtupleIColumn(0, 0, EventNumber);
	manager->FillNtupleDColumn(0, 1, (arre[EventNumber]) / 1000);
	manager->FillNtupleDColumn(0, 2, hitTime);
	manager->FillNtupleDColumn(0, 3, arrpx[EventNumber]);
	manager->FillNtupleDColumn(0, 4, arrpy[EventNumber]);
	manager->FillNtupleDColumn(0, 5, arrpz[EventNumber]);
	double inclination0 = std::acos(direction_z);
	double azimuth0 = std::atan2(direction_y, direction_x);
	manager->FillNtupleDColumn(0, 6, inclination0);
	manager->FillNtupleDColumn(0, 7, azimuth0);
	manager->AddNtupleRow(0);
}