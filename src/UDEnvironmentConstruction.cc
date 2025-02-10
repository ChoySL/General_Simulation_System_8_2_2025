#include "UDenvironmentconstruction.hh"

RPCConstruction::RPCConstruction()
{

	xWorldFull = 50 * m;
	yWorldFull = 50 * m;
	zWorldFull = 50 * m;

	xplateFull = 2 * m;
	yplateFull = 2 * m;
	zplateFull = 2 * m;

	xTruckFull = 12.19 * m;
	yTruckFull = 2.44 * m;
	zTruckFull = 2.59 * m;
	wheelRadius = 0.5 * m;

	xDetectorplateFull = 50 * m;
	yDetectorplateFull = 50 * m;
	zDetectorplateFull = 50 * m;

	xPixelFull = 6 * cm;
	yPixelFull = 1 * cm;
	zPixelFull = 12 * cm;
	numX = 155;
	numZ = 155;
	Xstrip = 600;
	Zstrip = 300;
	// heights of RPC plates from top/bottom of truck
	// height = {-3*m, 3*m};
	height = {-4 * m, -2 * m, 2 * m, 4 * m};

	DefineMaterials();
}
RPCConstruction::~RPCConstruction() {}

void RPCConstruction::DefineMaterials()
{

	G4NistManager *nist = G4NistManager::Instance();

	worldMat = nist->FindOrBuildMaterial("G4_AIR");
	C2H2F4 = new G4Material("C2H2F4", 1.977 * kg / m3, 3);
	C2H2F4->AddElement(nist->FindOrBuildElement("C"), 2);
	C2H2F4->AddElement(nist->FindOrBuildElement("H"), 2);
	C2H2F4->AddElement(nist->FindOrBuildElement("F"), 4);

	C4H10 = new G4Material("C4H10", 1.977 * kg / m3, 2);
	C4H10->AddElement(nist->FindOrBuildElement("C"), 4);
	C4H10->AddElement(nist->FindOrBuildElement("H"), 10);

	SF6 = new G4Material("SF6", 1.977 * kg / m3, 2);
	SF6->AddElement(nist->FindOrBuildElement("S"), 1);
	SF6->AddElement(nist->FindOrBuildElement("F"), 6);

	RPCG = new G4Material("RPCG", 183 * mg / cm3, 3);
	RPCG->AddMaterial(C2H2F4, 94.7 * perCent);
	RPCG->AddMaterial(C4H10, 5.0 * perCent);
	RPCG->AddMaterial(SF6, 0.3 * perCent);

	ammoniumNitrate = new G4Material("ammoniumNitrate", 1.72 * g / cm3, 3);
	ammoniumNitrate->AddElement(nist->FindOrBuildElement("N"), 2);
	ammoniumNitrate->AddElement(nist->FindOrBuildElement("H"), 4);
	ammoniumNitrate->AddElement(nist->FindOrBuildElement("O"), 3);
}

G4VPhysicalVolume *RPCConstruction::Construct()
{
	G4bool checkOverlaps = true;
	G4double DThickness = 0.5 * m;
	G4double Thickness = 1 * m;

	// Create world
	solidWorld = new G4Box("solidWorld", xWorldFull / 2, yWorldFull / 2, zWorldFull / 2);
	logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");
	physWorld = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicWorld, "physWorld", 0, false, 0, true);

	// Create truck

	solidTruck = new G4Box("solidTruck", xTruckFull / 2, yTruckFull / 2, zTruckFull / 2);
	logicTruck = new G4LogicalVolume(solidTruck, ammoniumNitrate, "logicTruck");
	physTruck = new G4PVPlacement(0, G4ThreeVector(0., -1 * m, 0.), logicTruck, "physTruck", logicWorld, false, 2, true);

	G4Box *soildplateup = new G4Box("soildplateup", 0.5 * xWorldFull, 0.5 * DThickness, 0.5 * yWorldFull);
	logicplateup = new G4LogicalVolume(soildplateup, worldMat, "logicplateup"); // RPCG
	G4VPhysicalVolume *physplateup = new G4PVPlacement(0, G4ThreeVector(0., 20 * m, 0.), logicplateup, "physplateup", logicWorld, false, 0, true);

	G4Box *soildplatedown = new G4Box("soildplatedown", 0.5 * xWorldFull, 0.5 * DThickness, 0.5 * yWorldFull); // 0.0001*DThickness
	logicplatedown = new G4LogicalVolume(soildplatedown, RPCG, "logicplatedown");
	G4VPhysicalVolume *physplatedown = new G4PVPlacement(0, G4ThreeVector(0., -20 * m, 0.), logicplatedown, "physplatedown", logicWorld, false, 0, true);

	/*
		G4Box *soildplatedup= new G4Box("soildplatedup", 0.5*xWorldFull, 0.5*DThickness,0.125*yWorldFull);
		logicplatedup = new G4LogicalVolume(soildplatedup, C2H2F4, "logicplatedup");
		G4VPhysicalVolume *physplatedup = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicplatedup, "physplatedup", logicplateup, true,0, true);
	*/
	/*
		//creating the e-field for plateup
		G4ElectricField* electricField = new G4UniformElectricField(G4ThreeVector(0.0, 4.5 * m, 0.0));
		G4EqMagElectricField* equation = new G4EqMagElectricField(electricField);
		G4FieldManager* fieldManager = new G4FieldManager(electricField);
		logicplateup->SetFieldManager(fieldManager, true);
	*/

	// Create a detector plate downside

	/*
	 */
	/*
		// Create wheel
		G4Rotate3D rotY(90*deg, G4ThreeVector(0,1,0));
		G4Translate3D trans1(G4ThreeVector(0., -yTruckFull/2-wheelRadius, zTruckFull/2-1*m));
		G4Translate3D trans2(G4ThreeVector(0., -yTruckFull/2-wheelRadius, -zTruckFull/2+1*m));
		solidWheel = new G4Tubs("solidWheel", 0., wheelRadius, xTruckFull/2, 0*deg, 360*deg);
		logicWheel = new G4LogicalVolume(solidWheel, worldMat, "logicWheel");
		physWheel = new G4PVPlacement((trans1)*(rotY), logicWheel, "physWheel", logicWorld, false, 0, true);
		physWheel = new G4PVPlacement((trans2)*(rotY), logicWheel, "physWheel", logicWorld, false, 0, true);


	// Creating the RPC Plate
		solidPlate = new G4Box("solidPlate", xPlateFull, yPlateFull, zPlateFull);
		logicPlate = new G4LogicalVolume(solidPixel, worldMat, "logicPlate");
		G4Translate3D transRPC(G4ThreeVector(xTrans, yTrans, zTrans));
		physPlate = new G4PVPlacement(transRPC, logicPixel, "physPixel", logicWorld, false, numZ*numX*k+numX*i+j, false);
	*/

	/*
		// Create RPC "pixels"
		solidPixel = new G4Box("solidPixel", xPixelFull/2, yPixelFull/2, zPixelFull/2);
		logicPixel = new G4LogicalVolume(solidPixel, worldMat, "logicPixel");

		for (G4int i = 0; i < numZ; i++) {
			for (G4int j = 0; j < numX; j++) {
				G4double xTrans = -(numX-1)/2.*xPixelFull + j * xPixelFull;
				G4double zTrans = -(numZ-1)/2.*zPixelFull + i * zPixelFull;
			for (G4int k = 0; k < 1; k++) {
				G4double yTrans = (height[k] > 0 ? yTruckFull/2+height[k] : -yTruckFull/2+height[k]);
				G4Translate3D transRPC(G4ThreeVector(xTrans, yTrans, zTrans));
				physPixel = new G4PVPlacement(transRPC, logicPixel, "physPixel", logicWorld, false, numZ*numX*k+numX*i+j, false);
				}
			}
		}
	*/

	return physWorld;
}

void RPCConstruction::ConstructSDandField()
{

	RPCDetector *rpcdetup = new RPCDetector("rpcdetup");
	SetSensitiveDetector(logicplateup, rpcdetup);
	RPCDetector *rpcdetdown = new RPCDetector("rpcdetdown");
	logicplatedown->SetSensitiveDetector(rpcdetdown);
}
