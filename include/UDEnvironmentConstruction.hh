#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH
#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "UDDataCollection.hh"
#include "G4GenericMessenger.hh"
#include "G4VisAttributes.hh"
#include "G4Color.hh"
#include "G4SDManager.hh"
#include <vector>
#include "G4ElectricField.hh"
#include "G4UniformElectricField.hh"
#include "G4EqMagElectricField.hh"
#include "G4FieldManager.hh"
#include "G4Electron.hh"
#include "G4ParticleGun.hh"
#include "G4FieldBuilder.hh"
#include "CLHEP/Units/SystemOfUnits.h"

class RPCConstruction : public G4VUserDetectorConstruction {
	public:
		RPCConstruction();
		~RPCConstruction();
		virtual G4VPhysicalVolume *Construct();   
		G4double getxWorldFull() const { return xWorldFull; }
		G4double getyWorldFull() const { return yWorldFull; }
		G4double getzWorldFull() const { return zWorldFull; }
		G4double getxRPCFull() const { return numX*xPixelFull; }
		G4double getzRPCFull() const { return numZ*zPixelFull; }
		G4double getyHighestRPC() const { return yTruckFull/2 + (*std::max_element(height.begin(), height.end())); }
	private:
		G4double xWorldFull, yWorldFull, zWorldFull;
		G4double xplateFull, yplateFull, zplateFull,xDetectorplateFull,yDetectorplateFull,zDetectorplateFull;
		G4double xTruckFull, yTruckFull, zTruckFull;
		G4double wheelRadius;
		G4double xPixelFull, yPixelFull, zPixelFull;
		std::vector<G4double> height;
		G4int numX, numZ,Xstrip,Zstrip;

		G4Tubs *solidWheel;
		G4LogicalVolume *logicWheel;
		G4VPhysicalVolume *physWheel;

		G4Box *solidWorld,*soildstripsup, *solidTruck, *solidPixel, *soildplate, *rpc;
		G4LogicalVolume *logicWorld, *logicTruck, *logicPixel,*logicstripsup,*logicplateup,*logicplatedown, *logicplateup1stlayerx_dir,*logicplateup1stlayerz_dir,*logicplateup1stlayerz_dir_divider,*logicplateup1stlayerx_dir_divider;
		G4VPhysicalVolume *physWorld, *physTruck, *physPixel,*physstripsup,*physplatedown,*physplateup, *physplateup1stlayerx_dir,*physplateup1stlayerz_dir,*physplateup1stlayerx_dir_divider,*physplateup1stlayerz_dir_divider;

		void DefineMaterials();
		G4Material *worldMat, *ammoniumNitrate, *RPCG, *Weathering_Steel,*C2H2F4,*SF6,*C4H10;
		virtual void ConstructSDandField();
		
};
#endif
