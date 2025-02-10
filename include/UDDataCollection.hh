#ifndef DETECTOR_HH
#define DETECTOR_HH

#include <iostream>
#include <cmath>
#include <math.h>
#include "G4VSensitiveDetector.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicsOrderedFreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "UDInCodeControllingSurface.hh"
# define M_PI 3.14159265358979323846
class RPCDetector: public G4VSensitiveDetector {

	public:
		RPCDetector(G4String);
		~RPCDetector();
	private:
		virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *);
		UDInCodeControllingSurface ControllingSurface;	
};

#endif