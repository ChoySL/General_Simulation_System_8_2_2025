#ifndef GENERATOR_HH
#define GENERATOR_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4Geantino.hh"
#include "G4IonTable.hh"
#include "UDEnvironmentConstruction.hh"
#include "G4RootAnalysisReader.hh"
#include "G4PhysicsOrderedFreeVector.hh"
#include "g4root_defs.hh"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "TFile.h"
#include "TNtuple.h"
#include "UDInCodeControllingSurface.hh"

class muonGenerator : public G4VUserPrimaryGeneratorAction {
	public:

		muonGenerator();
		~muonGenerator();
		virtual void GeneratePrimaries(G4Event *);

	private:
		G4ParticleGun *fParticleGun;
		UDInCodeControllingSurface ControllingSurface;	

};

#endif