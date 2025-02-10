#ifndef CREATENTUPLE_HH
#define CREATENTUPLE_HH

#include "G4UserRunAction.hh"
#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include <iostream>
#include <cmath>
#include <string>
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "UDInCodeControllingSurface.hh"
#include "UDDataCollection.hh"

class CreateNtuple : public G4UserRunAction {

	public:
	    
		CreateNtuple(int x,std::vector<G4String>& arr2);
		~CreateNtuple();

	private:
		UDInCodeControllingSurface ControllingSurface;	
};


#endif