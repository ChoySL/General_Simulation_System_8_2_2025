#ifndef RUNACTION_HH
#define RUNACTION_HH

#include "G4UImanager.hh"
#include "G4RunManager.hh"
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

class UDRunAction: public G4UserRunAction {

	public:
		virtual void BeginOfRunAction(const G4Run *);
		virtual void EndOfRunAction(const G4Run *);

	private:
		UDInCodeControllingSurface ControllingSurface;	
};


#endif