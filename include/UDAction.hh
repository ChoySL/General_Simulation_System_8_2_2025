#ifndef ACTION_HH
#define ACTION_HH

#include "G4VUserActionInitialization.hh"
#include "UDParticleGenerator.hh"
#include "UDcreateNtuple.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "UDRunAction.hh"
#include "UDEventAction.hh"
#include "UDSteppingAction.hh"

class Action : public G4VUserActionInitialization {
	public:
		Action();
		~Action();

		virtual void Build() const;
};

#endif