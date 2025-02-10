#include "UDphysicslist.hh"

PhysicsList::PhysicsList() {
    SetVerboseLevel(0);
    RegisterPhysics(new G4EmStandardPhysics()); 
    RegisterPhysics(new G4OpticalPhysics());  
    RegisterPhysics(new G4DecayPhysics());  
    RegisterPhysics(new G4RadioactiveDecayPhysics()); 
}

PhysicsList::~PhysicsList() {}