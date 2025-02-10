
#ifndef EVENTACATION_HH
#define EVENTACATION_HH

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4RunManager.hh"
#include <cstdio>
 
class UDEventAction : public G4UserEventAction {
public:

  UDEventAction();
 ~UDEventAction() override;
 void BeginOfEventAction(const G4Event *) override;
 void EndOfEventAction(const G4Event *) override;

 private:

};

#endif
