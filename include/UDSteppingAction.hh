#ifndef UDSTEPPINGACTION_HH
#define UDSTEPPINGACTION_HH

#include "G4AnalysisManager.hh"
#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include "UDInCodeControllingSurface.hh"
# define M_PI 3.14159265358979323846
// Forward declaration to avoid circular dependency
class UDEventAction;

class UDSteppingAction : public G4UserSteppingAction {
public:
    UDSteppingAction();
    virtual ~UDSteppingAction();

    virtual void UserSteppingAction(const G4Step* step) override;

private:
    UDInCodeControllingSurface ControllingSurface;	
// Pointer to EventAction for storing angles
};

#endif  // UDSTEPPINGACTION_HH