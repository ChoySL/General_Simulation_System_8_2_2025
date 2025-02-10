#include "UDEventAction.hh"

UDEventAction::UDEventAction() : G4UserEventAction() {}

UDEventAction::~UDEventAction() {}

void UDEventAction::BeginOfEventAction(const G4Event *Event)
{
    // std::printf("Event begins now\n");
    // int me = Event->GetEventID();
    // std::printf("Checker Event begins now:%d\n",me);
}

void UDEventAction::EndOfEventAction(const G4Event *)
{
    // std::printf("Event ends now\n");
}