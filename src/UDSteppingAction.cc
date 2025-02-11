#include "UDSteppingAction.hh"
#include "UDEventAction.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include <cmath>

UDSteppingAction::UDSteppingAction()
    : G4UserSteppingAction() {}//, fEventAction(eventAction) 

    UDSteppingAction::~UDSteppingAction() {}

void UDSteppingAction::UserSteppingAction(const G4Step* step) {
    G4Track* track = step->GetTrack();
    G4StepPoint *preStepPoint = step->GetPreStepPoint();
	G4StepPoint *postStepPoint = step->GetPostStepPoint();
    G4int eventNum = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    G4ThreeVector PremomMuon = preStepPoint->GetMomentum();
    G4ThreeVector postmomMuon = postStepPoint->GetMomentum();
    G4ThreeVector posMuon = preStepPoint->GetPosition();
    G4double PremuonMomentum = PremomMuon.mag()/1000;
    G4double postmuonMomentum = postmomMuon.mag()/1000;  
    G4double muonMomentum = PremuonMomentum-postmuonMomentum;
    int stepN = ControllingSurface.returnstep(eventNum);
    ControllingSurface.addstep(eventNum);
    if (track->GetParentID() == 0) {  // Only primary particles
        G4AnalysisManager *manager_stepping = G4AnalysisManager::Instance();
        std::vector< G4String > Detector_names;
        int net = ControllingSurface.GetNumSensitiveDetectors(Detector_names);
        std::vector<G4String> particle_list;
        ControllingSurface.Targeted_Particle_list(net,particle_list);
        int num_particles = particle_list.size();
        int pps=num_particles+net+3;
        G4ThreeVector preposMuon = preStepPoint->GetPosition();
        G4ThreeVector preMomentum = step->GetPreStepPoint()->GetMomentumDirection();
        G4double hitTime = preStepPoint->GetGlobalTime();
        G4double netvpre=preMomentum[0]*preMomentum[0]+preMomentum[1]*preMomentum[1]+preMomentum[2]*preMomentum[2];
        netvpre=sqrt(netvpre);
        G4ThreeVector momUDpre;
        momUDpre[0]=preMomentum[0]/netvpre;
        momUDpre[1]=preMomentum[1]/netvpre;
        momUDpre[2]=preMomentum[2]/netvpre;
        double inclination = std::acos(momUDpre[2]);
        double azimuth = std::atan2(momUDpre[1], momUDpre[0]);
        G4ThreeVector postMomentum = step->GetPostStepPoint()->GetMomentumDirection();
        G4double netvpost=postMomentum[0]*postMomentum[0]+postMomentum[1]*postMomentum[1]+postMomentum[2]*postMomentum[2];
        netvpost=sqrt(netvpost);
        G4ThreeVector momUDpost;
        momUDpost[0]=postMomentum[0]/netvpost;
        momUDpost[1]=postMomentum[1]/netvpost;
        momUDpost[2]=postMomentum[2]/netvpost;
        double inclination1 = std::acos(momUDpost[2]);
        double azimuth1 = std::atan2(momUDpost[1], momUDpost[0]);
        double inclination_diff=inclination1-inclination;
        double azimuth_diff=azimuth1-azimuth;
        
        manager_stepping->GetNtuple("Primany Particle Stepping");
        manager_stepping->FillNtupleIColumn(pps, 0, eventNum);
        manager_stepping->FillNtupleIColumn(pps, 1, stepN);
        manager_stepping->FillNtupleDColumn(pps, 2, muonMomentum);
        manager_stepping->FillNtupleDColumn(pps, 3, hitTime);
        manager_stepping->FillNtupleDColumn(pps, 4, preposMuon[0]);
        manager_stepping->FillNtupleDColumn(pps, 5, preposMuon[1]);
        manager_stepping->FillNtupleDColumn(pps, 6, preposMuon[2]);
        manager_stepping->FillNtupleDColumn(pps, 7, inclination1);
        manager_stepping->FillNtupleDColumn(pps, 8, azimuth1);
        manager_stepping->FillNtupleDColumn(pps, 9, inclination_diff);
        manager_stepping->FillNtupleDColumn(pps, 10, azimuth_diff);
        manager_stepping->AddNtupleRow(pps);

        //fEventAction->AddScatteringAngle(scatteringAngle);
    }
}