#ifndef UDINCODECONTROLLINGSURFACE_HH
#define UDINCODECONTROLLINGSURFACE_HH

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <math.h>
#include "G4LogicalVolumeStore.hh"
#include "G4SDManager.hh"
#include "G4LogicalVolume.hh"
#include "G4String.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4String.hh"

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TH1.h"
#include "TKey.h"
#include "TList.h"
#include "TCanvas.h"
#include <TNtuple.h>

# define M_PI 3.14159265358979323846

class UDInCodeControllingSurface {

public:
    void GetHit(G4String target);
    void RunStart();
    void addstep(int tar);
    int returnstep(int tar);
    void Getdeleted(std::vector<G4String> &outtype, G4String filename2);
    void CreateHistogram(G4String filet);
    void Getincident(double &inc);
    void StoreNumberOfParticles(int entries);
    void GetNumberOfParticles(int &outries);
    void Getfilename(G4String target, G4String &arr);
    void Particle_Dealer(int x,G4String particle,int& o);
    void Targeted_Particle_list(int net_number_detectors, std::vector<G4String> &outlist);
    void getfile(G4String output);
    void printfile();
    void Particles_no_finder(G4String particlet,int& k);
    void Particles_name_finder(G4String &Particle, int no);
    int GetNumSensitiveDetectors(std::vector<G4String> &arr2);

private:
    static std::vector<G4String> storedParticletype;
    static std::vector<G4int> steplist;  
    static G4String generator_file;
    static G4String Ntuple_file_temp;
    static G4String Ntuple_file;
    static G4int NumberOfInputParticles;
    static G4String Directory;

};

#endif // UDINCODECONTROLLINGSURFACE_HH

