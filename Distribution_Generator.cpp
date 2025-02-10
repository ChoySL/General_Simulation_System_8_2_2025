#include <TFile.h>
#include <TNtuple.h>
#include <string>

int main()
{

    int sw;
    int x = 0;
    while (x == 0)
    {
        int e2, nm;
        printf("Energy and no. of particles:\n");
        scanf("%d %d", &e2, &nm);
        string input = "C:/Geant4/geant4-v11.3.0-install/share/Geant4/Simulations/General_Simulation_System/Data_Input/";
        input = input + std::to_string(e2) + "_GeV_distribution_muon.root";
        std::cout << input << std::endl;

        TFile file(input.c_str(), "RECREATE");
        TNtuple *ntuple = new TNtuple("ntuple", "E", "e:px:py:pz:dx:dy:dz");

        for (int i = 1; i <= nm; i++)
        {

            double px2 = 0.0;
            double py2 = 5;
            double pz2 = 0.0;

            double dx2 = 0.0;
            double dy2 = -1;
            double dz2 = 0.0;
            ntuple->Fill(e2, px2, py2, pz2, dx2, dy2, dz2);
        }

        ntuple->Write();
        file.Close();
        e2 = 0;
        nm = 0;
    }
}