#include "UDInCodeControllingSurface.hh"

#define M_PI 3.14159265358979323846
std::vector<G4String> UDInCodeControllingSurface::storedParticletype;
G4String UDInCodeControllingSurface::generator_file;
G4String UDInCodeControllingSurface::Ntuple_file;
G4String UDInCodeControllingSurface::Ntuple_file_temp;
G4int UDInCodeControllingSurface::NumberOfInputParticles;
G4String UDInCodeControllingSurface::Directory;
double here_angle = 90;
int here_energy = 1;
int here_particle_amount = 1;
std::string here_modifier_output = "_";

void UDInCodeControllingSurface::Getfilename(G4String target, G4String &output)
{

    if (target == "generator")
    {

        output = "C:/Geant4/geant4-v11.3.0-install/share/Geant4/Simulations/";
        output = output + "General_Simulation_System/";
        Directory = output;
        output = output + "Event_Data_Input/";
        output = output + here_energy + "_GeV_distribution_muon.root";
        generator_file = output;
    }

    if (target == "Ntuple_Output")
    {

        output = Directory;
        output = output + "Event_Data_Output/";
        output = output + "output_" + here_energy + here_modifier_output + here_particle_amount + ".root";
        Ntuple_file = output;
    }

    if (target == "Ntuple_temp")
    {

        output = "Temportary/";
        output = output + "output1.root";
        Ntuple_file_temp = output;
    }
}
void UDInCodeControllingSurface::Getincident(double &inc)
{
    inc = (here_angle * M_PI) / 180.0;
}

void UDInCodeControllingSurface::GetHit(G4String target)
{

    bool found = false;

    for (int u = 0; u < storedParticletype.size(); u++)
    {
        if (target == storedParticletype[u])
        {
            found = true;
            // G4cout << "Found : " << target << G4endl;
            break; // Exit the loop if the element is found
        }
    }

    if (!found)
    {
        storedParticletype.push_back(target);
    }
}
void UDInCodeControllingSurface::StoreNumberOfParticles(int entries)
{
    NumberOfInputParticles = entries;
}
void UDInCodeControllingSurface::GetNumberOfParticles(int &outries)
{
    outries = NumberOfInputParticles;
}
void UDInCodeControllingSurface::Getdeleted(std::vector<G4String> &outtype, G4String filename2)
{
    outtype = storedParticletype;
    int s = storedParticletype.size();
    G4String filetrans = "C:/Geant4/geant4-v11.3.0-install/share/Geant4/Simulations/General_Simulation_Setup_testing_ground/build/Release/";
    G4String filetransout = filetrans + "Data_Output/Matter_lening/With_all/";
    filetransout = filetransout + ".root";
    filetrans = filetrans + Ntuple_file_temp;
    TFile file(filetrans.c_str(), "READ");
    TFile output(filename2.c_str(), "RECREATE");
    if (!file.IsOpen() || file.IsZombie())
    {
        // std::cerr << "Error opening file." << std::endl;
        return;
    }

    TList *keyList = file.GetListOfKeys();
    if (!keyList)
    {
        // std::cerr << "Error: No keys found in the file." << std::endl;
        file.Close();
        return;
    }

    TIter next(keyList);
    TKey *key;
    while ((key = (TKey *)next()))
    {
        // std::cout << "----------------------------------------" << std::endl;
        TTree *tree = dynamic_cast<TTree *>(file.Get(key->GetName()));
        if (tree)
        {
            Long64_t numEntries = tree->GetEntries();
            if (numEntries == 0)
            {
                // std::cout << "The tree is empty: " << key->GetName() << std::endl;
            }
            else
            {
                // std::cout << "The tree is not empty. It has " << numEntries << " entries: " << key->GetName() << std::endl;
                TTree *newTree = tree->CloneTree();
                output.cd();
                newTree->Write();
            }
        }
        else
        {
            // std::cerr << "Error: TTree not found with key name: " << key->GetName() << std::endl;
        }
    }
    output.Close();
    file.Close();
}

void UDInCodeControllingSurface::CreateHistogram(G4String filet)
{
    G4bool particle = false;
    TFile file(filet.c_str(), "UPDATE");

    if (!file.IsOpen() || file.IsZombie())
    {
        // std::cerr << "Error opening file." << std::endl;
        return;
    }

    TList *keyList = file.GetListOfKeys();
    if (!keyList)
    {
        // std::cerr << "Error: No keys found in the file." << std::endl;
        file.Close();
        return;
    }
    TIter next(keyList);
    TKey *key;
    while ((key = (TKey *)next()))
    {
        std::string treename = key->GetName();
        TTree *tree = dynamic_cast<TTree *>(file.Get(treename.c_str()));
        if (tree)
        {
            Long64_t numEntries = tree->GetEntries();
            int nett = numEntries;
            if (numEntries == 0)
            {
                // std::cout << "The tree is empty: " << key->GetName() << std::endl;
            }
            else
            {

                for (int j = 0; j < storedParticletype.size(); j++)
                {
                    std::string word = storedParticletype[j];
                    if (treename.find(word) != std::string::npos)
                    {
                        particle = true;
                    }
                }
                if (particle)
                {
                    Double_t dataValueI, dataValueA, dataValueID;
                    ///////////////////////////////////////////////////////////////////
                    TBranch *branchI = tree->GetBranch("Inclination [0,Pi]");
                    std::string Iname = treename.c_str();
                    Iname = Iname + " Inclination";
                    TH1F *histogram1 = new TH1F(Iname.c_str(), Iname.c_str(), 100, 1, 2);
                    if (branchI)
                    {
                        branchI->SetAddress(&dataValueI);
                    }
                    std::vector<Double_t> neterI;
                    ///////////////////////////////////////////////////////////////////
                    TBranch *branchID = tree->GetBranch("Inclination_diff [0,Pi]");
                    std::string IDname = treename.c_str();
                    IDname = IDname + " Inclination_diff";
                    TH1F *histogramID = new TH1F(IDname.c_str(), IDname.c_str(), 100, -2, 2);
                    if (branchID)
                    {
                        branchID->SetAddress(&dataValueID);
                    }
                    std::vector<Double_t> neterID;
                    ///////////////////////////////////////////////////////////////////
                    TBranch *branchA = tree->GetBranch("Azimuth [-Pi,Pi)");
                    std::string Aname = treename.c_str();
                    Aname = Aname + " Azimuth";
                    TH1F *histogram1a = new TH1F(Aname.c_str(), Aname.c_str(), 100, -2, -1.1);
                    if (branchA)
                    {
                        branchA->SetAddress(&dataValueA);
                    }
                    std::vector<Double_t> neterA;
                    ///////////////////////////////////////////////////////////////////

                    ///////////////////////////////////////////////////////////////////
                    double netI = 0.0;
                    double netA = 0.0;
                    double meanI = 0.0;
                    double meanA = 0.0;
                    double varianceII = 0.0;
                    double varianceAA = 0.0;
                    double varianceI = 0.0;
                    double varianceA = 0.0;
                    double SDI = 0.0;
                    double SDA = 0.0;

                    for (int i = 0; i < numEntries; i++)
                    {
                        tree->GetEntry(i);
                        neterI.push_back(dataValueI);
                        neterID.push_back(dataValueID);
                        // G4cout << "Inclination" << dataValueI <<G4endl;
                        netI += dataValueI;
                        varianceII += dataValueI * dataValueI;
                        neterA.push_back(dataValueA);
                        netA += dataValueA;
                        varianceAA += dataValueA * dataValueA;
                    }
                    meanI = netI / numEntries;
                    meanA = netA / numEntries;
                    varianceI = (varianceII / numEntries) - (meanI * meanI);
                    varianceA = (varianceAA / numEntries) - (meanA * meanA);
                    SDI = std::sqrt(varianceI);
                    SDA = std::sqrt(varianceA);

                    for (Long64_t i = 0; i < numEntries; i++)
                    {
                        histogram1->Fill(neterI[i]);
                        histogram1a->Fill(neterA[i]);
                        histogramID->Fill(neterID[i]);
                    }
                    histogram1->Write();
                    histogramID->Write();
                    histogram1a->Write();
                }
            }
        }
    }
    file.Close();
    printf("Completed\n");
    printf("Data Stored in");
    printf(Ntuple_file);
    printf("\n");
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int UDInCodeControllingSurface::GetNumSensitiveDetectors(std::vector<G4String> &arr2)
{
    G4LogicalVolumeStore *lvStore = G4LogicalVolumeStore::GetInstance();
    G4SDManager *sdManager = G4SDManager::GetSDMpointer();
    int numSensitiveDetectors = 0;

    for (size_t i = 0; i < lvStore->size(); i++)
    {
        G4LogicalVolume *logicVolume = (*lvStore)[i];

        if (logicVolume->GetSensitiveDetector() != nullptr)
        {
            numSensitiveDetectors++;
            G4String name22 = logicVolume->GetName();
            arr2.push_back(name22);
        }
    }

    return numSensitiveDetectors;
}

void UDInCodeControllingSurface::Targeted_Particle_list(int net_number_detectors, std::vector<G4String> &outlist)
{
    std::vector<G4String> list = {"B+", "B-", "B0", "Bc+", "Bc-", "Bs0", "D+", "D-", "D0",
                                  "Ds+", "Ds-", "GenericIon", "He3", "J/psi", "N(1440)+", "N(1440)0", "N(1520)+", "N(1520)0",
                                  "N(1535)+", "N(1535)0", "N(1650)+", "N(1650)0", "N(1675)+", "N(1675)0", "N(1680)+", "N(1680)0",
                                  "N(1700)+", "N(1700)0", "N(1710)+", "N(1710)0", "N(1720)+", "N(1720)0", "N(1900)+", "N(1900)0",
                                  "N(1990)+", "N(1990)0", "N(2090)+", "N(2090)0", "N(2190)+", "N(2190)0", "N(2220)+", "N(2220)0",
                                  "N(2250)+", "N(2250)0", "Upsilon", "a0(1450)+", "a0(1450)-", "a0(1450)0", "a0(980)+", "a0(980)-",
                                  "a0(980)0", "a1(1260)+", "a1(1260)-", "a1(1260)0", "a2(1320)+", "a2(1320)-", "a2(1320)0", "alpha",
                                  "anti_B0", "anti_Bs0", "anti_D0", "anti_He3", "anti_N(1440)+", "anti_N(1440)0", "anti_N(1520)+",
                                  "anti_N(1520)0", "anti_N(1535)+", "anti_N(1535)0", "anti_N(1650)+", "anti_N(1650)0", "anti_N(1675)+",
                                  "anti_N(1675)0", "anti_N(1680)+", "anti_N(1680)0", "anti_N(1700)+", "anti_N(1700)0", "anti_N(1710)+",
                                  "anti_N(1710)0", "anti_N(1720)+", "anti_N(1720)0", "anti_N(1900)+", "anti_N(1900)0", "anti_N(1990)+",
                                  "anti_N(1990)0", "anti_N(2090)+", "anti_N(2090)0", "anti_N(2190)+", "anti_N(2190)0", "anti_N(2220)+",
                                  "anti_N(2220)0", "anti_N(2250)+", "anti_N(2250)0", "anti_alpha", "anti_b_quark", "anti_bb1_diquark",
                                  "anti_bc0_diquark", "anti_bc1_diquark", "anti_bd0_diquark", "anti_bd1_diquark", "anti_bs0_diquark",
                                  "anti_bs1_diquark", "anti_bu0_diquark", "anti_bu1_diquark", "anti_c_quark", "anti_cc1_diquark",
                                  "anti_cd0_diquark", "anti_cd1_diquark", "anti_cs0_diquark", "anti_cs1_diquark", "anti_cu0_diquark",
                                  "anti_cu1_diquark", "anti_d_quark", "anti_dd1_diquark", "anti_delta(1600)+", "anti_delta(1600)++",
                                  "anti_delta(1600)-", "anti_delta(1600)0", "anti_delta(1620)+", "anti_delta(1620)++",
                                  "anti_delta(1620)-", "anti_delta(1620)0", "anti_delta(1700)+", "anti_delta(1700)++",
                                  "anti_delta(1700)-", "anti_delta(1700)0", "anti_delta(1900)+", "anti_delta(1900)++",
                                  "anti_delta(1900)-", "anti_delta(1900)0", "anti_delta(1905)+", "anti_delta(1905)++",
                                  "anti_delta(1905)-", "anti_delta(1905)0", "anti_delta(1910)+", "anti_delta(1910)++",
                                  "anti_delta(1910)-", "anti_delta(1910)0", "anti_delta(1920)+", "anti_delta(1920)++",
                                  "anti_delta(1920)-", "anti_delta(1920)0", "anti_delta(1930)+", "anti_delta(1930)++",
                                  "anti_delta(1930)-", "anti_delta(1930)0", "anti_delta(1950)+", "anti_delta(1950)++",
                                  "anti_delta(1950)-", "anti_delta(1950)0", "anti_delta+", "anti_delta++", "anti_delta-",
                                  "anti_delta0", "anti_deuteron", "anti_doublehyperH4", "anti_doublehyperdoubleneutron",
                                  "anti_hyperH4", "anti_hyperHe5", "anti_hyperalpha", "anti_hypertriton", "anti_k(1460)0",
                                  "anti_k0_star(1430)0", "anti_k1(1270)0", "anti_k1(1400)0", "anti_k2(1770)0", "anti_k2_star(1430)0",
                                  "anti_k2_star(1980)0", "anti_k3_star(1780)0", "anti_k_star(1410)0", "anti_k_star(1680)0",
                                  "anti_k_star0", "anti_kaon0", "anti_lambda", "anti_lambda(1405)", "anti_lambda(1520)",
                                  "anti_lambda(1600)", "anti_lambda(1670)", "anti_lambda(1690)", "anti_lambda(1800)",
                                  "anti_lambda(1810)", "anti_lambda(1820)", "anti_lambda(1830)", "anti_lambda(1890)",
                                  "anti_lambda(2100)", "anti_lambda(2110)", "anti_lambda_b", "anti_lambda_c+", "anti_neutron",
                                  "anti_nu_e", "anti_nu_mu", "anti_nu_tau", "anti_omega-", "anti_omega_b-", "anti_omega_c0",
                                  "anti_proton", "anti_s_quark", "anti_sd0_diquark", "anti_sd1_diquark", "anti_sigma(1385)+",
                                  "anti_sigma(1385)-", "anti_sigma(1385)0", "anti_sigma(1660)+", "anti_sigma(1660)-", "anti_sigma(1660)0",
                                  "anti_sigma(1670)+", "anti_sigma(1670)-", "anti_sigma(1670)0", "anti_sigma(1750)+", "anti_sigma(1750)-",
                                  "anti_sigma(1750)0", "anti_sigma(1775)+", "anti_sigma(1775)-", "anti_sigma(1775)0", "anti_sigma(1915)+",
                                  "anti_sigma(1915)-", "anti_sigma(1915)0", "anti_sigma(1940)+", "anti_sigma(1940)-", "anti_sigma(1940)0",
                                  "anti_sigma(2030)+", "anti_sigma(2030)-", "anti_sigma(2030)0", "anti_sigma+", "anti_sigma-", "anti_sigma0",
                                  "anti_sigma_b+", "anti_sigma_b-", "anti_sigma_b0", "anti_sigma_c+", "anti_sigma_c++", "anti_sigma_c0",
                                  "anti_ss1_diquark", "anti_su0_diquark", "anti_su1_diquark", "anti_t_quark", "anti_triton", "anti_u_quark",
                                  "anti_ud0_diquark", "anti_ud1_diquark", "anti_uu1_diquark", "anti_xi(1530)-", "anti_xi(1530)0",
                                  "anti_xi(1690)-", "anti_xi(1690)0", "anti_xi(1820)-", "anti_xi(1820)0", "anti_xi(1950)-",
                                  "anti_xi(1950)0", "anti_xi(2030)-", "anti_xi(2030)0", "anti_xi-", "anti_xi0", "anti_xi_b-", "anti_xi_b0",
                                  "anti_xi_c+", "anti_xi_c0", "b1(1235)+", "b1(1235)-", "b1(1235)0", "b_quark", "bb1_diquark", "bc0_diquark",
                                  "bc1_diquark", "bd0_diquark", "bd1_diquark", "bs0_diquark", "bs1_diquark", "bu0_diquark", "bu1_diquark",
                                  "c_quark", "cc1_diquark", "cd0_diquark", "cd1_diquark", "chargedgeantino", "cs0_diquark", "cs1_diquark",
                                  "cu0_diquark", "cu1_diquark", "d_quark", "dd1_diquark", "delta(1600)+", "delta(1600)++", "delta(1600)-",
                                  "delta(1600)0", "delta(1620)+", "delta(1620)++", "delta(1620)-", "delta(1620)0", "delta(1700)+",
                                  "delta(1700)++", "delta(1700)-", "delta(1700)0", "delta(1900)+", "delta(1900)++", "delta(1900)-",
                                  "delta(1900)0", "delta(1905)+", "delta(1905)++", "delta(1905)-", "delta(1905)0", "delta(1910)+",
                                  "delta(1910)++", "delta(1910)-", "delta(1910)0", "delta(1920)+", "delta(1920)++", "delta(1920)-",
                                  "delta(1920)0", "delta(1930)+", "delta(1930)++", "delta(1930)-", "delta(1930)0", "delta(1950)+",
                                  "delta(1950)++", "delta(1950)-", "delta(1950)0", "delta+", "delta++", "delta-", "delta0", "deuteron",
                                  "doublehyperH4", "doublehyperdoubleneutron", "e+", "e-", "eta", "eta(1295)", "eta(1405)", "eta(1475)",
                                  "eta2(1645)", "eta2(1870)", "eta_prime", "etac", "f0(1370)", "f0(1500)", "f0(1710)", "f0(500)", "f0(980)",
                                  "f1(1285)", "f1(1420)", "f2(1270)", "f2(1810)", "f2(2010)", "f2_prime(1525)", "gamma", "geantino", "gluon",
                                  "h1(1170)", "h1(1380)", "hyperH4", "hyperHe5", "hyperalpha", "hypertriton", "k(1460)+", "k(1460)-",
                                  "k(1460)0", "k0_star(1430)+", "k0_star(1430)-", "k0_star(1430)0", "k1(1270)+", "k1(1270)-", "k1(1270)0",
                                  "k1(1400)+", "k1(1400)-", "k1(1400)0", "k2(1770)+", "k2(1770)-", "k2(1770)0", "k2_star(1430)+",
                                  "k2_star(1430)-", "k2_star(1430)0", "k2_star(1980)+", "k2_star(1980)-", "k2_star(1980)0",
                                  "k3_star(1780)+", "k3_star(1780)-", "k3_star(1780)0", "k_star(1410)+", "k_star(1410)-", "k_star(1410)0",
                                  "k_star(1680)+", "k_star(1680)-", "k_star(1680)0", "k_star+", "k_star-", "k_star0", "kaon+", "kaon-", "kaon0",
                                  "kaon0L", "kaon0S", "lambda", "lambda(1405)", "lambda(1520)", "lambda(1600)", "lambda(1670)", "lambda(1690)",
                                  "lambda(1800)", "lambda(1810)", "lambda(1820)", "lambda(1830)", "lambda(1890)", "lambda(2100)",
                                  "lambda(2110)", "lambda_b", "lambda_c+", "mu+", "mu-", "neutron", "nu_e", "nu_mu", "nu_tau", "omega",
                                  "omega(1420)", "omega(1650)", "omega-", "omega3(1670)", "omega_b-", "omega_c0", "opticalphoton", "phi",
                                  "phi(1680)", "phi3(1850)", "pi(1300)+", "pi(1300)-", "pi(1300)0", "pi+", "pi-", "pi0", "pi2(1670)+",
                                  "pi2(1670)-", "pi2(1670)0", "proton", "rho(1450)+", "rho(1450)-", "rho(1450)0", "rho(1700)+", "rho(1700)-",
                                  "rho(1700)0", "rho+", "rho-", "rho0", "rho3(1690)+", "rho3(1690)-", "rho3(1690)0", "s_quark", "sd0_diquark",
                                  "sd1_diquark", "sigma(1385)+", "sigma(1385)-", "sigma(1385)0", "sigma(1660)+", "sigma(1660)-",
                                  "sigma(1660)0", "sigma(1670)+", "sigma(1670)-", "sigma(1670)0", "sigma(1750)+", "sigma(1750)-", "sigma(1750)0",
                                  "sigma(1775)+", "sigma(1775)-", "sigma(1775)0", "sigma(1915)+", "sigma(1915)-", "sigma(1915)0", "sigma(1940)+",
                                  "sigma(1940)-", "sigma(1940)0", "sigma(2030)+", "sigma(2030)-", "sigma(2030)0", "sigma+", "sigma-", "sigma0",
                                  "sigma_b+", "sigma_b-", "sigma_b0", "sigma_c+", "sigma_c++", "sigma_c0", "ss1_diquark", "su0_diquark", "su1_diquark",
                                  "t_quark", "tau+", "tau-", "triton", "u_quark", "ud0_diquark", "ud1_diquark", "uu1_diquark", "xi(1530)-",
                                  "xi(1530)0", "xi(1690)-", "xi(1690)0", "xi(1820)-", "xi(1820)0", "xi(1950)-", "xi(1950)0", "xi(2030)-",
                                  "xi(2030)0", "xi-", "xi0", "xi_b-", "xi_b0", "xi_c+", "xi_c0"}; //,"opticalphoton"

    outlist = list;
}

void UDInCodeControllingSurface::Particle_Dealer(int net_number_detectors, G4String particle, int &output)
{
    int j;
    if (particle == "B+")
    {
        j = 0;
    }
    if (particle == "B-")
    {
        j = 1;
    }
    if (particle == "B0")
    {
        j = 2;
    }
    if (particle == "Bc+")
    {
        j = 3;
    }
    if (particle == "Bc-")
    {
        j = 4;
    }
    if (particle == "Bs0")
    {
        j = 5;
    }
    if (particle == "D+")
    {
        j = 6;
    }
    if (particle == "D-")
    {
        j = 7;
    }
    if (particle == "D0")
    {
        j = 8;
    }
    if (particle == "Ds+")
    {
        j = 9;
    }
    if (particle == "Ds-")
    {
        j = 10;
    }
    if (particle == "GenericIon")
    {
        j = 11;
    }
    if (particle == "He3")
    {
        j = 12;
    }
    if (particle == "J/psi")
    {
        j = 13;
    }
    if (particle == "N(1440)+")
    {
        j = 14;
    }
    if (particle == "N(1440)0")
    {
        j = 15;
    }
    if (particle == "N(1520)+")
    {
        j = 16;
    }
    if (particle == "N(1520)0")
    {
        j = 17;
    }
    if (particle == "N(1535)+")
    {
        j = 18;
    }
    if (particle == "N(1535)0")
    {
        j = 19;
    }
    if (particle == "N(1650)+")
    {
        j = 20;
    }
    if (particle == "N(1650)0")
    {
        j = 21;
    }
    if (particle == "N(1675)+")
    {
        j = 22;
    }
    if (particle == "N(1675)0")
    {
        j = 23;
    }
    if (particle == "N(1680)+")
    {
        j = 24;
    }
    if (particle == "N(1680)0")
    {
        j = 25;
    }
    if (particle == "N(1700)+")
    {
        j = 26;
    }
    if (particle == "N(1700)0")
    {
        j = 27;
    }
    if (particle == "N(1710)+")
    {
        j = 28;
    }
    if (particle == "N(1710)0")
    {
        j = 29;
    }
    if (particle == "N(1720)+")
    {
        j = 30;
    }
    if (particle == "N(1720)0")
    {
        j = 31;
    }
    if (particle == "N(1900)+")
    {
        j = 32;
    }
    if (particle == "N(1900)0")
    {
        j = 33;
    }
    if (particle == "N(1990)+")
    {
        j = 34;
    }
    if (particle == "N(1990)0")
    {
        j = 35;
    }
    if (particle == "N(2090)+")
    {
        j = 36;
    }
    if (particle == "N(2090)0")
    {
        j = 37;
    }
    if (particle == "N(2190)+")
    {
        j = 38;
    }
    if (particle == "N(2190)0")
    {
        j = 39;
    }
    if (particle == "N(2220)+")
    {
        j = 40;
    }
    if (particle == "N(2220)0")
    {
        j = 41;
    }
    if (particle == "N(2250)+")
    {
        j = 42;
    }
    if (particle == "N(2250)0")
    {
        j = 43;
    }
    if (particle == "Upsilon")
    {
        j = 44;
    }
    if (particle == "a0(1450)+")
    {
        j = 45;
    }
    if (particle == "a0(1450)-")
    {
        j = 46;
    }
    if (particle == "a0(1450)0")
    {
        j = 47;
    }
    if (particle == "a0(980)+")
    {
        j = 48;
    }
    if (particle == "a0(980)-")
    {
        j = 49;
    }
    if (particle == "a0(980)0")
    {
        j = 50;
    }
    if (particle == "a1(1260)+")
    {
        j = 51;
    }
    if (particle == "a1(1260)-")
    {
        j = 52;
    }
    if (particle == "a1(1260)0")
    {
        j = 53;
    }
    if (particle == "a2(1320)+")
    {
        j = 54;
    }
    if (particle == "a2(1320)-")
    {
        j = 55;
    }
    if (particle == "a2(1320)0")
    {
        j = 56;
    }
    if (particle == "alpha")
    {
        j = 57;
    }
    if (particle == "anti_B0")
    {
        j = 58;
    }
    if (particle == "anti_Bs0")
    {
        j = 59;
    }
    if (particle == "anti_D0")
    {
        j = 60;
    }
    if (particle == "anti_He3")
    {
        j = 61;
    }
    if (particle == "anti_N(1440)+")
    {
        j = 62;
    }
    if (particle == "anti_N(1440)0")
    {
        j = 63;
    }
    if (particle == "anti_N(1520)+")
    {
        j = 64;
    }
    if (particle == "anti_N(1520)0")
    {
        j = 65;
    }
    if (particle == "anti_N(1535)+")
    {
        j = 66;
    }
    if (particle == "anti_N(1535)0")
    {
        j = 67;
    }
    if (particle == "anti_N(1650)+")
    {
        j = 68;
    }
    if (particle == "anti_N(1650)0")
    {
        j = 69;
    }
    if (particle == "anti_N(1675)+")
    {
        j = 70;
    }
    if (particle == "anti_N(1675)0")
    {
        j = 71;
    }
    if (particle == "anti_N(1680)+")
    {
        j = 72;
    }
    if (particle == "anti_N(1680)0")
    {
        j = 73;
    }
    if (particle == "anti_N(1700)+")
    {
        j = 74;
    }
    if (particle == "anti_N(1700)0")
    {
        j = 75;
    }
    if (particle == "anti_N(1710)+")
    {
        j = 76;
    }
    if (particle == "anti_N(1710)0")
    {
        j = 77;
    }
    if (particle == "anti_N(1720)+")
    {
        j = 78;
    }
    if (particle == "anti_N(1720)0")
    {
        j = 79;
    }
    if (particle == "anti_N(1900)+")
    {
        j = 80;
    }
    if (particle == "anti_N(1900)0")
    {
        j = 81;
    }
    if (particle == "anti_N(1990)+")
    {
        j = 82;
    }
    if (particle == "anti_N(1990)0")
    {
        j = 83;
    }
    if (particle == "anti_N(2090)+")
    {
        j = 84;
    }
    if (particle == "anti_N(2090)0")
    {
        j = 85;
    }
    if (particle == "anti_N(2190)+")
    {
        j = 86;
    }
    if (particle == "anti_N(2190)0")
    {
        j = 87;
    }
    if (particle == "anti_N(2220)+")
    {
        j = 88;
    }
    if (particle == "anti_N(2220)0")
    {
        j = 89;
    }
    if (particle == "anti_N(2250)+")
    {
        j = 90;
    }
    if (particle == "anti_N(2250)0")
    {
        j = 91;
    }
    if (particle == "anti_alpha")
    {
        j = 92;
    }
    if (particle == "anti_b_quark")
    {
        j = 93;
    }
    if (particle == "anti_bb1_diquark")
    {
        j = 94;
    }
    if (particle == "anti_bc0_diquark")
    {
        j = 95;
    }
    if (particle == "anti_bc1_diquark")
    {
        j = 96;
    }
    if (particle == "anti_bd0_diquark")
    {
        j = 97;
    }
    if (particle == "anti_bd1_diquark")
    {
        j = 98;
    }
    if (particle == "anti_bs0_diquark")
    {
        j = 99;
    }
    if (particle == "anti_bs1_diquark")
    {
        j = 100;
    }
    if (particle == "anti_bu0_diquark")
    {
        j = 101;
    }
    if (particle == "anti_bu1_diquark")
    {
        j = 102;
    }
    if (particle == "anti_c_quark")
    {
        j = 103;
    }
    if (particle == "anti_cc1_diquark")
    {
        j = 104;
    }
    if (particle == "anti_cd0_diquark")
    {
        j = 105;
    }
    if (particle == "anti_cd1_diquark")
    {
        j = 106;
    }
    if (particle == "anti_cs0_diquark")
    {
        j = 107;
    }
    if (particle == "anti_cs1_diquark")
    {
        j = 108;
    }
    if (particle == "anti_cu0_diquark")
    {
        j = 109;
    }
    if (particle == "anti_cu1_diquark")
    {
        j = 110;
    }
    if (particle == "anti_d_quark")
    {
        j = 111;
    }
    if (particle == "anti_dd1_diquark")
    {
        j = 112;
    }
    if (particle == "anti_delta(1600)+")
    {
        j = 113;
    }
    if (particle == "anti_delta(1600)++")
    {
        j = 114;
    }
    if (particle == "anti_delta(1600)-")
    {
        j = 115;
    }
    if (particle == "anti_delta(1600)0")
    {
        j = 116;
    }
    if (particle == "anti_delta(1620)+")
    {
        j = 117;
    }
    if (particle == "anti_delta(1620)++")
    {
        j = 118;
    }
    if (particle == "anti_delta(1620)-")
    {
        j = 119;
    }
    if (particle == "anti_delta(1620)0")
    {
        j = 120;
    }
    if (particle == "anti_delta(1700)+")
    {
        j = 121;
    }
    if (particle == "anti_delta(1700)++")
    {
        j = 122;
    }
    if (particle == "anti_delta(1700)-")
    {
        j = 123;
    }
    if (particle == "anti_delta(1700)0")
    {
        j = 124;
    }
    if (particle == "anti_delta(1900)+")
    {
        j = 125;
    }
    if (particle == "anti_delta(1900)++")
    {
        j = 126;
    }
    if (particle == "anti_delta(1900)-")
    {
        j = 127;
    }
    if (particle == "anti_delta(1900)0")
    {
        j = 128;
    }
    if (particle == "anti_delta(1905)+")
    {
        j = 129;
    }
    if (particle == "anti_delta(1905)++")
    {
        j = 130;
    }
    if (particle == "anti_delta(1905)-")
    {
        j = 131;
    }
    if (particle == "anti_delta(1905)0")
    {
        j = 132;
    }
    if (particle == "anti_delta(1910)+")
    {
        j = 133;
    }
    if (particle == "anti_delta(1910)++")
    {
        j = 134;
    }
    if (particle == "anti_delta(1910)-")
    {
        j = 135;
    }
    if (particle == "anti_delta(1910)0")
    {
        j = 136;
    }
    if (particle == "anti_delta(1920)+")
    {
        j = 137;
    }
    if (particle == "anti_delta(1920)++")
    {
        j = 138;
    }
    if (particle == "anti_delta(1920)-")
    {
        j = 139;
    }
    if (particle == "anti_delta(1920)0")
    {
        j = 140;
    }
    if (particle == "anti_delta(1930)+")
    {
        j = 141;
    }
    if (particle == "anti_delta(1930)++")
    {
        j = 142;
    }
    if (particle == "anti_delta(1930)-")
    {
        j = 143;
    }
    if (particle == "anti_delta(1930)0")
    {
        j = 144;
    }
    if (particle == "anti_delta(1950)+")
    {
        j = 145;
    }
    if (particle == "anti_delta(1950)++")
    {
        j = 146;
    }
    if (particle == "anti_delta(1950)-")
    {
        j = 147;
    }
    if (particle == "anti_delta(1950)0")
    {
        j = 148;
    }
    if (particle == "anti_delta+")
    {
        j = 149;
    }
    if (particle == "anti_delta++")
    {
        j = 150;
    }
    if (particle == "anti_delta-")
    {
        j = 151;
    }
    if (particle == "anti_delta0")
    {
        j = 152;
    }
    if (particle == "anti_deuteron")
    {
        j = 153;
    }
    if (particle == "anti_doublehyperH4")
    {
        j = 154;
    }
    if (particle == "anti_doublehyperdoubleneutron")
    {
        j = 155;
    }
    if (particle == "anti_hyperH4")
    {
        j = 156;
    }
    if (particle == "anti_hyperHe5")
    {
        j = 157;
    }
    if (particle == "anti_hyperalpha")
    {
        j = 158;
    }
    if (particle == "anti_hypertriton")
    {
        j = 159;
    }
    if (particle == "anti_k(1460)0")
    {
        j = 160;
    }
    if (particle == "anti_k0_star(1430)0")
    {
        j = 161;
    }
    if (particle == "anti_k1(1270)0")
    {
        j = 162;
    }
    if (particle == "anti_k1(1400)0")
    {
        j = 163;
    }
    if (particle == "anti_k2(1770)0")
    {
        j = 164;
    }
    if (particle == "anti_k2_star(1430)0")
    {
        j = 165;
    }
    if (particle == "anti_k2_star(1980)0")
    {
        j = 166;
    }
    if (particle == "anti_k3_star(1780)0")
    {
        j = 167;
    }
    if (particle == "anti_k_star(1410)0")
    {
        j = 168;
    }
    if (particle == "anti_k_star(1680)0")
    {
        j = 169;
    }
    if (particle == "anti_k_star0")
    {
        j = 170;
    }
    if (particle == "anti_kaon0")
    {
        j = 171;
    }
    if (particle == "anti_lambda")
    {
        j = 172;
    }
    if (particle == "anti_lambda(1405)")
    {
        j = 173;
    }
    if (particle == "anti_lambda(1520)")
    {
        j = 174;
    }
    if (particle == "anti_lambda(1600)")
    {
        j = 175;
    }
    if (particle == "anti_lambda(1670)")
    {
        j = 176;
    }
    if (particle == "anti_lambda(1690)")
    {
        j = 177;
    }
    if (particle == "anti_lambda(1800)")
    {
        j = 178;
    }
    if (particle == "anti_lambda(1810)")
    {
        j = 179;
    }
    if (particle == "anti_lambda(1820)")
    {
        j = 180;
    }
    if (particle == "anti_lambda(1830)")
    {
        j = 181;
    }
    if (particle == "anti_lambda(1890)")
    {
        j = 182;
    }
    if (particle == "anti_lambda(2100)")
    {
        j = 183;
    }
    if (particle == "anti_lambda(2110)")
    {
        j = 184;
    }
    if (particle == "anti_lambda_b")
    {
        j = 185;
    }
    if (particle == "anti_lambda_c+")
    {
        j = 186;
    }
    if (particle == "anti_neutron")
    {
        j = 187;
    }
    if (particle == "anti_nu_e")
    {
        j = 188;
    }
    if (particle == "anti_nu_mu")
    {
        j = 189;
    }
    if (particle == "anti_nu_tau")
    {
        j = 190;
    }
    if (particle == "anti_omega-")
    {
        j = 191;
    }
    if (particle == "anti_omega_b-")
    {
        j = 192;
    }
    if (particle == "anti_omega_c0")
    {
        j = 193;
    }
    if (particle == "anti_proton")
    {
        j = 194;
    }
    if (particle == "anti_s_quark")
    {
        j = 195;
    }
    if (particle == "anti_sd0_diquark")
    {
        j = 196;
    }
    if (particle == "anti_sd1_diquark")
    {
        j = 197;
    }
    if (particle == "anti_sigma(1385)+")
    {
        j = 198;
    }
    if (particle == "anti_sigma(1385)-")
    {
        j = 199;
    }
    if (particle == "anti_sigma(1385)0")
    {
        j = 200;
    }
    if (particle == "anti_sigma(1660)+")
    {
        j = 201;
    }
    if (particle == "anti_sigma(1660)-")
    {
        j = 202;
    }
    if (particle == "anti_sigma(1660)0")
    {
        j = 203;
    }
    if (particle == "anti_sigma(1670)+")
    {
        j = 204;
    }
    if (particle == "anti_sigma(1670)-")
    {
        j = 205;
    }
    if (particle == "anti_sigma(1670)0")
    {
        j = 206;
    }
    if (particle == "anti_sigma(1750)+")
    {
        j = 207;
    }
    if (particle == "anti_sigma(1750)-")
    {
        j = 208;
    }
    if (particle == "anti_sigma(1750)0")
    {
        j = 209;
    }
    if (particle == "anti_sigma(1775)+")
    {
        j = 210;
    }
    if (particle == "anti_sigma(1775)-")
    {
        j = 211;
    }
    if (particle == "anti_sigma(1775)0")
    {
        j = 212;
    }
    if (particle == "anti_sigma(1915)+")
    {
        j = 213;
    }
    if (particle == "anti_sigma(1915)-")
    {
        j = 214;
    }
    if (particle == "anti_sigma(1915)0")
    {
        j = 215;
    }
    if (particle == "anti_sigma(1940)+")
    {
        j = 216;
    }
    if (particle == "anti_sigma(1940)-")
    {
        j = 217;
    }
    if (particle == "anti_sigma(1940)0")
    {
        j = 218;
    }
    if (particle == "anti_sigma(2030)+")
    {
        j = 219;
    }
    if (particle == "anti_sigma(2030)-")
    {
        j = 220;
    }
    if (particle == "anti_sigma(2030)0")
    {
        j = 221;
    }
    if (particle == "anti_sigma+")
    {
        j = 222;
    }
    if (particle == "anti_sigma-")
    {
        j = 223;
    }
    if (particle == "anti_sigma0")
    {
        j = 224;
    }
    if (particle == "anti_sigma_b+")
    {
        j = 225;
    }
    if (particle == "anti_sigma_b-")
    {
        j = 226;
    }
    if (particle == "anti_sigma_b0")
    {
        j = 227;
    }
    if (particle == "anti_sigma_c+")
    {
        j = 228;
    }
    if (particle == "anti_sigma_c++")
    {
        j = 229;
    }
    if (particle == "anti_sigma_c0")
    {
        j = 230;
    }
    if (particle == "anti_ss1_diquark")
    {
        j = 231;
    }
    if (particle == "anti_su0_diquark")
    {
        j = 232;
    }
    if (particle == "anti_su1_diquark")
    {
        j = 233;
    }
    if (particle == "anti_t_quark")
    {
        j = 234;
    }
    if (particle == "anti_triton")
    {
        j = 235;
    }
    if (particle == "anti_u_quark")
    {
        j = 236;
    }
    if (particle == "anti_ud0_diquark")
    {
        j = 237;
    }
    if (particle == "anti_ud1_diquark")
    {
        j = 238;
    }
    if (particle == "anti_uu1_diquark")
    {
        j = 239;
    }
    if (particle == "anti_xi(1530)-")
    {
        j = 240;
    }
    if (particle == "anti_xi(1530)0")
    {
        j = 241;
    }
    if (particle == "anti_xi(1690)-")
    {
        j = 242;
    }
    if (particle == "anti_xi(1690)0")
    {
        j = 243;
    }
    if (particle == "anti_xi(1820)-")
    {
        j = 244;
    }
    if (particle == "anti_xi(1820)0")
    {
        j = 245;
    }
    if (particle == "anti_xi(1950)-")
    {
        j = 246;
    }
    if (particle == "anti_xi(1950)0")
    {
        j = 247;
    }
    if (particle == "anti_xi(2030)-")
    {
        j = 248;
    }
    if (particle == "anti_xi(2030)0")
    {
        j = 249;
    }
    if (particle == "anti_xi-")
    {
        j = 250;
    }
    if (particle == "anti_xi0")
    {
        j = 251;
    }
    if (particle == "anti_xi_b-")
    {
        j = 252;
    }
    if (particle == "anti_xi_b0")
    {
        j = 253;
    }
    if (particle == "anti_xi_c+")
    {
        j = 254;
    }
    if (particle == "anti_xi_c0")
    {
        j = 255;
    }
    if (particle == "b1(1235)+")
    {
        j = 256;
    }
    if (particle == "b1(1235)-")
    {
        j = 257;
    }
    if (particle == "b1(1235)0")
    {
        j = 258;
    }
    if (particle == "b_quark")
    {
        j = 259;
    }
    if (particle == "bb1_diquark")
    {
        j = 260;
    }
    if (particle == "bc0_diquark")
    {
        j = 261;
    }
    if (particle == "bc1_diquark")
    {
        j = 262;
    }
    if (particle == "bd0_diquark")
    {
        j = 263;
    }
    if (particle == "bd1_diquark")
    {
        j = 264;
    }
    if (particle == "bs0_diquark")
    {
        j = 265;
    }
    if (particle == "bs1_diquark")
    {
        j = 266;
    }
    if (particle == "bu0_diquark")
    {
        j = 267;
    }
    if (particle == "bu1_diquark")
    {
        j = 268;
    }
    if (particle == "c_quark")
    {
        j = 269;
    }
    if (particle == "cc1_diquark")
    {
        j = 270;
    }
    if (particle == "cd0_diquark")
    {
        j = 271;
    }
    if (particle == "cd1_diquark")
    {
        j = 272;
    }
    if (particle == "chargedgeantino")
    {
        j = 273;
    }
    if (particle == "cs0_diquark")
    {
        j = 274;
    }
    if (particle == "cs1_diquark")
    {
        j = 275;
    }
    if (particle == "cu0_diquark")
    {
        j = 276;
    }
    if (particle == "cu1_diquark")
    {
        j = 277;
    }
    if (particle == "d_quark")
    {
        j = 278;
    }
    if (particle == "dd1_diquark")
    {
        j = 279;
    }
    if (particle == "delta(1600)+")
    {
        j = 280;
    }
    if (particle == "delta(1600)++")
    {
        j = 281;
    }
    if (particle == "delta(1600)-")
    {
        j = 282;
    }
    if (particle == "delta(1600)0")
    {
        j = 283;
    }
    if (particle == "delta(1620)+")
    {
        j = 284;
    }
    if (particle == "delta(1620)++")
    {
        j = 285;
    }
    if (particle == "delta(1620)-")
    {
        j = 286;
    }
    if (particle == "delta(1620)0")
    {
        j = 287;
    }
    if (particle == "delta(1700)+")
    {
        j = 288;
    }
    if (particle == "delta(1700)++")
    {
        j = 289;
    }
    if (particle == "delta(1700)-")
    {
        j = 290;
    }
    if (particle == "delta(1700)0")
    {
        j = 291;
    }
    if (particle == "delta(1900)+")
    {
        j = 292;
    }
    if (particle == "delta(1900)++")
    {
        j = 293;
    }
    if (particle == "delta(1900)-")
    {
        j = 294;
    }
    if (particle == "delta(1900)0")
    {
        j = 295;
    }
    if (particle == "delta(1905)+")
    {
        j = 296;
    }
    if (particle == "delta(1905)++")
    {
        j = 297;
    }
    if (particle == "delta(1905)-")
    {
        j = 298;
    }
    if (particle == "delta(1905)0")
    {
        j = 299;
    }
    if (particle == "delta(1910)+")
    {
        j = 300;
    }
    if (particle == "delta(1910)++")
    {
        j = 301;
    }
    if (particle == "delta(1910)-")
    {
        j = 302;
    }
    if (particle == "delta(1910)0")
    {
        j = 303;
    }
    if (particle == "delta(1920)+")
    {
        j = 304;
    }
    if (particle == "delta(1920)++")
    {
        j = 305;
    }
    if (particle == "delta(1920)-")
    {
        j = 306;
    }
    if (particle == "delta(1920)0")
    {
        j = 307;
    }
    if (particle == "delta(1930)+")
    {
        j = 308;
    }
    if (particle == "delta(1930)++")
    {
        j = 309;
    }
    if (particle == "delta(1930)-")
    {
        j = 310;
    }
    if (particle == "delta(1930)0")
    {
        j = 311;
    }
    if (particle == "delta(1950)+")
    {
        j = 312;
    }
    if (particle == "delta(1950)++")
    {
        j = 313;
    }
    if (particle == "delta(1950)-")
    {
        j = 314;
    }
    if (particle == "delta(1950)0")
    {
        j = 315;
    }
    if (particle == "delta+")
    {
        j = 316;
    }
    if (particle == "delta++")
    {
        j = 317;
    }
    if (particle == "delta-")
    {
        j = 318;
    }
    if (particle == "delta0")
    {
        j = 319;
    }
    if (particle == "deuteron")
    {
        j = 320;
    }
    if (particle == "doublehyperH4")
    {
        j = 321;
    }
    if (particle == "doublehyperdoubleneutron")
    {
        j = 322;
    }
    if (particle == "e+")
    {
        j = 323;
    }
    if (particle == "e-")
    {
        j = 324;
    }
    if (particle == "eta")
    {
        j = 325;
    }
    if (particle == "eta(1295)")
    {
        j = 326;
    }
    if (particle == "eta(1405)")
    {
        j = 327;
    }
    if (particle == "eta(1475)")
    {
        j = 328;
    }
    if (particle == "eta2(1645)")
    {
        j = 329;
    }
    if (particle == "eta2(1870)")
    {
        j = 330;
    }
    if (particle == "eta_prime")
    {
        j = 331;
    }
    if (particle == "etac")
    {
        j = 332;
    }
    if (particle == "f0(1370)")
    {
        j = 333;
    }
    if (particle == "f0(1500)")
    {
        j = 334;
    }
    if (particle == "f0(1710)")
    {
        j = 335;
    }
    if (particle == "f0(500)")
    {
        j = 336;
    }
    if (particle == "f0(980)")
    {
        j = 337;
    }
    if (particle == "f1(1285)")
    {
        j = 338;
    }
    if (particle == "f1(1420)")
    {
        j = 339;
    }
    if (particle == "f2(1270)")
    {
        j = 340;
    }
    if (particle == "f2(1810)")
    {
        j = 341;
    }
    if (particle == "f2(2010)")
    {
        j = 342;
    }
    if (particle == "f2_prime(1525)")
    {
        j = 343;
    }
    if (particle == "gamma")
    {
        j = 344;
    }
    if (particle == "geantino")
    {
        j = 345;
    }
    if (particle == "gluon")
    {
        j = 346;
    }
    if (particle == "h1(1170)")
    {
        j = 347;
    }
    if (particle == "h1(1380)")
    {
        j = 348;
    }
    if (particle == "hyperH4")
    {
        j = 349;
    }
    if (particle == "hyperHe5")
    {
        j = 350;
    }
    if (particle == "hyperalpha")
    {
        j = 351;
    }
    if (particle == "hypertriton")
    {
        j = 352;
    }
    if (particle == "k(1460)+")
    {
        j = 353;
    }
    if (particle == "k(1460)-")
    {
        j = 354;
    }
    if (particle == "k(1460)0")
    {
        j = 355;
    }
    if (particle == "k0_star(1430)+")
    {
        j = 356;
    }
    if (particle == "k0_star(1430)-")
    {
        j = 357;
    }
    if (particle == "k0_star(1430)0")
    {
        j = 358;
    }
    if (particle == "k1(1270)+")
    {
        j = 359;
    }
    if (particle == "k1(1270)-")
    {
        j = 360;
    }
    if (particle == "k1(1270)0")
    {
        j = 361;
    }
    if (particle == "k1(1400)+")
    {
        j = 362;
    }
    if (particle == "k1(1400)-")
    {
        j = 363;
    }
    if (particle == "k1(1400)0")
    {
        j = 364;
    }
    if (particle == "k2(1770)+")
    {
        j = 365;
    }
    if (particle == "k2(1770)-")
    {
        j = 366;
    }
    if (particle == "k2(1770)0")
    {
        j = 367;
    }
    if (particle == "k2_star(1430)+")
    {
        j = 368;
    }
    if (particle == "k2_star(1430)-")
    {
        j = 369;
    }
    if (particle == "k2_star(1430)0")
    {
        j = 370;
    }
    if (particle == "k2_star(1980)+")
    {
        j = 371;
    }
    if (particle == "k2_star(1980)-")
    {
        j = 372;
    }
    if (particle == "k2_star(1980)0")
    {
        j = 373;
    }
    if (particle == "k3_star(1780)+")
    {
        j = 374;
    }
    if (particle == "k3_star(1780)-")
    {
        j = 375;
    }
    if (particle == "k3_star(1780)0")
    {
        j = 376;
    }
    if (particle == "k_star(1410)+")
    {
        j = 377;
    }
    if (particle == "k_star(1410)-")
    {
        j = 378;
    }
    if (particle == "k_star(1410)0")
    {
        j = 379;
    }
    if (particle == "k_star(1680)+")
    {
        j = 380;
    }
    if (particle == "k_star(1680)-")
    {
        j = 381;
    }
    if (particle == "k_star(1680)0")
    {
        j = 382;
    }
    if (particle == "k_star+")
    {
        j = 383;
    }
    if (particle == "k_star-")
    {
        j = 384;
    }
    if (particle == "k_star0")
    {
        j = 385;
    }
    if (particle == "kaon+")
    {
        j = 386;
    }
    if (particle == "kaon-")
    {
        j = 387;
    }
    if (particle == "kaon0")
    {
        j = 388;
    }
    if (particle == "kaon0L")
    {
        j = 389;
    }
    if (particle == "kaon0S")
    {
        j = 390;
    }
    if (particle == "lambda")
    {
        j = 391;
    }
    if (particle == "lambda(1405)")
    {
        j = 392;
    }
    if (particle == "lambda(1520)")
    {
        j = 393;
    }
    if (particle == "lambda(1600)")
    {
        j = 394;
    }
    if (particle == "lambda(1670)")
    {
        j = 395;
    }
    if (particle == "lambda(1690)")
    {
        j = 396;
    }
    if (particle == "lambda(1800)")
    {
        j = 397;
    }
    if (particle == "lambda(1810)")
    {
        j = 398;
    }
    if (particle == "lambda(1820)")
    {
        j = 399;
    }
    if (particle == "lambda(1830)")
    {
        j = 400;
    }
    if (particle == "lambda(1890)")
    {
        j = 401;
    }
    if (particle == "lambda(2100)")
    {
        j = 402;
    }
    if (particle == "lambda(2110)")
    {
        j = 403;
    }
    if (particle == "lambda_b")
    {
        j = 404;
    }
    if (particle == "lambda_c+")
    {
        j = 405;
    }
    if (particle == "mu+")
    {
        j = 406;
    }
    if (particle == "mu-")
    {
        j = 407;
    }
    if (particle == "neutron")
    {
        j = 408;
    }
    if (particle == "nu_e")
    {
        j = 409;
    }
    if (particle == "nu_mu")
    {
        j = 410;
    }
    if (particle == "nu_tau")
    {
        j = 411;
    }
    if (particle == "omega")
    {
        j = 412;
    }
    if (particle == "omega(1420)")
    {
        j = 413;
    }
    if (particle == "omega(1650)")
    {
        j = 414;
    }
    if (particle == "omega-")
    {
        j = 415;
    }
    if (particle == "omega3(1670)")
    {
        j = 416;
    }
    if (particle == "omega_b-")
    {
        j = 417;
    }
    if (particle == "omega_c0")
    {
        j = 418;
    }
    if (particle == "opticalphoton")
    {
        j = 419;
    }
    if (particle == "phi")
    {
        j = 420;
    }
    if (particle == "phi(1680)")
    {
        j = 421;
    }
    if (particle == "phi3(1850)")
    {
        j = 422;
    }
    if (particle == "pi(1300)+")
    {
        j = 423;
    }
    if (particle == "pi(1300)-")
    {
        j = 424;
    }
    if (particle == "pi(1300)0")
    {
        j = 425;
    }
    if (particle == "pi+")
    {
        j = 426;
    }
    if (particle == "pi-")
    {
        j = 427;
    }
    if (particle == "pi0")
    {
        j = 428;
    }
    if (particle == "pi2(1670)+")
    {
        j = 429;
    }
    if (particle == "pi2(1670)-")
    {
        j = 430;
    }
    if (particle == "pi2(1670)0")
    {
        j = 431;
    }
    if (particle == "proton")
    {
        j = 432;
    }
    if (particle == "rho(1450)+")
    {
        j = 433;
    }
    if (particle == "rho(1450)-")
    {
        j = 434;
    }
    if (particle == "rho(1450)0")
    {
        j = 435;
    }
    if (particle == "rho(1700)+")
    {
        j = 436;
    }
    if (particle == "rho(1700)-")
    {
        j = 437;
    }
    if (particle == "rho(1700)0")
    {
        j = 438;
    }
    if (particle == "rho+")
    {
        j = 439;
    }
    if (particle == "rho-")
    {
        j = 440;
    }
    if (particle == "rho0")
    {
        j = 441;
    }
    if (particle == "rho3(1690)+")
    {
        j = 442;
    }
    if (particle == "rho3(1690)-")
    {
        j = 443;
    }
    if (particle == "rho3(1690)0")
    {
        j = 444;
    }
    if (particle == "s_quark")
    {
        j = 445;
    }
    if (particle == "sd0_diquark")
    {
        j = 446;
    }
    if (particle == "sd1_diquark")
    {
        j = 447;
    }
    if (particle == "sigma(1385)+")
    {
        j = 448;
    }
    if (particle == "sigma(1385)-")
    {
        j = 449;
    }
    if (particle == "sigma(1385)0")
    {
        j = 450;
    }
    if (particle == "sigma(1660)+")
    {
        j = 451;
    }
    if (particle == "sigma(1660)-")
    {
        j = 452;
    }
    if (particle == "sigma(1660)0")
    {
        j = 453;
    }
    if (particle == "sigma(1670)+")
    {
        j = 454;
    }
    if (particle == "sigma(1670)-")
    {
        j = 455;
    }
    if (particle == "sigma(1670)0")
    {
        j = 456;
    }
    if (particle == "sigma(1750)+")
    {
        j = 457;
    }
    if (particle == "sigma(1750)-")
    {
        j = 458;
    }
    if (particle == "sigma(1750)0")
    {
        j = 459;
    }
    if (particle == "sigma(1775)+")
    {
        j = 460;
    }
    if (particle == "sigma(1775)-")
    {
        j = 461;
    }
    if (particle == "sigma(1775)0")
    {
        j = 462;
    }
    if (particle == "sigma(1915)+")
    {
        j = 463;
    }
    if (particle == "sigma(1915)-")
    {
        j = 464;
    }
    if (particle == "sigma(1915)0")
    {
        j = 465;
    }
    if (particle == "sigma(1940)+")
    {
        j = 466;
    }
    if (particle == "sigma(1940)-")
    {
        j = 467;
    }
    if (particle == "sigma(1940)0")
    {
        j = 468;
    }
    if (particle == "sigma(2030)+")
    {
        j = 469;
    }
    if (particle == "sigma(2030)-")
    {
        j = 470;
    }
    if (particle == "sigma(2030)0")
    {
        j = 471;
    }
    if (particle == "sigma+")
    {
        j = 472;
    }
    if (particle == "sigma-")
    {
        j = 473;
    }
    if (particle == "sigma0")
    {
        j = 474;
    }
    if (particle == "sigma_b+")
    {
        j = 475;
    }
    if (particle == "sigma_b-")
    {
        j = 476;
    }
    if (particle == "sigma_b0")
    {
        j = 477;
    }
    if (particle == "sigma_c+")
    {
        j = 478;
    }
    if (particle == "sigma_c++")
    {
        j = 479;
    }
    if (particle == "sigma_c0")
    {
        j = 480;
    }
    if (particle == "ss1_diquark")
    {
        j = 481;
    }
    if (particle == "su0_diquark")
    {
        j = 482;
    }
    if (particle == "su1_diquark")
    {
        j = 483;
    }
    if (particle == "t_quark")
    {
        j = 484;
    }
    if (particle == "tau+")
    {
        j = 485;
    }
    if (particle == "tau-")
    {
        j = 486;
    }
    if (particle == "triton")
    {
        j = 487;
    }
    if (particle == "u_quark")
    {
        j = 488;
    }
    if (particle == "ud0_diquark")
    {
        j = 489;
    }
    if (particle == "ud1_diquark")
    {
        j = 490;
    }
    if (particle == "uu1_diquark")
    {
        j = 491;
    }
    if (particle == "xi(1530)-")
    {
        j = 492;
    }
    if (particle == "xi(1530)0")
    {
        j = 493;
    }
    if (particle == "xi(1690)-")
    {
        j = 494;
    }
    if (particle == "xi(1690)0")
    {
        j = 495;
    }
    if (particle == "xi(1820)-")
    {
        j = 496;
    }
    if (particle == "xi(1820)0")
    {
        j = 497;
    }
    if (particle == "xi(1950)-")
    {
        j = 498;
    }
    if (particle == "xi(1950)0")
    {
        j = 499;
    }
    if (particle == "xi(2030)-")
    {
        j = 500;
    }
    if (particle == "xi(2030)0")
    {
        j = 501;
    }
    if (particle == "xi-")
    {
        j = 502;
    }
    if (particle == "xi0")
    {
        j = 503;
    }
    if (particle == "xi_b-")
    {
        j = 504;
    }
    if (particle == "xi_b0")
    {
        j = 505;
    }
    if (particle == "xi_c+")
    {
        j = 506;
    }
    if (particle == "xi_c0")
    {
        j = 507;
    }
    output = net_number_detectors + 2 + j;
    /*
   net_number_detectors=net_number_detectors + 1;
   if(particle=="mu+"||particle=="mu-"){
       output = net_number_detectors + 1;

       }else if(particle=="e+"||particle=="e-"){

       output = net_number_detectors+2;

       }else{//others non targeted

       output = net_number_detectors+3;

       }
    */
}

void UDInCodeControllingSurface::Particles_no_finder(G4String particlet, int &no)
{
    if (particlet == "B+")
    {
        no = 0;
    }
    if (particlet == "B-")
    {
        no = 1;
    }
    if (particlet == "B0")
    {
        no = 2;
    }
    if (particlet == "Bc+")
    {
        no = 3;
    }
    if (particlet == "Bc-")
    {
        no = 4;
    }
    if (particlet == "Bs0")
    {
        no = 5;
    }
    if (particlet == "D+")
    {
        no = 6;
    }
    if (particlet == "D-")
    {
        no = 7;
    }
    if (particlet == "D0")
    {
        no = 8;
    }
    if (particlet == "Ds+")
    {
        no = 9;
    }
    if (particlet == "Ds-")
    {
        no = 10;
    }
    if (particlet == "GenericIon")
    {
        no = 11;
    }
    if (particlet == "He3")
    {
        no = 12;
    }
    if (particlet == "J/psi")
    {
        no = 13;
    }
    if (particlet == "N(1440)+")
    {
        no = 14;
    }
    if (particlet == "N(1440)0")
    {
        no = 15;
    }
    if (particlet == "N(1520)+")
    {
        no = 16;
    }
    if (particlet == "N(1520)0")
    {
        no = 17;
    }
    if (particlet == "N(1535)+")
    {
        no = 18;
    }
    if (particlet == "N(1535)0")
    {
        no = 19;
    }
    if (particlet == "N(1650)+")
    {
        no = 20;
    }
    if (particlet == "N(1650)0")
    {
        no = 21;
    }
    if (particlet == "N(1675)+")
    {
        no = 22;
    }
    if (particlet == "N(1675)0")
    {
        no = 23;
    }
    if (particlet == "N(1680)+")
    {
        no = 24;
    }
    if (particlet == "N(1680)0")
    {
        no = 25;
    }
    if (particlet == "N(1700)+")
    {
        no = 26;
    }
    if (particlet == "N(1700)0")
    {
        no = 27;
    }
    if (particlet == "N(1710)+")
    {
        no = 28;
    }
    if (particlet == "N(1710)0")
    {
        no = 29;
    }
    if (particlet == "N(1720)+")
    {
        no = 30;
    }
    if (particlet == "N(1720)0")
    {
        no = 31;
    }
    if (particlet == "N(1900)+")
    {
        no = 32;
    }
    if (particlet == "N(1900)0")
    {
        no = 33;
    }
    if (particlet == "N(1990)+")
    {
        no = 34;
    }
    if (particlet == "N(1990)0")
    {
        no = 35;
    }
    if (particlet == "N(2090)+")
    {
        no = 36;
    }
    if (particlet == "N(2090)0")
    {
        no = 37;
    }
    if (particlet == "N(2190)+")
    {
        no = 38;
    }
    if (particlet == "N(2190)0")
    {
        no = 39;
    }
    if (particlet == "N(2220)+")
    {
        no = 40;
    }
    if (particlet == "N(2220)0")
    {
        no = 41;
    }
    if (particlet == "N(2250)+")
    {
        no = 42;
    }
    if (particlet == "N(2250)0")
    {
        no = 43;
    }
    if (particlet == "Upsilon")
    {
        no = 44;
    }
    if (particlet == "a0(1450)+")
    {
        no = 45;
    }
    if (particlet == "a0(1450)-")
    {
        no = 46;
    }
    if (particlet == "a0(1450)0")
    {
        no = 47;
    }
    if (particlet == "a0(980)+")
    {
        no = 48;
    }
    if (particlet == "a0(980)-")
    {
        no = 49;
    }
    if (particlet == "a0(980)0")
    {
        no = 50;
    }
    if (particlet == "a1(1260)+")
    {
        no = 51;
    }
    if (particlet == "a1(1260)-")
    {
        no = 52;
    }
    if (particlet == "a1(1260)0")
    {
        no = 53;
    }
    if (particlet == "a2(1320)+")
    {
        no = 54;
    }
    if (particlet == "a2(1320)-")
    {
        no = 55;
    }
    if (particlet == "a2(1320)0")
    {
        no = 56;
    }
    if (particlet == "alpha")
    {
        no = 57;
    }
    if (particlet == "anti_B0")
    {
        no = 58;
    }
    if (particlet == "anti_Bs0")
    {
        no = 59;
    }
    if (particlet == "anti_D0")
    {
        no = 60;
    }
    if (particlet == "anti_He3")
    {
        no = 61;
    }
    if (particlet == "anti_N(1440)+")
    {
        no = 62;
    }
    if (particlet == "anti_N(1440)0")
    {
        no = 63;
    }
    if (particlet == "anti_N(1520)+")
    {
        no = 64;
    }
    if (particlet == "anti_N(1520)0")
    {
        no = 65;
    }
    if (particlet == "anti_N(1535)+")
    {
        no = 66;
    }
    if (particlet == "anti_N(1535)0")
    {
        no = 67;
    }
    if (particlet == "anti_N(1650)+")
    {
        no = 68;
    }
    if (particlet == "anti_N(1650)0")
    {
        no = 69;
    }
    if (particlet == "anti_N(1675)+")
    {
        no = 70;
    }
    if (particlet == "anti_N(1675)0")
    {
        no = 71;
    }
    if (particlet == "anti_N(1680)+")
    {
        no = 72;
    }
    if (particlet == "anti_N(1680)0")
    {
        no = 73;
    }
    if (particlet == "anti_N(1700)+")
    {
        no = 74;
    }
    if (particlet == "anti_N(1700)0")
    {
        no = 75;
    }
    if (particlet == "anti_N(1710)+")
    {
        no = 76;
    }
    if (particlet == "anti_N(1710)0")
    {
        no = 77;
    }
    if (particlet == "anti_N(1720)+")
    {
        no = 78;
    }
    if (particlet == "anti_N(1720)0")
    {
        no = 79;
    }
    if (particlet == "anti_N(1900)+")
    {
        no = 80;
    }
    if (particlet == "anti_N(1900)0")
    {
        no = 81;
    }
    if (particlet == "anti_N(1990)+")
    {
        no = 82;
    }
    if (particlet == "anti_N(1990)0")
    {
        no = 83;
    }
    if (particlet == "anti_N(2090)+")
    {
        no = 84;
    }
    if (particlet == "anti_N(2090)0")
    {
        no = 85;
    }
    if (particlet == "anti_N(2190)+")
    {
        no = 86;
    }
    if (particlet == "anti_N(2190)0")
    {
        no = 87;
    }
    if (particlet == "anti_N(2220)+")
    {
        no = 88;
    }
    if (particlet == "anti_N(2220)0")
    {
        no = 89;
    }
    if (particlet == "anti_N(2250)+")
    {
        no = 90;
    }
    if (particlet == "anti_N(2250)0")
    {
        no = 91;
    }
    if (particlet == "anti_alpha")
    {
        no = 92;
    }
    if (particlet == "anti_b_quark")
    {
        no = 93;
    }
    if (particlet == "anti_bb1_diquark")
    {
        no = 94;
    }
    if (particlet == "anti_bc0_diquark")
    {
        no = 95;
    }
    if (particlet == "anti_bc1_diquark")
    {
        no = 96;
    }
    if (particlet == "anti_bd0_diquark")
    {
        no = 97;
    }
    if (particlet == "anti_bd1_diquark")
    {
        no = 98;
    }
    if (particlet == "anti_bs0_diquark")
    {
        no = 99;
    }
    if (particlet == "anti_bs1_diquark")
    {
        no = 100;
    }
    if (particlet == "anti_bu0_diquark")
    {
        no = 101;
    }
    if (particlet == "anti_bu1_diquark")
    {
        no = 102;
    }
    if (particlet == "anti_c_quark")
    {
        no = 103;
    }
    if (particlet == "anti_cc1_diquark")
    {
        no = 104;
    }
    if (particlet == "anti_cd0_diquark")
    {
        no = 105;
    }
    if (particlet == "anti_cd1_diquark")
    {
        no = 106;
    }
    if (particlet == "anti_cs0_diquark")
    {
        no = 107;
    }
    if (particlet == "anti_cs1_diquark")
    {
        no = 108;
    }
    if (particlet == "anti_cu0_diquark")
    {
        no = 109;
    }
    if (particlet == "anti_cu1_diquark")
    {
        no = 110;
    }
    if (particlet == "anti_d_quark")
    {
        no = 111;
    }
    if (particlet == "anti_dd1_diquark")
    {
        no = 112;
    }
    if (particlet == "anti_delta(1600)+")
    {
        no = 113;
    }
    if (particlet == "anti_delta(1600)++")
    {
        no = 114;
    }
    if (particlet == "anti_delta(1600)-")
    {
        no = 115;
    }
    if (particlet == "anti_delta(1600)0")
    {
        no = 116;
    }
    if (particlet == "anti_delta(1620)+")
    {
        no = 117;
    }
    if (particlet == "anti_delta(1620)++")
    {
        no = 118;
    }
    if (particlet == "anti_delta(1620)-")
    {
        no = 119;
    }
    if (particlet == "anti_delta(1620)0")
    {
        no = 120;
    }
    if (particlet == "anti_delta(1700)+")
    {
        no = 121;
    }
    if (particlet == "anti_delta(1700)++")
    {
        no = 122;
    }
    if (particlet == "anti_delta(1700)-")
    {
        no = 123;
    }
    if (particlet == "anti_delta(1700)0")
    {
        no = 124;
    }
    if (particlet == "anti_delta(1900)+")
    {
        no = 125;
    }
    if (particlet == "anti_delta(1900)++")
    {
        no = 126;
    }
    if (particlet == "anti_delta(1900)-")
    {
        no = 127;
    }
    if (particlet == "anti_delta(1900)0")
    {
        no = 128;
    }
    if (particlet == "anti_delta(1905)+")
    {
        no = 129;
    }
    if (particlet == "anti_delta(1905)++")
    {
        no = 130;
    }
    if (particlet == "anti_delta(1905)-")
    {
        no = 131;
    }
    if (particlet == "anti_delta(1905)0")
    {
        no = 132;
    }
    if (particlet == "anti_delta(1910)+")
    {
        no = 133;
    }
    if (particlet == "anti_delta(1910)++")
    {
        no = 134;
    }
    if (particlet == "anti_delta(1910)-")
    {
        no = 135;
    }
    if (particlet == "anti_delta(1910)0")
    {
        no = 136;
    }
    if (particlet == "anti_delta(1920)+")
    {
        no = 137;
    }
    if (particlet == "anti_delta(1920)++")
    {
        no = 138;
    }
    if (particlet == "anti_delta(1920)-")
    {
        no = 139;
    }
    if (particlet == "anti_delta(1920)0")
    {
        no = 140;
    }
    if (particlet == "anti_delta(1930)+")
    {
        no = 141;
    }
    if (particlet == "anti_delta(1930)++")
    {
        no = 142;
    }
    if (particlet == "anti_delta(1930)-")
    {
        no = 143;
    }
    if (particlet == "anti_delta(1930)0")
    {
        no = 144;
    }
    if (particlet == "anti_delta(1950)+")
    {
        no = 145;
    }
    if (particlet == "anti_delta(1950)++")
    {
        no = 146;
    }
    if (particlet == "anti_delta(1950)-")
    {
        no = 147;
    }
    if (particlet == "anti_delta(1950)0")
    {
        no = 148;
    }
    if (particlet == "anti_delta+")
    {
        no = 149;
    }
    if (particlet == "anti_delta++")
    {
        no = 150;
    }
    if (particlet == "anti_delta-")
    {
        no = 151;
    }
    if (particlet == "anti_delta0")
    {
        no = 152;
    }
    if (particlet == "anti_deuteron")
    {
        no = 153;
    }
    if (particlet == "anti_doublehyperH4")
    {
        no = 154;
    }
    if (particlet == "anti_doublehyperdoubleneutron")
    {
        no = 155;
    }
    if (particlet == "anti_hyperH4")
    {
        no = 156;
    }
    if (particlet == "anti_hyperHe5")
    {
        no = 157;
    }
    if (particlet == "anti_hyperalpha")
    {
        no = 158;
    }
    if (particlet == "anti_hypertriton")
    {
        no = 159;
    }
    if (particlet == "anti_k(1460)0")
    {
        no = 160;
    }
    if (particlet == "anti_k0_star(1430)0")
    {
        no = 161;
    }
    if (particlet == "anti_k1(1270)0")
    {
        no = 162;
    }
    if (particlet == "anti_k1(1400)0")
    {
        no = 163;
    }
    if (particlet == "anti_k2(1770)0")
    {
        no = 164;
    }
    if (particlet == "anti_k2_star(1430)0")
    {
        no = 165;
    }
    if (particlet == "anti_k2_star(1980)0")
    {
        no = 166;
    }
    if (particlet == "anti_k3_star(1780)0")
    {
        no = 167;
    }
    if (particlet == "anti_k_star(1410)0")
    {
        no = 168;
    }
    if (particlet == "anti_k_star(1680)0")
    {
        no = 169;
    }
    if (particlet == "anti_k_star0")
    {
        no = 170;
    }
    if (particlet == "anti_kaon0")
    {
        no = 171;
    }
    if (particlet == "anti_lambda")
    {
        no = 172;
    }
    if (particlet == "anti_lambda(1405)")
    {
        no = 173;
    }
    if (particlet == "anti_lambda(1520)")
    {
        no = 174;
    }
    if (particlet == "anti_lambda(1600)")
    {
        no = 175;
    }
    if (particlet == "anti_lambda(1670)")
    {
        no = 176;
    }
    if (particlet == "anti_lambda(1690)")
    {
        no = 177;
    }
    if (particlet == "anti_lambda(1800)")
    {
        no = 178;
    }
    if (particlet == "anti_lambda(1810)")
    {
        no = 179;
    }
    if (particlet == "anti_lambda(1820)")
    {
        no = 180;
    }
    if (particlet == "anti_lambda(1830)")
    {
        no = 181;
    }
    if (particlet == "anti_lambda(1890)")
    {
        no = 182;
    }
    if (particlet == "anti_lambda(2100)")
    {
        no = 183;
    }
    if (particlet == "anti_lambda(2110)")
    {
        no = 184;
    }
    if (particlet == "anti_lambda_b")
    {
        no = 185;
    }
    if (particlet == "anti_lambda_c+")
    {
        no = 186;
    }
    if (particlet == "anti_neutron")
    {
        no = 187;
    }
    if (particlet == "anti_nu_e")
    {
        no = 188;
    }
    if (particlet == "anti_nu_mu")
    {
        no = 189;
    }
    if (particlet == "anti_nu_tau")
    {
        no = 190;
    }
    if (particlet == "anti_omega-")
    {
        no = 191;
    }
    if (particlet == "anti_omega_b-")
    {
        no = 192;
    }
    if (particlet == "anti_omega_c0")
    {
        no = 193;
    }
    if (particlet == "anti_proton")
    {
        no = 194;
    }
    if (particlet == "anti_s_quark")
    {
        no = 195;
    }
    if (particlet == "anti_sd0_diquark")
    {
        no = 196;
    }
    if (particlet == "anti_sd1_diquark")
    {
        no = 197;
    }
    if (particlet == "anti_sigma(1385)+")
    {
        no = 198;
    }
    if (particlet == "anti_sigma(1385)-")
    {
        no = 199;
    }
    if (particlet == "anti_sigma(1385)0")
    {
        no = 200;
    }
    if (particlet == "anti_sigma(1660)+")
    {
        no = 201;
    }
    if (particlet == "anti_sigma(1660)-")
    {
        no = 202;
    }
    if (particlet == "anti_sigma(1660)0")
    {
        no = 203;
    }
    if (particlet == "anti_sigma(1670)+")
    {
        no = 204;
    }
    if (particlet == "anti_sigma(1670)-")
    {
        no = 205;
    }
    if (particlet == "anti_sigma(1670)0")
    {
        no = 206;
    }
    if (particlet == "anti_sigma(1750)+")
    {
        no = 207;
    }
    if (particlet == "anti_sigma(1750)-")
    {
        no = 208;
    }
    if (particlet == "anti_sigma(1750)0")
    {
        no = 209;
    }
    if (particlet == "anti_sigma(1775)+")
    {
        no = 210;
    }
    if (particlet == "anti_sigma(1775)-")
    {
        no = 211;
    }
    if (particlet == "anti_sigma(1775)0")
    {
        no = 212;
    }
    if (particlet == "anti_sigma(1915)+")
    {
        no = 213;
    }
    if (particlet == "anti_sigma(1915)-")
    {
        no = 214;
    }
    if (particlet == "anti_sigma(1915)0")
    {
        no = 215;
    }
    if (particlet == "anti_sigma(1940)+")
    {
        no = 216;
    }
    if (particlet == "anti_sigma(1940)-")
    {
        no = 217;
    }
    if (particlet == "anti_sigma(1940)0")
    {
        no = 218;
    }
    if (particlet == "anti_sigma(2030)+")
    {
        no = 219;
    }
    if (particlet == "anti_sigma(2030)-")
    {
        no = 220;
    }
    if (particlet == "anti_sigma(2030)0")
    {
        no = 221;
    }
    if (particlet == "anti_sigma+")
    {
        no = 222;
    }
    if (particlet == "anti_sigma-")
    {
        no = 223;
    }
    if (particlet == "anti_sigma0")
    {
        no = 224;
    }
    if (particlet == "anti_sigma_b+")
    {
        no = 225;
    }
    if (particlet == "anti_sigma_b-")
    {
        no = 226;
    }
    if (particlet == "anti_sigma_b0")
    {
        no = 227;
    }
    if (particlet == "anti_sigma_c+")
    {
        no = 228;
    }
    if (particlet == "anti_sigma_c++")
    {
        no = 229;
    }
    if (particlet == "anti_sigma_c0")
    {
        no = 230;
    }
    if (particlet == "anti_ss1_diquark")
    {
        no = 231;
    }
    if (particlet == "anti_su0_diquark")
    {
        no = 232;
    }
    if (particlet == "anti_su1_diquark")
    {
        no = 233;
    }
    if (particlet == "anti_t_quark")
    {
        no = 234;
    }
    if (particlet == "anti_triton")
    {
        no = 235;
    }
    if (particlet == "anti_u_quark")
    {
        no = 236;
    }
    if (particlet == "anti_ud0_diquark")
    {
        no = 237;
    }
    if (particlet == "anti_ud1_diquark")
    {
        no = 238;
    }
    if (particlet == "anti_uu1_diquark")
    {
        no = 239;
    }
    if (particlet == "anti_xi(1530)-")
    {
        no = 240;
    }
    if (particlet == "anti_xi(1530)0")
    {
        no = 241;
    }
    if (particlet == "anti_xi(1690)-")
    {
        no = 242;
    }
    if (particlet == "anti_xi(1690)0")
    {
        no = 243;
    }
    if (particlet == "anti_xi(1820)-")
    {
        no = 244;
    }
    if (particlet == "anti_xi(1820)0")
    {
        no = 245;
    }
    if (particlet == "anti_xi(1950)-")
    {
        no = 246;
    }
    if (particlet == "anti_xi(1950)0")
    {
        no = 247;
    }
    if (particlet == "anti_xi(2030)-")
    {
        no = 248;
    }
    if (particlet == "anti_xi(2030)0")
    {
        no = 249;
    }
    if (particlet == "anti_xi-")
    {
        no = 250;
    }
    if (particlet == "anti_xi0")
    {
        no = 251;
    }
    if (particlet == "anti_xi_b-")
    {
        no = 252;
    }
    if (particlet == "anti_xi_b0")
    {
        no = 253;
    }
    if (particlet == "anti_xi_c+")
    {
        no = 254;
    }
    if (particlet == "anti_xi_c0")
    {
        no = 255;
    }
    if (particlet == "b1(1235)+")
    {
        no = 256;
    }
    if (particlet == "b1(1235)-")
    {
        no = 257;
    }
    if (particlet == "b1(1235)0")
    {
        no = 258;
    }
    if (particlet == "b_quark")
    {
        no = 259;
    }
    if (particlet == "bb1_diquark")
    {
        no = 260;
    }
    if (particlet == "bc0_diquark")
    {
        no = 261;
    }
    if (particlet == "bc1_diquark")
    {
        no = 262;
    }
    if (particlet == "bd0_diquark")
    {
        no = 263;
    }
    if (particlet == "bd1_diquark")
    {
        no = 264;
    }
    if (particlet == "bs0_diquark")
    {
        no = 265;
    }
    if (particlet == "bs1_diquark")
    {
        no = 266;
    }
    if (particlet == "bu0_diquark")
    {
        no = 267;
    }
    if (particlet == "bu1_diquark")
    {
        no = 268;
    }
    if (particlet == "c_quark")
    {
        no = 269;
    }
    if (particlet == "cc1_diquark")
    {
        no = 270;
    }
    if (particlet == "cd0_diquark")
    {
        no = 271;
    }
    if (particlet == "cd1_diquark")
    {
        no = 272;
    }
    if (particlet == "chargedgeantino")
    {
        no = 273;
    }
    if (particlet == "cs0_diquark")
    {
        no = 274;
    }
    if (particlet == "cs1_diquark")
    {
        no = 275;
    }
    if (particlet == "cu0_diquark")
    {
        no = 276;
    }
    if (particlet == "cu1_diquark")
    {
        no = 277;
    }
    if (particlet == "d_quark")
    {
        no = 278;
    }
    if (particlet == "dd1_diquark")
    {
        no = 279;
    }
    if (particlet == "delta(1600)+")
    {
        no = 280;
    }
    if (particlet == "delta(1600)++")
    {
        no = 281;
    }
    if (particlet == "delta(1600)-")
    {
        no = 282;
    }
    if (particlet == "delta(1600)0")
    {
        no = 283;
    }
    if (particlet == "delta(1620)+")
    {
        no = 284;
    }
    if (particlet == "delta(1620)++")
    {
        no = 285;
    }
    if (particlet == "delta(1620)-")
    {
        no = 286;
    }
    if (particlet == "delta(1620)0")
    {
        no = 287;
    }
    if (particlet == "delta(1700)+")
    {
        no = 288;
    }
    if (particlet == "delta(1700)++")
    {
        no = 289;
    }
    if (particlet == "delta(1700)-")
    {
        no = 290;
    }
    if (particlet == "delta(1700)0")
    {
        no = 291;
    }
    if (particlet == "delta(1900)+")
    {
        no = 292;
    }
    if (particlet == "delta(1900)++")
    {
        no = 293;
    }
    if (particlet == "delta(1900)-")
    {
        no = 294;
    }
    if (particlet == "delta(1900)0")
    {
        no = 295;
    }
    if (particlet == "delta(1905)+")
    {
        no = 296;
    }
    if (particlet == "delta(1905)++")
    {
        no = 297;
    }
    if (particlet == "delta(1905)-")
    {
        no = 298;
    }
    if (particlet == "delta(1905)0")
    {
        no = 299;
    }
    if (particlet == "delta(1910)+")
    {
        no = 300;
    }
    if (particlet == "delta(1910)++")
    {
        no = 301;
    }
    if (particlet == "delta(1910)-")
    {
        no = 302;
    }
    if (particlet == "delta(1910)0")
    {
        no = 303;
    }
    if (particlet == "delta(1920)+")
    {
        no = 304;
    }
    if (particlet == "delta(1920)++")
    {
        no = 305;
    }
    if (particlet == "delta(1920)-")
    {
        no = 306;
    }
    if (particlet == "delta(1920)0")
    {
        no = 307;
    }
    if (particlet == "delta(1930)+")
    {
        no = 308;
    }
    if (particlet == "delta(1930)++")
    {
        no = 309;
    }
    if (particlet == "delta(1930)-")
    {
        no = 310;
    }
    if (particlet == "delta(1930)0")
    {
        no = 311;
    }
    if (particlet == "delta(1950)+")
    {
        no = 312;
    }
    if (particlet == "delta(1950)++")
    {
        no = 313;
    }
    if (particlet == "delta(1950)-")
    {
        no = 314;
    }
    if (particlet == "delta(1950)0")
    {
        no = 315;
    }
    if (particlet == "delta+")
    {
        no = 316;
    }
    if (particlet == "delta++")
    {
        no = 317;
    }
    if (particlet == "delta-")
    {
        no = 318;
    }
    if (particlet == "delta0")
    {
        no = 319;
    }
    if (particlet == "deuteron")
    {
        no = 320;
    }
    if (particlet == "doublehyperH4")
    {
        no = 321;
    }
    if (particlet == "doublehyperdoubleneutron")
    {
        no = 322;
    }
    if (particlet == "e+")
    {
        no = 323;
    }
    if (particlet == "e-")
    {
        no = 324;
    }
    if (particlet == "eta")
    {
        no = 325;
    }
    if (particlet == "eta(1295)")
    {
        no = 326;
    }
    if (particlet == "eta(1405)")
    {
        no = 327;
    }
    if (particlet == "eta(1475)")
    {
        no = 328;
    }
    if (particlet == "eta2(1645)")
    {
        no = 329;
    }
    if (particlet == "eta2(1870)")
    {
        no = 330;
    }
    if (particlet == "eta_prime")
    {
        no = 331;
    }
    if (particlet == "etac")
    {
        no = 332;
    }
    if (particlet == "f0(1370)")
    {
        no = 333;
    }
    if (particlet == "f0(1500)")
    {
        no = 334;
    }
    if (particlet == "f0(1710)")
    {
        no = 335;
    }
    if (particlet == "f0(500)")
    {
        no = 336;
    }
    if (particlet == "f0(980)")
    {
        no = 337;
    }
    if (particlet == "f1(1285)")
    {
        no = 338;
    }
    if (particlet == "f1(1420)")
    {
        no = 339;
    }
    if (particlet == "f2(1270)")
    {
        no = 340;
    }
    if (particlet == "f2(1810)")
    {
        no = 341;
    }
    if (particlet == "f2(2010)")
    {
        no = 342;
    }
    if (particlet == "f2_prime(1525)")
    {
        no = 343;
    }
    if (particlet == "gamma")
    {
        no = 344;
    }
    if (particlet == "geantino")
    {
        no = 345;
    }
    if (particlet == "gluon")
    {
        no = 346;
    }
    if (particlet == "h1(1170)")
    {
        no = 347;
    }
    if (particlet == "h1(1380)")
    {
        no = 348;
    }
    if (particlet == "hyperH4")
    {
        no = 349;
    }
    if (particlet == "hyperHe5")
    {
        no = 350;
    }
    if (particlet == "hyperalpha")
    {
        no = 351;
    }
    if (particlet == "hypertriton")
    {
        no = 352;
    }
    if (particlet == "k(1460)+")
    {
        no = 353;
    }
    if (particlet == "k(1460)-")
    {
        no = 354;
    }
    if (particlet == "k(1460)0")
    {
        no = 355;
    }
    if (particlet == "k0_star(1430)+")
    {
        no = 356;
    }
    if (particlet == "k0_star(1430)-")
    {
        no = 357;
    }
    if (particlet == "k0_star(1430)0")
    {
        no = 358;
    }
    if (particlet == "k1(1270)+")
    {
        no = 359;
    }
    if (particlet == "k1(1270)-")
    {
        no = 360;
    }
    if (particlet == "k1(1270)0")
    {
        no = 361;
    }
    if (particlet == "k1(1400)+")
    {
        no = 362;
    }
    if (particlet == "k1(1400)-")
    {
        no = 363;
    }
    if (particlet == "k1(1400)0")
    {
        no = 364;
    }
    if (particlet == "k2(1770)+")
    {
        no = 365;
    }
    if (particlet == "k2(1770)-")
    {
        no = 366;
    }
    if (particlet == "k2(1770)0")
    {
        no = 367;
    }
    if (particlet == "k2_star(1430)+")
    {
        no = 368;
    }
    if (particlet == "k2_star(1430)-")
    {
        no = 369;
    }
    if (particlet == "k2_star(1430)0")
    {
        no = 370;
    }
    if (particlet == "k2_star(1980)+")
    {
        no = 371;
    }
    if (particlet == "k2_star(1980)-")
    {
        no = 372;
    }
    if (particlet == "k2_star(1980)0")
    {
        no = 373;
    }
    if (particlet == "k3_star(1780)+")
    {
        no = 374;
    }
    if (particlet == "k3_star(1780)-")
    {
        no = 375;
    }
    if (particlet == "k3_star(1780)0")
    {
        no = 376;
    }
    if (particlet == "k_star(1410)+")
    {
        no = 377;
    }
    if (particlet == "k_star(1410)-")
    {
        no = 378;
    }
    if (particlet == "k_star(1410)0")
    {
        no = 379;
    }
    if (particlet == "k_star(1680)+")
    {
        no = 380;
    }
    if (particlet == "k_star(1680)-")
    {
        no = 381;
    }
    if (particlet == "k_star(1680)0")
    {
        no = 382;
    }
    if (particlet == "k_star+")
    {
        no = 383;
    }
    if (particlet == "k_star-")
    {
        no = 384;
    }
    if (particlet == "k_star0")
    {
        no = 385;
    }
    if (particlet == "kaon+")
    {
        no = 386;
    }
    if (particlet == "kaon-")
    {
        no = 387;
    }
    if (particlet == "kaon0")
    {
        no = 388;
    }
    if (particlet == "kaon0L")
    {
        no = 389;
    }
    if (particlet == "kaon0S")
    {
        no = 390;
    }
    if (particlet == "lambda")
    {
        no = 391;
    }
    if (particlet == "lambda(1405)")
    {
        no = 392;
    }
    if (particlet == "lambda(1520)")
    {
        no = 393;
    }
    if (particlet == "lambda(1600)")
    {
        no = 394;
    }
    if (particlet == "lambda(1670)")
    {
        no = 395;
    }
    if (particlet == "lambda(1690)")
    {
        no = 396;
    }
    if (particlet == "lambda(1800)")
    {
        no = 397;
    }
    if (particlet == "lambda(1810)")
    {
        no = 398;
    }
    if (particlet == "lambda(1820)")
    {
        no = 399;
    }
    if (particlet == "lambda(1830)")
    {
        no = 400;
    }
    if (particlet == "lambda(1890)")
    {
        no = 401;
    }
    if (particlet == "lambda(2100)")
    {
        no = 402;
    }
    if (particlet == "lambda(2110)")
    {
        no = 403;
    }
    if (particlet == "lambda_b")
    {
        no = 404;
    }
    if (particlet == "lambda_c+")
    {
        no = 405;
    }
    if (particlet == "mu+")
    {
        no = 406;
    }
    if (particlet == "mu-")
    {
        no = 407;
    }
    if (particlet == "neutron")
    {
        no = 408;
    }
    if (particlet == "nu_e")
    {
        no = 409;
    }
    if (particlet == "nu_mu")
    {
        no = 410;
    }
    if (particlet == "nu_tau")
    {
        no = 411;
    }
    if (particlet == "omega")
    {
        no = 412;
    }
    if (particlet == "omega(1420)")
    {
        no = 413;
    }
    if (particlet == "omega(1650)")
    {
        no = 414;
    }
    if (particlet == "omega-")
    {
        no = 415;
    }
    if (particlet == "omega3(1670)")
    {
        no = 416;
    }
    if (particlet == "omega_b-")
    {
        no = 417;
    }
    if (particlet == "omega_c0")
    {
        no = 418;
    }
    if (particlet == "opticalphoton")
    {
        no = 419;
    }
    if (particlet == "phi")
    {
        no = 420;
    }
    if (particlet == "phi(1680)")
    {
        no = 421;
    }
    if (particlet == "phi3(1850)")
    {
        no = 422;
    }
    if (particlet == "pi(1300)+")
    {
        no = 423;
    }
    if (particlet == "pi(1300)-")
    {
        no = 424;
    }
    if (particlet == "pi(1300)0")
    {
        no = 425;
    }
    if (particlet == "pi+")
    {
        no = 426;
    }
    if (particlet == "pi-")
    {
        no = 427;
    }
    if (particlet == "pi0")
    {
        no = 428;
    }
    if (particlet == "pi2(1670)+")
    {
        no = 429;
    }
    if (particlet == "pi2(1670)-")
    {
        no = 430;
    }
    if (particlet == "pi2(1670)0")
    {
        no = 431;
    }
    if (particlet == "proton")
    {
        no = 432;
    }
    if (particlet == "rho(1450)+")
    {
        no = 433;
    }
    if (particlet == "rho(1450)-")
    {
        no = 434;
    }
    if (particlet == "rho(1450)0")
    {
        no = 435;
    }
    if (particlet == "rho(1700)+")
    {
        no = 436;
    }
    if (particlet == "rho(1700)-")
    {
        no = 437;
    }
    if (particlet == "rho(1700)0")
    {
        no = 438;
    }
    if (particlet == "rho+")
    {
        no = 439;
    }
    if (particlet == "rho-")
    {
        no = 440;
    }
    if (particlet == "rho0")
    {
        no = 441;
    }
    if (particlet == "rho3(1690)+")
    {
        no = 442;
    }
    if (particlet == "rho3(1690)-")
    {
        no = 443;
    }
    if (particlet == "rho3(1690)0")
    {
        no = 444;
    }
    if (particlet == "s_quark")
    {
        no = 445;
    }
    if (particlet == "sd0_diquark")
    {
        no = 446;
    }
    if (particlet == "sd1_diquark")
    {
        no = 447;
    }
    if (particlet == "sigma(1385)+")
    {
        no = 448;
    }
    if (particlet == "sigma(1385)-")
    {
        no = 449;
    }
    if (particlet == "sigma(1385)0")
    {
        no = 450;
    }
    if (particlet == "sigma(1660)+")
    {
        no = 451;
    }
    if (particlet == "sigma(1660)-")
    {
        no = 452;
    }
    if (particlet == "sigma(1660)0")
    {
        no = 453;
    }
    if (particlet == "sigma(1670)+")
    {
        no = 454;
    }
    if (particlet == "sigma(1670)-")
    {
        no = 455;
    }
    if (particlet == "sigma(1670)0")
    {
        no = 456;
    }
    if (particlet == "sigma(1750)+")
    {
        no = 457;
    }
    if (particlet == "sigma(1750)-")
    {
        no = 458;
    }
    if (particlet == "sigma(1750)0")
    {
        no = 459;
    }
    if (particlet == "sigma(1775)+")
    {
        no = 460;
    }
    if (particlet == "sigma(1775)-")
    {
        no = 461;
    }
    if (particlet == "sigma(1775)0")
    {
        no = 462;
    }
    if (particlet == "sigma(1915)+")
    {
        no = 463;
    }
    if (particlet == "sigma(1915)-")
    {
        no = 464;
    }
    if (particlet == "sigma(1915)0")
    {
        no = 465;
    }
    if (particlet == "sigma(1940)+")
    {
        no = 466;
    }
    if (particlet == "sigma(1940)-")
    {
        no = 467;
    }
    if (particlet == "sigma(1940)0")
    {
        no = 468;
    }
    if (particlet == "sigma(2030)+")
    {
        no = 469;
    }
    if (particlet == "sigma(2030)-")
    {
        no = 470;
    }
    if (particlet == "sigma(2030)0")
    {
        no = 471;
    }
    if (particlet == "sigma+")
    {
        no = 472;
    }
    if (particlet == "sigma-")
    {
        no = 473;
    }
    if (particlet == "sigma0")
    {
        no = 474;
    }
    if (particlet == "sigma_b+")
    {
        no = 475;
    }
    if (particlet == "sigma_b-")
    {
        no = 476;
    }
    if (particlet == "sigma_b0")
    {
        no = 477;
    }
    if (particlet == "sigma_c+")
    {
        no = 478;
    }
    if (particlet == "sigma_c++")
    {
        no = 479;
    }
    if (particlet == "sigma_c0")
    {
        no = 480;
    }
    if (particlet == "ss1_diquark")
    {
        no = 481;
    }
    if (particlet == "su0_diquark")
    {
        no = 482;
    }
    if (particlet == "su1_diquark")
    {
        no = 483;
    }
    if (particlet == "t_quark")
    {
        no = 484;
    }
    if (particlet == "tau+")
    {
        no = 485;
    }
    if (particlet == "tau-")
    {
        no = 486;
    }
    if (particlet == "triton")
    {
        no = 487;
    }
    if (particlet == "u_quark")
    {
        no = 488;
    }
    if (particlet == "ud0_diquark")
    {
        no = 489;
    }
    if (particlet == "ud1_diquark")
    {
        no = 490;
    }
    if (particlet == "uu1_diquark")
    {
        no = 491;
    }
    if (particlet == "xi(1530)-")
    {
        no = 492;
    }
    if (particlet == "xi(1530)0")
    {
        no = 493;
    }
    if (particlet == "xi(1690)-")
    {
        no = 494;
    }
    if (particlet == "xi(1690)0")
    {
        no = 495;
    }
    if (particlet == "xi(1820)-")
    {
        no = 496;
    }
    if (particlet == "xi(1820)0")
    {
        no = 497;
    }
    if (particlet == "xi(1950)-")
    {
        no = 498;
    }
    if (particlet == "xi(1950)0")
    {
        no = 499;
    }
    if (particlet == "xi(2030)-")
    {
        no = 500;
    }
    if (particlet == "xi(2030)0")
    {
        no = 501;
    }
    if (particlet == "xi-")
    {
        no = 502;
    }
    if (particlet == "xi0")
    {
        no = 503;
    }
    if (particlet == "xi_b-")
    {
        no = 504;
    }
    if (particlet == "xi_b0")
    {
        no = 505;
    }
    if (particlet == "xi_c+")
    {
        no = 506;
    }
    if (particlet == "xi_c0")
    {
        no = 507;
    }
}

void UDInCodeControllingSurface::Particles_name_finder(G4String &Particle, int no)
{

    if (no == 0)
    {
        Particle = "B+";
    }
    if (no == 1)
    {
        Particle = "B-";
    }
    if (no == 2)
    {
        Particle = "B0";
    }
    if (no == 3)
    {
        Particle = "Bc+";
    }
    if (no == 4)
    {
        Particle = "Bc-";
    }
    if (no == 5)
    {
        Particle = "Bs0";
    }
    if (no == 6)
    {
        Particle = "D+";
    }
    if (no == 7)
    {
        Particle = "D-";
    }
    if (no == 8)
    {
        Particle = "D0";
    }
    if (no == 9)
    {
        Particle = "Ds+";
    }
    if (no == 10)
    {
        Particle = "Ds-";
    }
    if (no == 11)
    {
        Particle = "GenericIon";
    }
    if (no == 12)
    {
        Particle = "He3";
    }
    if (no == 13)
    {
        Particle = "J/psi";
    }
    if (no == 14)
    {
        Particle = "N(1440)+";
    }
    if (no == 15)
    {
        Particle = "N(1440)0";
    }
    if (no == 16)
    {
        Particle = "N(1520)+";
    }
    if (no == 17)
    {
        Particle = "N(1520)0";
    }
    if (no == 18)
    {
        Particle = "N(1535)+";
    }
    if (no == 19)
    {
        Particle = "N(1535)0";
    }
    if (no == 20)
    {
        Particle = "N(1650)+";
    }
    if (no == 21)
    {
        Particle = "N(1650)0";
    }
    if (no == 22)
    {
        Particle = "N(1675)+";
    }
    if (no == 23)
    {
        Particle = "N(1675)0";
    }
    if (no == 24)
    {
        Particle = "N(1680)+";
    }
    if (no == 25)
    {
        Particle = "N(1680)0";
    }
    if (no == 26)
    {
        Particle = "N(1700)+";
    }
    if (no == 27)
    {
        Particle = "N(1700)0";
    }
    if (no == 28)
    {
        Particle = "N(1710)+";
    }
    if (no == 29)
    {
        Particle = "N(1710)0";
    }
    if (no == 30)
    {
        Particle = "N(1720)+";
    }
    if (no == 31)
    {
        Particle = "N(1720)0";
    }
    if (no == 32)
    {
        Particle = "N(1900)+";
    }
    if (no == 33)
    {
        Particle = "N(1900)0";
    }
    if (no == 34)
    {
        Particle = "N(1990)+";
    }
    if (no == 35)
    {
        Particle = "N(1990)0";
    }
    if (no == 36)
    {
        Particle = "N(2090)+";
    }
    if (no == 37)
    {
        Particle = "N(2090)0";
    }
    if (no == 38)
    {
        Particle = "N(2190)+";
    }
    if (no == 39)
    {
        Particle = "N(2190)0";
    }
    if (no == 40)
    {
        Particle = "N(2220)+";
    }
    if (no == 41)
    {
        Particle = "N(2220)0";
    }
    if (no == 42)
    {
        Particle = "N(2250)+";
    }
    if (no == 43)
    {
        Particle = "N(2250)0";
    }
    if (no == 44)
    {
        Particle = "Upsilon";
    }
    if (no == 45)
    {
        Particle = "a0(1450)+";
    }
    if (no == 46)
    {
        Particle = "a0(1450)-";
    }
    if (no == 47)
    {
        Particle = "a0(1450)0";
    }
    if (no == 48)
    {
        Particle = "a0(980)+";
    }
    if (no == 49)
    {
        Particle = "a0(980)-";
    }
    if (no == 50)
    {
        Particle = "a0(980)0";
    }
    if (no == 51)
    {
        Particle = "a1(1260)+";
    }
    if (no == 52)
    {
        Particle = "a1(1260)-";
    }
    if (no == 53)
    {
        Particle = "a1(1260)0";
    }
    if (no == 54)
    {
        Particle = "a2(1320)+";
    }
    if (no == 55)
    {
        Particle = "a2(1320)-";
    }
    if (no == 56)
    {
        Particle = "a2(1320)0";
    }
    if (no == 57)
    {
        Particle = "alpha";
    }
    if (no == 58)
    {
        Particle = "anti_B0";
    }
    if (no == 59)
    {
        Particle = "anti_Bs0";
    }
    if (no == 60)
    {
        Particle = "anti_D0";
    }
    if (no == 61)
    {
        Particle = "anti_He3";
    }
    if (no == 62)
    {
        Particle = "anti_N(1440)+";
    }
    if (no == 63)
    {
        Particle = "anti_N(1440)0";
    }
    if (no == 64)
    {
        Particle = "anti_N(1520)+";
    }
    if (no == 65)
    {
        Particle = "anti_N(1520)0";
    }
    if (no == 66)
    {
        Particle = "anti_N(1535)+";
    }
    if (no == 67)
    {
        Particle = "anti_N(1535)0";
    }
    if (no == 68)
    {
        Particle = "anti_N(1650)+";
    }
    if (no == 69)
    {
        Particle = "anti_N(1650)0";
    }
    if (no == 70)
    {
        Particle = "anti_N(1675)+";
    }
    if (no == 71)
    {
        Particle = "anti_N(1675)0";
    }
    if (no == 72)
    {
        Particle = "anti_N(1680)+";
    }
    if (no == 73)
    {
        Particle = "anti_N(1680)0";
    }
    if (no == 74)
    {
        Particle = "anti_N(1700)+";
    }
    if (no == 75)
    {
        Particle = "anti_N(1700)0";
    }
    if (no == 76)
    {
        Particle = "anti_N(1710)+";
    }
    if (no == 77)
    {
        Particle = "anti_N(1710)0";
    }
    if (no == 78)
    {
        Particle = "anti_N(1720)+";
    }
    if (no == 79)
    {
        Particle = "anti_N(1720)0";
    }
    if (no == 80)
    {
        Particle = "anti_N(1900)+";
    }
    if (no == 81)
    {
        Particle = "anti_N(1900)0";
    }
    if (no == 82)
    {
        Particle = "anti_N(1990)+";
    }
    if (no == 83)
    {
        Particle = "anti_N(1990)0";
    }
    if (no == 84)
    {
        Particle = "anti_N(2090)+";
    }
    if (no == 85)
    {
        Particle = "anti_N(2090)0";
    }
    if (no == 86)
    {
        Particle = "anti_N(2190)+";
    }
    if (no == 87)
    {
        Particle = "anti_N(2190)0";
    }
    if (no == 88)
    {
        Particle = "anti_N(2220)+";
    }
    if (no == 89)
    {
        Particle = "anti_N(2220)0";
    }
    if (no == 90)
    {
        Particle = "anti_N(2250)+";
    }
    if (no == 91)
    {
        Particle = "anti_N(2250)0";
    }
    if (no == 92)
    {
        Particle = "anti_alpha";
    }
    if (no == 93)
    {
        Particle = "anti_b_quark";
    }
    if (no == 94)
    {
        Particle = "anti_bb1_diquark";
    }
    if (no == 95)
    {
        Particle = "anti_bc0_diquark";
    }
    if (no == 96)
    {
        Particle = "anti_bc1_diquark";
    }
    if (no == 97)
    {
        Particle = "anti_bd0_diquark";
    }
    if (no == 98)
    {
        Particle = "anti_bd1_diquark";
    }
    if (no == 99)
    {
        Particle = "anti_bs0_diquark";
    }
    if (no == 100)
    {
        Particle = "anti_bs1_diquark";
    }
    if (no == 101)
    {
        Particle = "anti_bu0_diquark";
    }
    if (no == 102)
    {
        Particle = "anti_bu1_diquark";
    }
    if (no == 103)
    {
        Particle = "anti_c_quark";
    }
    if (no == 104)
    {
        Particle = "anti_cc1_diquark";
    }
    if (no == 105)
    {
        Particle = "anti_cd0_diquark";
    }
    if (no == 106)
    {
        Particle = "anti_cd1_diquark";
    }
    if (no == 107)
    {
        Particle = "anti_cs0_diquark";
    }
    if (no == 108)
    {
        Particle = "anti_cs1_diquark";
    }
    if (no == 109)
    {
        Particle = "anti_cu0_diquark";
    }
    if (no == 110)
    {
        Particle = "anti_cu1_diquark";
    }
    if (no == 111)
    {
        Particle = "anti_d_quark";
    }
    if (no == 112)
    {
        Particle = "anti_dd1_diquark";
    }
    if (no == 113)
    {
        Particle = "anti_delta(1600)+";
    }
    if (no == 114)
    {
        Particle = "anti_delta(1600)++";
    }
    if (no == 115)
    {
        Particle = "anti_delta(1600)-";
    }
    if (no == 116)
    {
        Particle = "anti_delta(1600)0";
    }
    if (no == 117)
    {
        Particle = "anti_delta(1620)+";
    }
    if (no == 118)
    {
        Particle = "anti_delta(1620)++";
    }
    if (no == 119)
    {
        Particle = "anti_delta(1620)-";
    }
    if (no == 120)
    {
        Particle = "anti_delta(1620)0";
    }
    if (no == 121)
    {
        Particle = "anti_delta(1700)+";
    }
    if (no == 122)
    {
        Particle = "anti_delta(1700)++";
    }
    if (no == 123)
    {
        Particle = "anti_delta(1700)-";
    }
    if (no == 124)
    {
        Particle = "anti_delta(1700)0";
    }
    if (no == 125)
    {
        Particle = "anti_delta(1900)+";
    }
    if (no == 126)
    {
        Particle = "anti_delta(1900)++";
    }
    if (no == 127)
    {
        Particle = "anti_delta(1900)-";
    }
    if (no == 128)
    {
        Particle = "anti_delta(1900)0";
    }
    if (no == 129)
    {
        Particle = "anti_delta(1905)+";
    }
    if (no == 130)
    {
        Particle = "anti_delta(1905)++";
    }
    if (no == 131)
    {
        Particle = "anti_delta(1905)-";
    }
    if (no == 132)
    {
        Particle = "anti_delta(1905)0";
    }
    if (no == 133)
    {
        Particle = "anti_delta(1910)+";
    }
    if (no == 134)
    {
        Particle = "anti_delta(1910)++";
    }
    if (no == 135)
    {
        Particle = "anti_delta(1910)-";
    }
    if (no == 136)
    {
        Particle = "anti_delta(1910)0";
    }
    if (no == 137)
    {
        Particle = "anti_delta(1920)+";
    }
    if (no == 138)
    {
        Particle = "anti_delta(1920)++";
    }
    if (no == 139)
    {
        Particle = "anti_delta(1920)-";
    }
    if (no == 140)
    {
        Particle = "anti_delta(1920)0";
    }
    if (no == 141)
    {
        Particle = "anti_delta(1930)+";
    }
    if (no == 142)
    {
        Particle = "anti_delta(1930)++";
    }
    if (no == 143)
    {
        Particle = "anti_delta(1930)-";
    }
    if (no == 144)
    {
        Particle = "anti_delta(1930)0";
    }
    if (no == 145)
    {
        Particle = "anti_delta(1950)+";
    }
    if (no == 146)
    {
        Particle = "anti_delta(1950)++";
    }
    if (no == 147)
    {
        Particle = "anti_delta(1950)-";
    }
    if (no == 148)
    {
        Particle = "anti_delta(1950)0";
    }
    if (no == 149)
    {
        Particle = "anti_delta+";
    }
    if (no == 150)
    {
        Particle = "anti_delta++";
    }
    if (no == 151)
    {
        Particle = "anti_delta-";
    }
    if (no == 152)
    {
        Particle = "anti_delta0";
    }
    if (no == 153)
    {
        Particle = "anti_deuteron";
    }
    if (no == 154)
    {
        Particle = "anti_doublehyperH4";
    }
    if (no == 155)
    {
        Particle = "anti_doublehyperdoubleneutron";
    }
    if (no == 156)
    {
        Particle = "anti_hyperH4";
    }
    if (no == 157)
    {
        Particle = "anti_hyperHe5";
    }
    if (no == 158)
    {
        Particle = "anti_hyperalpha";
    }
    if (no == 159)
    {
        Particle = "anti_hypertriton";
    }
    if (no == 160)
    {
        Particle = "anti_k(1460)0";
    }
    if (no == 161)
    {
        Particle = "anti_k0_star(1430)0";
    }
    if (no == 162)
    {
        Particle = "anti_k1(1270)0";
    }
    if (no == 163)
    {
        Particle = "anti_k1(1400)0";
    }
    if (no == 164)
    {
        Particle = "anti_k2(1770)0";
    }
    if (no == 165)
    {
        Particle = "anti_k2_star(1430)0";
    }
    if (no == 166)
    {
        Particle = "anti_k2_star(1980)0";
    }
    if (no == 167)
    {
        Particle = "anti_k3_star(1780)0";
    }
    if (no == 168)
    {
        Particle = "anti_k_star(1410)0";
    }
    if (no == 169)
    {
        Particle = "anti_k_star(1680)0";
    }
    if (no == 170)
    {
        Particle = "anti_k_star0";
    }
    if (no == 171)
    {
        Particle = "anti_kaon0";
    }
    if (no == 172)
    {
        Particle = "anti_lambda";
    }
    if (no == 173)
    {
        Particle = "anti_lambda(1405)";
    }
    if (no == 174)
    {
        Particle = "anti_lambda(1520)";
    }
    if (no == 175)
    {
        Particle = "anti_lambda(1600)";
    }
    if (no == 176)
    {
        Particle = "anti_lambda(1670)";
    }
    if (no == 177)
    {
        Particle = "anti_lambda(1690)";
    }
    if (no == 178)
    {
        Particle = "anti_lambda(1800)";
    }
    if (no == 179)
    {
        Particle = "anti_lambda(1810)";
    }
    if (no == 180)
    {
        Particle = "anti_lambda(1820)";
    }
    if (no == 181)
    {
        Particle = "anti_lambda(1830)";
    }
    if (no == 182)
    {
        Particle = "anti_lambda(1890)";
    }
    if (no == 183)
    {
        Particle = "anti_lambda(2100)";
    }
    if (no == 184)
    {
        Particle = "anti_lambda(2110)";
    }
    if (no == 185)
    {
        Particle = "anti_lambda_b";
    }
    if (no == 186)
    {
        Particle = "anti_lambda_c+";
    }
    if (no == 187)
    {
        Particle = "anti_neutron";
    }
    if (no == 188)
    {
        Particle = "anti_nu_e";
    }
    if (no == 189)
    {
        Particle = "anti_nu_mu";
    }
    if (no == 190)
    {
        Particle = "anti_nu_tau";
    }
    if (no == 191)
    {
        Particle = "anti_omega-";
    }
    if (no == 192)
    {
        Particle = "anti_omega_b-";
    }
    if (no == 193)
    {
        Particle = "anti_omega_c0";
    }
    if (no == 194)
    {
        Particle = "anti_proton";
    }
    if (no == 195)
    {
        Particle = "anti_s_quark";
    }
    if (no == 196)
    {
        Particle = "anti_sd0_diquark";
    }
    if (no == 197)
    {
        Particle = "anti_sd1_diquark";
    }
    if (no == 198)
    {
        Particle = "anti_sigma(1385)+";
    }
    if (no == 199)
    {
        Particle = "anti_sigma(1385)-";
    }
    if (no == 200)
    {
        Particle = "anti_sigma(1385)0";
    }
    if (no == 201)
    {
        Particle = "anti_sigma(1660)+";
    }
    if (no == 202)
    {
        Particle = "anti_sigma(1660)-";
    }
    if (no == 203)
    {
        Particle = "anti_sigma(1660)0";
    }
    if (no == 204)
    {
        Particle = "anti_sigma(1670)+";
    }
    if (no == 205)
    {
        Particle = "anti_sigma(1670)-";
    }
    if (no == 206)
    {
        Particle = "anti_sigma(1670)0";
    }
    if (no == 207)
    {
        Particle = "anti_sigma(1750)+";
    }
    if (no == 208)
    {
        Particle = "anti_sigma(1750)-";
    }
    if (no == 209)
    {
        Particle = "anti_sigma(1750)0";
    }
    if (no == 210)
    {
        Particle = "anti_sigma(1775)+";
    }
    if (no == 211)
    {
        Particle = "anti_sigma(1775)-";
    }
    if (no == 212)
    {
        Particle = "anti_sigma(1775)0";
    }
    if (no == 213)
    {
        Particle = "anti_sigma(1915)+";
    }
    if (no == 214)
    {
        Particle = "anti_sigma(1915)-";
    }
    if (no == 215)
    {
        Particle = "anti_sigma(1915)0";
    }
    if (no == 216)
    {
        Particle = "anti_sigma(1940)+";
    }
    if (no == 217)
    {
        Particle = "anti_sigma(1940)-";
    }
    if (no == 218)
    {
        Particle = "anti_sigma(1940)0";
    }
    if (no == 219)
    {
        Particle = "anti_sigma(2030)+";
    }
    if (no == 220)
    {
        Particle = "anti_sigma(2030)-";
    }
    if (no == 221)
    {
        Particle = "anti_sigma(2030)0";
    }
    if (no == 222)
    {
        Particle = "anti_sigma+";
    }
    if (no == 223)
    {
        Particle = "anti_sigma-";
    }
    if (no == 224)
    {
        Particle = "anti_sigma0";
    }
    if (no == 225)
    {
        Particle = "anti_sigma_b+";
    }
    if (no == 226)
    {
        Particle = "anti_sigma_b-";
    }
    if (no == 227)
    {
        Particle = "anti_sigma_b0";
    }
    if (no == 228)
    {
        Particle = "anti_sigma_c+";
    }
    if (no == 229)
    {
        Particle = "anti_sigma_c++";
    }
    if (no == 230)
    {
        Particle = "anti_sigma_c0";
    }
    if (no == 231)
    {
        Particle = "anti_ss1_diquark";
    }
    if (no == 232)
    {
        Particle = "anti_su0_diquark";
    }
    if (no == 233)
    {
        Particle = "anti_su1_diquark";
    }
    if (no == 234)
    {
        Particle = "anti_t_quark";
    }
    if (no == 235)
    {
        Particle = "anti_triton";
    }
    if (no == 236)
    {
        Particle = "anti_u_quark";
    }
    if (no == 237)
    {
        Particle = "anti_ud0_diquark";
    }
    if (no == 238)
    {
        Particle = "anti_ud1_diquark";
    }
    if (no == 239)
    {
        Particle = "anti_uu1_diquark";
    }
    if (no == 240)
    {
        Particle = "anti_xi(1530)-";
    }
    if (no == 241)
    {
        Particle = "anti_xi(1530)0";
    }
    if (no == 242)
    {
        Particle = "anti_xi(1690)-";
    }
    if (no == 243)
    {
        Particle = "anti_xi(1690)0";
    }
    if (no == 244)
    {
        Particle = "anti_xi(1820)-";
    }
    if (no == 245)
    {
        Particle = "anti_xi(1820)0";
    }
    if (no == 246)
    {
        Particle = "anti_xi(1950)-";
    }
    if (no == 247)
    {
        Particle = "anti_xi(1950)0";
    }
    if (no == 248)
    {
        Particle = "anti_xi(2030)-";
    }
    if (no == 249)
    {
        Particle = "anti_xi(2030)0";
    }
    if (no == 250)
    {
        Particle = "anti_xi-";
    }
    if (no == 251)
    {
        Particle = "anti_xi0";
    }
    if (no == 252)
    {
        Particle = "anti_xi_b-";
    }
    if (no == 253)
    {
        Particle = "anti_xi_b0";
    }
    if (no == 254)
    {
        Particle = "anti_xi_c+";
    }
    if (no == 255)
    {
        Particle = "anti_xi_c0";
    }
    if (no == 256)
    {
        Particle = "b1(1235)+";
    }
    if (no == 257)
    {
        Particle = "b1(1235)-";
    }
    if (no == 258)
    {
        Particle = "b1(1235)0";
    }
    if (no == 259)
    {
        Particle = "b_quark";
    }
    if (no == 260)
    {
        Particle = "bb1_diquark";
    }
    if (no == 261)
    {
        Particle = "bc0_diquark";
    }
    if (no == 262)
    {
        Particle = "bc1_diquark";
    }
    if (no == 263)
    {
        Particle = "bd0_diquark";
    }
    if (no == 264)
    {
        Particle = "bd1_diquark";
    }
    if (no == 265)
    {
        Particle = "bs0_diquark";
    }
    if (no == 266)
    {
        Particle = "bs1_diquark";
    }
    if (no == 267)
    {
        Particle = "bu0_diquark";
    }
    if (no == 268)
    {
        Particle = "bu1_diquark";
    }
    if (no == 269)
    {
        Particle = "c_quark";
    }
    if (no == 270)
    {
        Particle = "cc1_diquark";
    }
    if (no == 271)
    {
        Particle = "cd0_diquark";
    }
    if (no == 272)
    {
        Particle = "cd1_diquark";
    }
    if (no == 273)
    {
        Particle = "chargedgeantino";
    }
    if (no == 274)
    {
        Particle = "cs0_diquark";
    }
    if (no == 275)
    {
        Particle = "cs1_diquark";
    }
    if (no == 276)
    {
        Particle = "cu0_diquark";
    }
    if (no == 277)
    {
        Particle = "cu1_diquark";
    }
    if (no == 278)
    {
        Particle = "d_quark";
    }
    if (no == 279)
    {
        Particle = "dd1_diquark";
    }
    if (no == 280)
    {
        Particle = "delta(1600)+";
    }
    if (no == 281)
    {
        Particle = "delta(1600)++";
    }
    if (no == 282)
    {
        Particle = "delta(1600)-";
    }
    if (no == 283)
    {
        Particle = "delta(1600)0";
    }
    if (no == 284)
    {
        Particle = "delta(1620)+";
    }
    if (no == 285)
    {
        Particle = "delta(1620)++";
    }
    if (no == 286)
    {
        Particle = "delta(1620)-";
    }
    if (no == 287)
    {
        Particle = "delta(1620)0";
    }
    if (no == 288)
    {
        Particle = "delta(1700)+";
    }
    if (no == 289)
    {
        Particle = "delta(1700)++";
    }
    if (no == 290)
    {
        Particle = "delta(1700)-";
    }
    if (no == 291)
    {
        Particle = "delta(1700)0";
    }
    if (no == 292)
    {
        Particle = "delta(1900)+";
    }
    if (no == 293)
    {
        Particle = "delta(1900)++";
    }
    if (no == 294)
    {
        Particle = "delta(1900)-";
    }
    if (no == 295)
    {
        Particle = "delta(1900)0";
    }
    if (no == 296)
    {
        Particle = "delta(1905)+";
    }
    if (no == 297)
    {
        Particle = "delta(1905)++";
    }
    if (no == 298)
    {
        Particle = "delta(1905)-";
    }
    if (no == 299)
    {
        Particle = "delta(1905)0";
    }
    if (no == 300)
    {
        Particle = "delta(1910)+";
    }
    if (no == 301)
    {
        Particle = "delta(1910)++";
    }
    if (no == 302)
    {
        Particle = "delta(1910)-";
    }
    if (no == 303)
    {
        Particle = "delta(1910)0";
    }
    if (no == 304)
    {
        Particle = "delta(1920)+";
    }
    if (no == 305)
    {
        Particle = "delta(1920)++";
    }
    if (no == 306)
    {
        Particle = "delta(1920)-";
    }
    if (no == 307)
    {
        Particle = "delta(1920)0";
    }
    if (no == 308)
    {
        Particle = "delta(1930)+";
    }
    if (no == 309)
    {
        Particle = "delta(1930)++";
    }
    if (no == 310)
    {
        Particle = "delta(1930)-";
    }
    if (no == 311)
    {
        Particle = "delta(1930)0";
    }
    if (no == 312)
    {
        Particle = "delta(1950)+";
    }
    if (no == 313)
    {
        Particle = "delta(1950)++";
    }
    if (no == 314)
    {
        Particle = "delta(1950)-";
    }
    if (no == 315)
    {
        Particle = "delta(1950)0";
    }
    if (no == 316)
    {
        Particle = "delta+";
    }
    if (no == 317)
    {
        Particle = "delta++";
    }
    if (no == 318)
    {
        Particle = "delta-";
    }
    if (no == 319)
    {
        Particle = "delta0";
    }
    if (no == 320)
    {
        Particle = "deuteron";
    }
    if (no == 321)
    {
        Particle = "doublehyperH4";
    }
    if (no == 322)
    {
        Particle = "doublehyperdoubleneutron";
    }
    if (no == 323)
    {
        Particle = "e+";
    }
    if (no == 324)
    {
        Particle = "e-";
    }
    if (no == 325)
    {
        Particle = "eta";
    }
    if (no == 326)
    {
        Particle = "eta(1295)";
    }
    if (no == 327)
    {
        Particle = "eta(1405)";
    }
    if (no == 328)
    {
        Particle = "eta(1475)";
    }
    if (no == 329)
    {
        Particle = "eta2(1645)";
    }
    if (no == 330)
    {
        Particle = "eta2(1870)";
    }
    if (no == 331)
    {
        Particle = "eta_prime";
    }
    if (no == 332)
    {
        Particle = "etac";
    }
    if (no == 333)
    {
        Particle = "f0(1370)";
    }
    if (no == 334)
    {
        Particle = "f0(1500)";
    }
    if (no == 335)
    {
        Particle = "f0(1710)";
    }
    if (no == 336)
    {
        Particle = "f0(500)";
    }
    if (no == 337)
    {
        Particle = "f0(980)";
    }
    if (no == 338)
    {
        Particle = "f1(1285)";
    }
    if (no == 339)
    {
        Particle = "f1(1420)";
    }
    if (no == 340)
    {
        Particle = "f2(1270)";
    }
    if (no == 341)
    {
        Particle = "f2(1810)";
    }
    if (no == 342)
    {
        Particle = "f2(2010)";
    }
    if (no == 343)
    {
        Particle = "f2_prime(1525)";
    }
    if (no == 344)
    {
        Particle = "gamma";
    }
    if (no == 345)
    {
        Particle = "geantino";
    }
    if (no == 346)
    {
        Particle = "gluon";
    }
    if (no == 347)
    {
        Particle = "h1(1170)";
    }
    if (no == 348)
    {
        Particle = "h1(1380)";
    }
    if (no == 349)
    {
        Particle = "hyperH4";
    }
    if (no == 350)
    {
        Particle = "hyperHe5";
    }
    if (no == 351)
    {
        Particle = "hyperalpha";
    }
    if (no == 352)
    {
        Particle = "hypertriton";
    }
    if (no == 353)
    {
        Particle = "k(1460)+";
    }
    if (no == 354)
    {
        Particle = "k(1460)-";
    }
    if (no == 355)
    {
        Particle = "k(1460)0";
    }
    if (no == 356)
    {
        Particle = "k0_star(1430)+";
    }
    if (no == 357)
    {
        Particle = "k0_star(1430)-";
    }
    if (no == 358)
    {
        Particle = "k0_star(1430)0";
    }
    if (no == 359)
    {
        Particle = "k1(1270)+";
    }
    if (no == 360)
    {
        Particle = "k1(1270)-";
    }
    if (no == 361)
    {
        Particle = "k1(1270)0";
    }
    if (no == 362)
    {
        Particle = "k1(1400)+";
    }
    if (no == 363)
    {
        Particle = "k1(1400)-";
    }
    if (no == 364)
    {
        Particle = "k1(1400)0";
    }
    if (no == 365)
    {
        Particle = "k2(1770)+";
    }
    if (no == 366)
    {
        Particle = "k2(1770)-";
    }
    if (no == 367)
    {
        Particle = "k2(1770)0";
    }
    if (no == 368)
    {
        Particle = "k2_star(1430)+";
    }
    if (no == 369)
    {
        Particle = "k2_star(1430)-";
    }
    if (no == 370)
    {
        Particle = "k2_star(1430)0";
    }
    if (no == 371)
    {
        Particle = "k2_star(1980)+";
    }
    if (no == 372)
    {
        Particle = "k2_star(1980)-";
    }
    if (no == 373)
    {
        Particle = "k2_star(1980)0";
    }
    if (no == 374)
    {
        Particle = "k3_star(1780)+";
    }
    if (no == 375)
    {
        Particle = "k3_star(1780)-";
    }
    if (no == 376)
    {
        Particle = "k3_star(1780)0";
    }
    if (no == 377)
    {
        Particle = "k_star(1410)+";
    }
    if (no == 378)
    {
        Particle = "k_star(1410)-";
    }
    if (no == 379)
    {
        Particle = "k_star(1410)0";
    }
    if (no == 380)
    {
        Particle = "k_star(1680)+";
    }
    if (no == 381)
    {
        Particle = "k_star(1680)-";
    }
    if (no == 382)
    {
        Particle = "k_star(1680)0";
    }
    if (no == 383)
    {
        Particle = "k_star+";
    }
    if (no == 384)
    {
        Particle = "k_star-";
    }
    if (no == 385)
    {
        Particle = "k_star0";
    }
    if (no == 386)
    {
        Particle = "kaon+";
    }
    if (no == 387)
    {
        Particle = "kaon-";
    }
    if (no == 388)
    {
        Particle = "kaon0";
    }
    if (no == 389)
    {
        Particle = "kaon0L";
    }
    if (no == 390)
    {
        Particle = "kaon0S";
    }
    if (no == 391)
    {
        Particle = "lambda";
    }
    if (no == 392)
    {
        Particle = "lambda(1405)";
    }
    if (no == 393)
    {
        Particle = "lambda(1520)";
    }
    if (no == 394)
    {
        Particle = "lambda(1600)";
    }
    if (no == 395)
    {
        Particle = "lambda(1670)";
    }
    if (no == 396)
    {
        Particle = "lambda(1690)";
    }
    if (no == 397)
    {
        Particle = "lambda(1800)";
    }
    if (no == 398)
    {
        Particle = "lambda(1810)";
    }
    if (no == 399)
    {
        Particle = "lambda(1820)";
    }
    if (no == 400)
    {
        Particle = "lambda(1830)";
    }
    if (no == 401)
    {
        Particle = "lambda(1890)";
    }
    if (no == 402)
    {
        Particle = "lambda(2100)";
    }
    if (no == 403)
    {
        Particle = "lambda(2110)";
    }
    if (no == 404)
    {
        Particle = "lambda_b";
    }
    if (no == 405)
    {
        Particle = "lambda_c+";
    }
    if (no == 406)
    {
        Particle = "mu+";
    }
    if (no == 407)
    {
        Particle = "mu-";
    }
    if (no == 408)
    {
        Particle = "neutron";
    }
    if (no == 409)
    {
        Particle = "nu_e";
    }
    if (no == 410)
    {
        Particle = "nu_mu";
    }
    if (no == 411)
    {
        Particle = "nu_tau";
    }
    if (no == 412)
    {
        Particle = "omega";
    }
    if (no == 413)
    {
        Particle = "omega(1420)";
    }
    if (no == 414)
    {
        Particle = "omega(1650)";
    }
    if (no == 415)
    {
        Particle = "omega-";
    }
    if (no == 416)
    {
        Particle = "omega3(1670)";
    }
    if (no == 417)
    {
        Particle = "omega_b-";
    }
    if (no == 418)
    {
        Particle = "omega_c0";
    }
    if (no == 419)
    {
        Particle = "opticalphoton";
    }
    if (no == 420)
    {
        Particle = "phi";
    }
    if (no == 421)
    {
        Particle = "phi(1680)";
    }
    if (no == 422)
    {
        Particle = "phi3(1850)";
    }
    if (no == 423)
    {
        Particle = "pi(1300)+";
    }
    if (no == 424)
    {
        Particle = "pi(1300)-";
    }
    if (no == 425)
    {
        Particle = "pi(1300)0";
    }
    if (no == 426)
    {
        Particle = "pi+";
    }
    if (no == 427)
    {
        Particle = "pi-";
    }
    if (no == 428)
    {
        Particle = "pi0";
    }
    if (no == 429)
    {
        Particle = "pi2(1670)+";
    }
    if (no == 430)
    {
        Particle = "pi2(1670)-";
    }
    if (no == 431)
    {
        Particle = "pi2(1670)0";
    }
    if (no == 432)
    {
        Particle = "proton";
    }
    if (no == 433)
    {
        Particle = "rho(1450)+";
    }
    if (no == 434)
    {
        Particle = "rho(1450)-";
    }
    if (no == 435)
    {
        Particle = "rho(1450)0";
    }
    if (no == 436)
    {
        Particle = "rho(1700)+";
    }
    if (no == 437)
    {
        Particle = "rho(1700)-";
    }
    if (no == 438)
    {
        Particle = "rho(1700)0";
    }
    if (no == 439)
    {
        Particle = "rho+";
    }
    if (no == 440)
    {
        Particle = "rho-";
    }
    if (no == 441)
    {
        Particle = "rho0";
    }
    if (no == 442)
    {
        Particle = "rho3(1690)+";
    }
    if (no == 443)
    {
        Particle = "rho3(1690)-";
    }
    if (no == 444)
    {
        Particle = "rho3(1690)0";
    }
    if (no == 445)
    {
        Particle = "s_quark";
    }
    if (no == 446)
    {
        Particle = "sd0_diquark";
    }
    if (no == 447)
    {
        Particle = "sd1_diquark";
    }
    if (no == 448)
    {
        Particle = "sigma(1385)+";
    }
    if (no == 449)
    {
        Particle = "sigma(1385)-";
    }
    if (no == 450)
    {
        Particle = "sigma(1385)0";
    }
    if (no == 451)
    {
        Particle = "sigma(1660)+";
    }
    if (no == 452)
    {
        Particle = "sigma(1660)-";
    }
    if (no == 453)
    {
        Particle = "sigma(1660)0";
    }
    if (no == 454)
    {
        Particle = "sigma(1670)+";
    }
    if (no == 455)
    {
        Particle = "sigma(1670)-";
    }
    if (no == 456)
    {
        Particle = "sigma(1670)0";
    }
    if (no == 457)
    {
        Particle = "sigma(1750)+";
    }
    if (no == 458)
    {
        Particle = "sigma(1750)-";
    }
    if (no == 459)
    {
        Particle = "sigma(1750)0";
    }
    if (no == 460)
    {
        Particle = "sigma(1775)+";
    }
    if (no == 461)
    {
        Particle = "sigma(1775)-";
    }
    if (no == 462)
    {
        Particle = "sigma(1775)0";
    }
    if (no == 463)
    {
        Particle = "sigma(1915)+";
    }
    if (no == 464)
    {
        Particle = "sigma(1915)-";
    }
    if (no == 465)
    {
        Particle = "sigma(1915)0";
    }
    if (no == 466)
    {
        Particle = "sigma(1940)+";
    }
    if (no == 467)
    {
        Particle = "sigma(1940)-";
    }
    if (no == 468)
    {
        Particle = "sigma(1940)0";
    }
    if (no == 469)
    {
        Particle = "sigma(2030)+";
    }
    if (no == 470)
    {
        Particle = "sigma(2030)-";
    }
    if (no == 471)
    {
        Particle = "sigma(2030)0";
    }
    if (no == 472)
    {
        Particle = "sigma+";
    }
    if (no == 473)
    {
        Particle = "sigma-";
    }
    if (no == 474)
    {
        Particle = "sigma0";
    }
    if (no == 475)
    {
        Particle = "sigma_b+";
    }
    if (no == 476)
    {
        Particle = "sigma_b-";
    }
    if (no == 477)
    {
        Particle = "sigma_b0";
    }
    if (no == 478)
    {
        Particle = "sigma_c+";
    }
    if (no == 479)
    {
        Particle = "sigma_c++";
    }
    if (no == 480)
    {
        Particle = "sigma_c0";
    }
    if (no == 481)
    {
        Particle = "ss1_diquark";
    }
    if (no == 482)
    {
        Particle = "su0_diquark";
    }
    if (no == 483)
    {
        Particle = "su1_diquark";
    }
    if (no == 484)
    {
        Particle = "t_quark";
    }
    if (no == 485)
    {
        Particle = "tau+";
    }
    if (no == 486)
    {
        Particle = "tau-";
    }
    if (no == 487)
    {
        Particle = "triton";
    }
    if (no == 488)
    {
        Particle = "u_quark";
    }
    if (no == 489)
    {
        Particle = "ud0_diquark";
    }
    if (no == 490)
    {
        Particle = "ud1_diquark";
    }
    if (no == 491)
    {
        Particle = "uu1_diquark";
    }
    if (no == 492)
    {
        Particle = "xi(1530)-";
    }
    if (no == 493)
    {
        Particle = "xi(1530)0";
    }
    if (no == 494)
    {
        Particle = "xi(1690)-";
    }
    if (no == 495)
    {
        Particle = "xi(1690)0";
    }
    if (no == 496)
    {
        Particle = "xi(1820)-";
    }
    if (no == 497)
    {
        Particle = "xi(1820)0";
    }
    if (no == 498)
    {
        Particle = "xi(1950)-";
    }
    if (no == 499)
    {
        Particle = "xi(1950)0";
    }
    if (no == 500)
    {
        Particle = "xi(2030)-";
    }
    if (no == 501)
    {
        Particle = "xi(2030)0";
    }
    if (no == 502)
    {
        Particle = "xi-";
    }
    if (no == 503)
    {
        Particle = "xi0";
    }
    if (no == 504)
    {
        Particle = "xi_b-";
    }
    if (no == 505)
    {
        Particle = "xi_b0";
    }
    if (no == 506)
    {
        Particle = "xi_c+";
    }
    if (no == 507)
    {
        Particle = "xi_c0";
    }
}