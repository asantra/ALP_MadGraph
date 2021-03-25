#define drawGridRootFiles_cxx
#include "drawGridRootFiles.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <chrono>  // for high_resolution_clock
#include <map>
#include <string>
#include "TLorentzVector.h"
#include <fstream>
#include "TMath.h"
#include "TF1.h"
#include "TRandom3.h"
#include <fstream>
#include <typeinfo>
#include <sys/stat.h> 
#include <sys/types.h>
#include <sstream>
#include "math.h"

using namespace std;

/// global variable regarding the output filename which depends on the process

//TString outputFileSuffix = "photons_from_gbeam_JETI40_025fs_6500nm";
//TString outputFileSuffix = "photons_from_gbeam_phase2_120fs_10000nm";
TString outputFileSuffix = "photons_from_gbeam_phase2_120fs_10000nm_tungsten";

double openingAngle(TLorentzVector a, TLorentzVector b){
    double angle = a.Angle(b.Vect());
    return angle;
}

double dPhiCalc(double phiLead, double phiTrail){
  double dphi = fabs(phiLead - phiTrail);
  if(dphi > TMath::Pi()) dphi = TMath::Pi()*2. - dphi;
  return dphi;
}


double dThetaCalc(double thetaLead, double thetaTrail){
  double dtheta = (thetaLead - thetaTrail);
  return dtheta;
}

double thetaCal(double eta){
    double theta = 2*atan(exp(-eta));
    return theta*180/TMath::Pi();
}

double getTau(std::string mALP, std::string version){
      double tauA              = -999.0;
      double width             = -999.0;
      std::ifstream myfile("/storage/agrp/arkas/ALPFiles/LHEFiles/xSecValuesFromRunTagFiles/mass"+mALP+"GeV/run_01_tag_1_banner_"+mALP+"GeV_beam1.0GeV_v"+version+".txt", ios::in);
      
      //std::cout << "taufilename: /storage/agrp/arkas/ALPFiles/LHEFiles/xSecValuesFromRunTagFiles/run_01_tag_1_banner_"+mALP+"GeV_beam1.0GeV_v"+version+".txt" << std::endl;
      
      std::string str;
      string decayLine, pdgId; // for storing each word
      while (std::getline(myfile, str)){   
        if (str.find("DECAY  9000005") != std::string::npos){
            istringstream ss(str);
            ss >> decayLine >> pdgId >> width;
            break;
        }
      }
      if(width!=-999.0)
          tauA = 1.0/width;
      else{
          return -1;
      }
      return tauA;
}


// //// get the xsec from the text file
double getXsecFromFile(std::string mALP, double beamEnergy, std::string version){
    
    std::string beamEnergyStr = Form("%.1f",beamEnergy);
    //std::cout << "xsecfilename: /storage/agrp/arkas/ALPFiles/LHEFiles/xSecValuesFromRunTagFiles/run_01_tag_1_banner_"+mALP+"GeV_beam"+beamEnergyStr+"GeV_v"+version+".txt" << std::endl;
    std::ifstream myfile("/storage/agrp/arkas/ALPFiles/LHEFiles/xSecValuesFromRunTagFiles/mass"+mALP+"GeV/run_01_tag_1_banner_"+mALP+"GeV_beam"+beamEnergyStr+"GeV_v"+version+".txt", ios::in);
    
    std::string str;
    std::string hash, word1, word2, pb, colon;
    double number=-99999999999.0; // for storing each word
    while (std::getline(myfile, str))
    {   
        if (str.find("Integrated weight") != std::string::npos){
            istringstream ss(str);
            ss >> hash >> word1 >> word2 >> pb >> colon >> number;
            break;
        }
    }
    return number;
}

/// this is weight factor coming from the signal samples made by Anthony
double getWeightFactor(double mALP){
    double factor = 0.0;
    /// getting the xsec and nPhoton weights from a text file.
    string line;
    ifstream myfile("/storage/agrp/arkas/ALPFiles/GridAnalyzerDirectory/outputPhotonNumbers_"+outputFileSuffix+".txt");
    string dummyLine;
    //std::getline(myfile, dummyLine); // go over the first line which is not numbers
    double energy, nPhoton, dwde;
    
    /// loop over the numbered lines in the text file
    if (myfile.is_open()){
        while (myfile >> energy >> nPhoton >> dwde){
           if(energy > mALP)factor += dwde;
        }
        myfile.close();
    }
    return factor;
}

//// xsec and width factor coming from lambda
double getLambdaFactor(std::string oneOverLambda){
    double lambda       = 1./std::stof(oneOverLambda);
    double fInMG        = 4*lambda; /// in MG system, the coupling f is 4*lambda
    double fUsed        = 4*pow(10,4);
    double lambdaFactor = pow(fUsed,2)/pow(fInMG,2);
    return lambdaFactor;
}

//// the main code snippet
void drawGridRootFiles::Loop(std::string mass, std::string version, std::string cutPhoton, std::string oneOverLambda)
{
    /// this is to generate random numbers with proper seed
    bool debug = false;
//   In a ROOT session, you can do:
//      root> .L drawGridRootFiles.C
//      root> drawGridRootFiles t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
    if (fChain == 0) return;
    // start the clock
    auto start = std::chrono::steady_clock::now();
    
    /// create the histogram of distribution
//     TFile *fOut                         = new TFile("combinedOutputFileOneEvent.root", "RECREATE");
    TString photonNameFile= "";
    double photonThreshold=-999.0;
    if(cutPhoton=="0.5"){
        photonNameFile  = "_0.5GeVPhotonCut";
        photonThreshold = 0.5;
    }
    else{
        photonNameFile  = "_0.0GeVPhotonCut";
        photonThreshold = 0.0;
    }
    
    /// open the directory for the acceptance files; if the directory doesn't exist, create it
    string dirname       = "/storage/agrp/arkas/ALPFiles/GridAnalyzerDirectory/acceptanceFiles";
    TString dirnameT     = (TString)dirname;
    int status           = mkdir(dirname.c_str(),0777);
    
    string dirnameMass   = "/storage/agrp/arkas/ALPFiles/GridAnalyzerDirectory/acceptanceFiles/mass"+mass+"GeV";
    TString dirnameMassT = (TString)dirnameMass;
    int status2          = mkdir(dirnameMass.c_str(),0777);
    
    /// the text file containing all the numbers
    ofstream acceptanceFiles(dirnameT+"/PhotonAcceptanceTextFiles_"+outputFileSuffix+"_ALPAllMass"+photonNameFile+"_V"+version.c_str()+".txt", ios::app);
    
    /// the output root file containing all the distributions
    TFile *fOut                         = new TFile(dirnameMassT+"/PhotonAcceptanceFiles_"+outputFileSuffix+"_ALPMass"+mass.c_str()+"GeV_oneOverLambda"+oneOverLambda.c_str()+photonNameFile+"_V"+version.c_str()+".root", "RECREATE");
    fOut->cd();
    
    ///1D histogram
    TH1F *h_IncomingPhotonE                   = new TH1F("h_IncomingPhotonE", "Incoming photon energy; E [GeV]; Events", 210, 0, 10.5);
    TH1F *h_IncomingPhotonPhi                 = new TH1F("h_IncomingPhotonPhi", "Incoming photon #phi; #phi; Events", 64, -3.2, 3.2);
    TH1F *h_IncomingPhotonEta                 = new TH1F("h_IncomingPhotonEta", "Incoming photon #eta; #eta; Events", 3000, -1500.0, 1500.0);
    TH1F *h_OutgoingPhotonE_NoCut             = new TH1F("h_OutgoingPhotonE_NoCut", "Outgoing photon energy, without cuts; E [GeV]; Events", 500, 0, 20);
    TH1F *h_OutgoingPhotonE                   = new TH1F("h_OutgoingPhotonE", "Outgoing photon energy; E [GeV]; Events", 500, 0, 20);
    TH1F *h_OutgoingPhotonOpeningAngle_NoCut  = new TH1F("h_OutgoingPhotonOpeningAngle_NoCut", "Outgoing photon opening angle, no cut; Angle [radian]; Events", 200, 0, 3.4);
    TH1F *h_OutgoingPhotonOpeningAngle        = new TH1F("h_OutgoingPhotonOpeningAngle", "Outgoing photon opening angle; Angle [radian]; Events", 200, 0, 3.4);
    TH1F *h_OutgoingPhotonOpeningTheta        = new TH1F("h_OutgoingPhotonOpeningTheta", "Outgoing photon opening theta; #Theta [radian]; Events", 400, -3.4, 3.4);
    TH1F *h_DistanceInDetector                = new TH1F("h_DistanceInDetector", "Distance of 2 photons on detector; d [m]; Events", 1000, 0, 100);
    TH1F *h_DistanceInDetectorLessBins        = new TH1F("h_DistanceInDetectorLessBins", "Distance of 2 photons on detector; d [m]; Events", 100, 0, 5.0);
    TH1F *h_rALP                              = new TH1F("h_rALP", "h_rALP",600,-2,23);
    TH1F *h_probALP                           = new TH1F("h_probALP", "h_rALP",200,0,1);
    TH1F *h_RatioLA                           = new TH1F("h_RatioLA", "h_RatioLA",100,0,2);
    TH1F *h_ALPDecayTime                      = new TH1F("h_ALPDecayTime", "decay time of ALPs assuming speed of light, E > 0.5 GeV; time [ns]; Events", 100, 0, 25);
    TH1F *h_ALPArrivalTime                    = new TH1F("h_ALPArrivalTime", "arrival time of ALPs assuming speed of light, E > 0.5 GeV; time [ns]; Events", 100, 0, 25);
    TH1F *h_ALPDecayTimeEGt0p1                = new TH1F("h_ALPDecayTimeEGt0p1", "decay time of ALPs assuming speed of light, E > 0.1 GeV; time [ns]; Events", 100, 0, 25);
    TH1F *h_ALPArrivalTimeEGt0p1              = new TH1F("h_ALPArrivalTimeEGt0p1", "arrival time of ALPs assuming speed of light, E > 0.1 GeV; time [ns]; Events", 100, 0, 25);
    
    // TH2
    TH2F *h_PhotonY_PhotonX                   = new TH2F("h_PhotonY_PhotonX", "photon Y vs photon X, photon E > 0.0 GeV; X [m]; Y [m]", 200, -1.0, 1.0, 200, -1.0, 1.0);
    TH2F *h_PhotonY_PhotonX_E0p5GeV           = new TH2F("h_PhotonY_PhotonX_E0p5GeV", "photon Y vs photon X, photon E > 0.5 GeV; X [m]; Y [m]", 200, -1.0, 1.0, 200, -1.0, 1.0);
    TH2F *h_ALPY_ALPX                         = new TH2F("h_ALPY_ALPX", "alp Y vs alp X, decayed photon E > 0.0 GeV; X [m]; Y [m]", 500, -0.35, 0.35, 500, -0.35, 0.35);
    TH2F *h_ALPY_ALPX_E0p5GeV                 = new TH2F("h_ALPY_ALPX_E0p5GeV", "alp Y vs alp X, decayed photon E > 0.5 GeV; X [m]; Y [m]", 500, -0.35, 0.35, 500, -0.35, 0.35);
    TH2F *h_ALPRxy_Z_E0p5GeV                  = new TH2F("h_ALPRxy_Z_E0p5GeV", "alp R(x,y) vs alp decay Z, decayed photon E > 0.5 GeV; Z [m]; R_{xy} [m]", 300, 0.5, 3.5, 500, 0.0, 5.0);
    
    
    
    Long64_t nentries = fChain->GetEntriesFast();
    /// map to store the weightInput from the signal
    map<double, vector<double> > weightInput;
    
    
    /// getting the xsec and nPhoton weights from a text file.
    string line;
    ifstream myfile("/storage/agrp/arkas/ALPFiles/GridAnalyzerDirectory/outputPhotonNumbers_"+outputFileSuffix+".txt");
    string dummyLine;
    //std::getline(myfile, dummyLine); // go over the first line which is not numbers
    double energy, nPhoton, dwde;
    
    /// loop over the numbered lines in the text file
    if (myfile.is_open()){
        while (myfile >> energy >> nPhoton >> dwde){
            vector<double> photonDetails;
            photonDetails.push_back(nPhoton);
            photonDetails.push_back(dwde);
            weightInput.insert(std::pair<double, vector<double> >(energy, photonDetails));
        }
        myfile.close();
    }
    
    double mALP                = std::stof(mass); /// taking different ALP mass values
    double oneOverLambdaDouble = std::stof(oneOverLambda);
    //double tauAMG              = getTau(mass, "8")*getLambdaFactor(oneOverLambda); /// modify the tauA for differnt lambda
    double tauATheory          = 1./(pow(mALP,3)*pow(oneOverLambdaDouble,2)/(64*TMath::Pi())); // width = m^3/(64*pi*lambda^2) theoretically
    double width               = pow(mALP,3)*pow(oneOverLambdaDouble,2)/(64*TMath::Pi());
    double tauA                = tauATheory;
    double ctauA               = tauA*(0.1975)*pow(10,-15); // c is 1 in natural units, the additional factor is the conversion factor from GeV^{-1} to meters
    std::cout << "For mass " << mass << " GeV and version " << version << " the c*tauA is " << ctauA << std::endl;
    if(tauA==-1){
        std::cout << "Something is wrong in tauA file for " << mass << " GeV mass version 8, ... Exiting!" << std::endl;
        return;
    }
    
    
    //// store the xsec in a map for all the photon energy of the same mALP
    map<double, double> xSecValueMap;
    
    for(int j=0; j < 175; ++j){
        float photonBeamEnergy = 0.1*j;
        if(photonBeamEnergy < mALP)continue;
        double xsec = getXsecFromFile(mass, photonBeamEnergy, "8")*getLambdaFactor(oneOverLambda);
        xSecValueMap.insert(std::pair<double, double>(photonBeamEnergy, xsec));
        if(debug)std::cout << "photon energy: " << photonBeamEnergy << " xsec: " << xsec << std::endl;
    }
    
    
    
    Long64_t nbytes(0), nb(0);
    int    failEnergyCut(0);
    double weightAllPassEachBeam[175]      = {0.0};
    double weightNoEnergyPassEachBeam[175] = {0.0};
    double weightTotalALPEachBeam[175]     = {0.0};
    double inPhotonBeam[175]               = {0.0};
    double xsecFromMG[175]                 = {0.0};
    double totalExpFactor[175]             = {0.0};
    double expFactorFullRange[175]         = {0.0};
    double totalWeightAllPass[175]         = {0.0};
    double totalWeightNoEnergyPass[175]    = {0.0};
    double totalALP[175]                   = {0.0};
    double nTotalEvents                    = 10000.0;
    double incomingBeamNormWeight          = getWeightFactor(mALP);
    std::cout << "The incoming beam normalization weight: " << incomingBeamNormWeight << std::endl;
    
    
    double rhoT               = 11.35;
    double X0                 = 0.5612;
    double AT                 = 207;
    double M0                 = 1.661*pow(10,-24);
    double weightHistogram    = 1.0;
    double expFactor          = 1.0;
    double dwdeForEachBeam    = 1.0;
    
    
    
    ///loop over all the entries
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      //if(jentry>0)break;
      if(jentry%10000==0)std::cout << "processed: " << jentry << std::endl;
      // loop over particles
      
      double incomingPhotonBeam = -999.0;
      ///select only the incoming photon beam
      for(int nParticle=0; nParticle < Particle_;++nParticle){
        if(Particle_Status[nParticle] == -1 && Particle_PID[nParticle] == 22)
            incomingPhotonBeam = Particle_E[nParticle];
      }
      
      /// need one specific photon beam
      /// if(std::abs(incomingPhotonBeam-2.2)>0.01)continue;
      if(debug){
        std::cout << "----------------------" << std::endl;
        std::cout << "This is event number: " << jentry << std::endl;
        std::cout << "Incoming photon beam  " << incomingPhotonBeam << " GeV" << std::endl;
      }
      
      int incomingPhotonOrder(-999);
      vector<int> outgoingPhotonOrder;
      
      for (int nParticle = 0; nParticle < Particle_; ++nParticle){
          /// select the incoming photon
          if(Particle_Status[nParticle] == -1 && Particle_PID[nParticle] == 22){
              /// selecting the weight for the incoming photon beam
              map<double, vector<double> >::iterator itr; 
              incomingPhotonOrder = nParticle;
              for (itr = weightInput.begin(); itr != weightInput.end(); ++itr) {
                 if(std::abs(itr->first - Particle_E[nParticle]) < 0.04 ){
                    //cout << itr->first << " : " << Particle_E[nParticle] << endl;
                     /// getting the weights, 30k is the total number of incoming photons
                    weightHistogram = itr->second.at(0)/nTotalEvents;
                    dwdeForEachBeam = itr->second.at(1);
                 }
              }
              
              /// fill the incoming photon distributions
              h_IncomingPhotonE->Fill(Particle_E[nParticle], weightHistogram); 
              h_IncomingPhotonPhi->Fill(Particle_Phi[nParticle], weightHistogram);
              h_IncomingPhotonEta->Fill(Particle_Eta[nParticle], weightHistogram);
              
          }
          /// select the outgoing photon
          /// photon without any cuts
          if(Particle_Status[nParticle] == 1 && Particle_PID[nParticle] == 22){
              h_OutgoingPhotonE_NoCut->Fill(Particle_E[nParticle], weightHistogram);
          }
          
          if(Particle_Status[nParticle] == 1 && Particle_PID[nParticle] == 22 && Particle_E[nParticle] > photonThreshold){
              h_OutgoingPhotonE->Fill(Particle_E[nParticle], weightHistogram);
              outgoingPhotonOrder.push_back(nParticle);
          }
       } // particle loop ends
       
       /// select high end beam
       /// if(Particle_E[incomingPhotonOrder] < 10.0)continue;
       /// only select events with two photons having the threshold cut, otherwise no need to go further
       if(outgoingPhotonOrder.size() < 2)continue;
       
       ///  select two outgoing photons
       vector<TLorentzVector> photon;
       
       for(size_t i = 0; i < outgoingPhotonOrder.size(); ++i){
           TLorentzVector photonSample;
           photonSample.SetPtEtaPhiE(Particle_PT[outgoingPhotonOrder.at(i)], Particle_Eta[outgoingPhotonOrder.at(i)], Particle_Phi[outgoingPhotonOrder.at(i)], Particle_E[outgoingPhotonOrder.at(i)]);
           photon.push_back(photonSample);
       }
       
       /// opening angle
       double openingangle = std::abs(openingAngle(photon.at(0), photon.at(1)));
       h_OutgoingPhotonOpeningAngle_NoCut->Fill(openingangle, weightHistogram);
       
       if(photon.at(0).E() > photonThreshold && photon.at(1).E() > photonThreshold)h_OutgoingPhotonOpeningAngle->Fill(openingangle, weightHistogram);
       
       /// get the properties of ALP by adding the two photons it decay into.
       TLorentzVector alp = photon.at(0) + photon.at(1);
       
       double alpEta      = alp.Eta();
       double alpTheta    = thetaCal(alpEta); // theta in degrees
       
       //// from here make the acceptance plots
       alpTheta = (TMath::Pi()/180.0)*alpTheta; // converting theta to radian
    
       
       //// add the parameters needed for ALP travel
       
       double pA         = alp.P(); // taking momentum from MG //sqrt(pow(Particle_E[incomingPhotonOrder],2) - pow(mALP,2));
       double LA         = ctauA*pA/mALP;
       double LS         = 0.5;
       double LD         = 3.0;
       expFactor  = (exp(-LS/LA) - exp(-(LD+LS)/LA));
       
       if(debug){
            std::cout << "The width:        " << width << std::endl; 
            std::cout << "The tauA:         " << tauA << std::endl;
            std::cout << "The tauATheory:   " << tauATheory << std::endl;
            std::cout << "The ALP momentum: " << pA << std::endl;
            std::cout << "The pA/mALP is:   " << pA/mALP << std::endl;
            std::cout << "ctauA:            " << ctauA << std::endl; 
            std::cout << "The LA value for: " << mALP << " GeV is " << LA << std::endl;
            std::cout << "The expFactor:    " << expFactor << std::endl;
       }
       
       // get random number from a function
       
       TRandom3 *rnd = new TRandom3();
       rnd->SetSeed();

       // draw the random number r 100 times
       double photonAllPass(0.0), photonNoEnergyPass(0.0), photonAll(0.0);
       double theta[2] = {-999.0,-999.0};
       theta[0] = TMath::Pi()/180.0*thetaCal(Particle_Eta[outgoingPhotonOrder.at(0)]);
       theta[1] = TMath::Pi()/180.0*thetaCal(Particle_Eta[outgoingPhotonOrder.at(1)]);
       
       double deltaTheta = dThetaCalc(theta[0], theta[1]);
       if(Particle_E[outgoingPhotonOrder.at(0)] > 0.5 && Particle_E[outgoingPhotonOrder.at(1)] > 0.5)
            h_OutgoingPhotonOpeningTheta->Fill(deltaTheta, weightHistogram);
       
       
       double totalR  = 0.0;
       int zAxisFail  = 0;
       int radiusFail = 0;
       int bothFail   = 0;
       /// loop to randomize the decay point of ALP
       int m = 0;
       for(m=0;m<1000;++m){
           double prob = rnd->Exp(LA); // means you get exp( -t/tau )
           double r = prob;
           h_rALP->Fill(r);
           h_probALP->Fill(prob);
           totalR += r;
           
           
           
           
           /// block to get occupancy plots after each meter.
           /// get x,y,z of the ALP
           double initialZ = 0.5; /// the position where the dump starts
           r = 1.0/cos(alpTheta); /// fixed the position of decay of the alp, end of the dump in z
           double xALP     = r*sin(alpTheta)*cos(alp.Phi());
           double yALP     = r*sin(alpTheta)*sin(alp.Phi());
           double zALP     = r*cos(alpTheta);
           
           bool zAxisPass = (zALP >= LS && zALP <= (LS+LD));
           
           /// working with tungsten dump
           double LTungstenS               = 1.5;
           double LTungstenD               = 2.5;
           double alpVelocity              = abs(alp.P()/alp.E()); // this is beta
           float timeALPDecay              = initialZ*3.33+r*3.33/alpVelocity; // in ns, light travels 1 m in 3.33 ns, so alp will travel beta m in 3.33 ns, 0.5 m is the position of the begining of the dump
           float rStillToReachDetector1    = ((LTungstenS+LTungstenD)-(zALP+initialZ))/cos(theta[0]);
           float timeStillToReachDetector1 = rStillToReachDetector1*3.33; // in ns, this is photon
           float rStillToReachDetector2    = ((LTungstenS+LTungstenD)-(zALP+initialZ))/cos(theta[1]);
           float timeStillToReachDetector2 = rStillToReachDetector2*3.33; // in ns, this is photon
           
           
           /// fill the timing histogram
           if(Particle_E[outgoingPhotonOrder.at(0)] > 0.5 && Particle_E[outgoingPhotonOrder.at(1)] > 0.5)h_ALPDecayTime->Fill(timeALPDecay, weightHistogram/1000.0);
           if(Particle_E[outgoingPhotonOrder.at(0)] > 0.1 && Particle_E[outgoingPhotonOrder.at(1)] > 0.1)h_ALPDecayTimeEGt0p1->Fill(timeALPDecay, weightHistogram/1000.0);
           
           if(timeStillToReachDetector1 < 0 || timeStillToReachDetector2 < 0)continue;
           
           double x1mTungsten[2], y1mTungsten[2];
           x1mTungsten[0] = xALP + ((LTungstenS+LTungstenD)-(zALP+initialZ))*tan(theta[0])*cos(Particle_Phi[outgoingPhotonOrder.at(0)]); 
           x1mTungsten[1] = xALP + ((LTungstenS+LTungstenD)-(zALP+initialZ))*tan(theta[1])*cos(Particle_Phi[outgoingPhotonOrder.at(1)]);
           y1mTungsten[0] = yALP + ((LTungstenS+LTungstenD)-(zALP+initialZ))*tan(theta[0])*sin(Particle_Phi[outgoingPhotonOrder.at(0)]);
           y1mTungsten[1] = yALP + ((LTungstenS+LTungstenD)-(zALP+initialZ))*tan(theta[1])*sin(Particle_Phi[outgoingPhotonOrder.at(1)]);
           
           
           //// for alps decaying in the correct window
           if((zALP+initialZ) >= LTungstenS && (zALP+initialZ) <= (LTungstenS+LTungstenD)){
              /// leading photon 
              /// photons need to be captured by the detector
              if( (sqrt(pow(x1mTungsten[0],2)+pow(y1mTungsten[0],2)) <= 1.0) && (sqrt(pow(x1mTungsten[1],2)+pow(y1mTungsten[1],2)) <= 1.0) ){
                if(Particle_E[outgoingPhotonOrder.at(0)] > 0.5 && Particle_E[outgoingPhotonOrder.at(1)] > 0.5){
                    h_ALPArrivalTime->Fill(timeALPDecay+timeStillToReachDetector1, weightHistogram/1000.0);
                    h_ALPArrivalTime->Fill(timeALPDecay+timeStillToReachDetector2, weightHistogram/1000.0);
                }
                if(Particle_E[outgoingPhotonOrder.at(0)] > 0.1 && Particle_E[outgoingPhotonOrder.at(1)] > 0.1){
                    h_ALPArrivalTimeEGt0p1->Fill(timeALPDecay+timeStillToReachDetector1, weightHistogram/1000.0);
                    h_ALPArrivalTimeEGt0p1->Fill(timeALPDecay+timeStillToReachDetector2, weightHistogram/1000.0);
                }
              }
           }
           
           
           
           
           //// back to case where we have dump at origin 0
           double x1m[2], y1m[2];
           x1m[0] = xALP + (LD+LS - zALP)*tan(theta[0])*cos(Particle_Phi[outgoingPhotonOrder.at(0)]); 
           x1m[1] = xALP + (LD+LS - zALP)*tan(theta[1])*cos(Particle_Phi[outgoingPhotonOrder.at(1)]);
           y1m[0] = yALP + (LD+LS - zALP)*tan(theta[0])*sin(Particle_Phi[outgoingPhotonOrder.at(0)]);
           y1m[1] = yALP + (LD+LS - zALP)*tan(theta[1])*sin(Particle_Phi[outgoingPhotonOrder.at(1)]);
           
           if(zALP >= LS && zALP <= LS+0.05){
               h_PhotonY_PhotonX->Fill(x1m[0], y1m[0], weightHistogram/1000.0);
               h_PhotonY_PhotonX->Fill(x1m[1], y1m[1], weightHistogram/1000.0);
               h_ALPY_ALPX->Fill(xALP, yALP, weightHistogram/1000.0);
               
               if(Particle_E[outgoingPhotonOrder.at(0)] > 0.5 && Particle_E[outgoingPhotonOrder.at(1)] > 0.5){
                   h_PhotonY_PhotonX_E0p5GeV->Fill(x1m[0], y1m[0], weightHistogram/1000.0);
                   h_PhotonY_PhotonX_E0p5GeV->Fill(x1m[1], y1m[1], weightHistogram/1000.0);
                   h_ALPY_ALPX_E0p5GeV->Fill(xALP, yALP, weightHistogram/1000.0);
               }
           }
           
           if(zAxisPass){
              if(Particle_E[outgoingPhotonOrder.at(0)] > 0.5 && Particle_E[outgoingPhotonOrder.at(1)] > 0.5)h_ALPRxy_Z_E0p5GeV->Fill(zALP, r*sin(alpTheta), weightHistogram/1000.0);
           }
           
           float radius1          = sqrt(pow(x1m[0],2)+pow(y1m[0],2));
           float radius2          = sqrt(pow(x1m[1],2)+pow(y1m[1],2));
           
           //// fill this histogram only once per ALP
           if(m==0){
                if(Particle_E[outgoingPhotonOrder.at(0)] > 0.5 && Particle_E[outgoingPhotonOrder.at(1)] > 0.5){
                    h_DistanceInDetector->Fill(radius1, weightHistogram);
                    h_DistanceInDetector->Fill(radius2, weightHistogram); 
                    h_DistanceInDetectorLessBins->Fill(radius1, weightHistogram);
                    h_DistanceInDetectorLessBins->Fill(radius2, weightHistogram);
                }
           }
           
           bool photon1RadiusPass = (sqrt(pow(x1m[0],2)+pow(y1m[0],2)) <= 1.0);
           bool photon2RadiusPass = (sqrt(pow(x1m[1],2)+pow(y1m[1],2)) <= 1.0);
           bool photon1EnergyPass = (Particle_E[outgoingPhotonOrder.at(0)] > 0.5);
           bool photon2EnergyPass = (Particle_E[outgoingPhotonOrder.at(1)] > 0.5);
           
           //if(debug)std::cout << "radius 1 " << radius1 << " radius2 " << radius2 << " energy1: " << Particle_E[outgoingPhotonOrder.at(0)] << " energy2: " << Particle_E[outgoingPhotonOrder.at(1)] <<  std::endl;
           
           if(!zAxisPass)zAxisFail++;
           if(!(photon1RadiusPass && photon2RadiusPass))radiusFail++;
           if(!(zAxisPass && photon1RadiusPass && photon2RadiusPass))bothFail++;
           
           if(zAxisPass && photon1RadiusPass && photon2RadiusPass && photon1EnergyPass && photon2EnergyPass)photonAllPass++; /// removing the energy cut for now ---
           if(zAxisPass && photon1RadiusPass && photon2RadiusPass && (photon1EnergyPass || photon2EnergyPass))photonNoEnergyPass++;
           
           photonAll++;
       }
       
       double weightAllPass      = photonAllPass/photonAll;
       double weightNoEnergyPass = photonNoEnergyPass/photonAll;
       double aveR               = totalR/photonAll;
       /// save the acceptance for each incoming photon beam
       int beamEntry             = jentry/int(nTotalEvents);
       
       totalWeightAllPass[beamEntry]       += weightAllPass*weightHistogram*nTotalEvents;    /// adding the weight due to incoming photon spectrum
       totalWeightNoEnergyPass[beamEntry]  += weightNoEnergyPass*weightHistogram*nTotalEvents; /// adding the weight due to incoming photon spectrum
       totalALP[beamEntry]                 += weightHistogram*nTotalEvents;
       
       /// exponential factor
       totalExpFactor[beamEntry]           += expFactor*weightHistogram*nTotalEvents;
       
       /// storing the acceptance for each beam chunk
       if(totalALP[beamEntry]!=0){
           weightAllPassEachBeam[beamEntry]       = totalWeightAllPass[beamEntry]/totalALP[beamEntry]; // *incomingBeamNormWeight
           weightNoEnergyPassEachBeam[beamEntry]  = totalWeightNoEnergyPass[beamEntry]/totalALP[beamEntry]; // *incomingBeamNormWeight
           expFactorFullRange[beamEntry]          = totalExpFactor[beamEntry]/totalALP[beamEntry];
       }
       else{
           weightAllPassEachBeam[beamEntry]       = 0;
           weightNoEnergyPassEachBeam[beamEntry]  = 0;
       }
       
       /// just the number of photons expected per beam chunk
       weightTotalALPEachBeam[beamEntry]    = weightHistogram*nTotalEvents;
       inPhotonBeam[beamEntry]              = Particle_E[incomingPhotonOrder];
       /// xsec for each beam chunk, take it only for the first time so that the reading is not done for all values
       if(xsecFromMG[beamEntry] == 0){
           map<double, double>::iterator xsecitr;
           for (xsecitr = xSecValueMap.begin(); xsecitr != xSecValueMap.end(); ++xsecitr) {
               if(std::abs(xsecitr->first - Particle_E[incomingPhotonOrder]) < 0.04 ){
                  xsecFromMG[beamEntry]  = xsecitr->second;
               }
           }
       }
       
       if(debug){
            std::cout << "Mass: " << mALP << " photon beam: " << beamEntry << " xsecFromMap: " << xsecFromMG[beamEntry] << " xSecFromFile: " << getXsecFromFile(mass, Particle_E[incomingPhotonOrder], "8")*getLambdaFactor(oneOverLambda) << std::endl;
       }
       
       if(xsecFromMG[beamEntry] < 0){
           std::cout << "Something is wrong in xsec file for " << mass << " GeV mass and beam energy: " << Particle_E[incomingPhotonOrder] << " .... Exiting!" << std::endl;
           return;
       }
       
       
       if(debug && (jentry+1)%int(nTotalEvents)==0){
           std::cout << "beam entry :                        " << beamEntry << std::endl;
           std::cout << "From array weightAllPassEachBeam:   " << weightAllPassEachBeam[beamEntry] << std::endl;
           std::cout << "From double totalWeightAllPass:     " << totalWeightAllPass[beamEntry] << std::endl;
           std::cout << "From double totalALP:               " << totalALP[beamEntry] << std::endl;
           std::cout << "From xsec:                          " << xsecFromMG[beamEntry] << std::endl;
           std::cout << "getLambdaFactor:                    " << getLambdaFactor(oneOverLambda) << std::endl;
       }
       
       
       
       if((Particle_E[outgoingPhotonOrder.at(0)] < 0.5 || Particle_E[outgoingPhotonOrder.at(1)] < 0.5))failEnergyCut++;
       
       if(debug){
            std::cout << "From the " << m << " trial, failed trial for zAxis cut : " << zAxisFail << std::endl;
            std::cout << "From the " << m << " trial, failed trial for radius cut: " << radiusFail << std::endl;
            std::cout << "From the " << m << " trial, failed trial for both cut  : " << bothFail << std::endl;
            std::cout << "From the " << m << " trial, all pass for both cut      : " << photonAllPass << std::endl;
            std::cout << "From the " << m << " trial, no energy pass             : " << photonNoEnergyPass << std::endl;
            std::cout << "From the " << m << " trial, weightAllPass              : " << weightAllPass << std::endl;
            std::cout << "From the " << m << " trial, weightNoEnergyPass         : " << weightNoEnergyPass << std::endl;
            std::cout << "Incoming beam norm weight                        : " << incomingBeamNormWeight << std::endl;
            std::cout << "The average r value                              : " << aveR << std::endl;
            std::cout << "From xsec                                        : " << xsecFromMG[beamEntry] << std::endl;
       }
       
       h_RatioLA->Fill(LA/aveR);
       
    }/// loop over all events end
    
    /// multiply the prefactor
    double preFactor      = (rhoT*X0)/(AT*M0); 
    /// loop over the acceptance arrays for each incoming photon beam
    double sumXsecTimesWeight(0.0), sumXsecTimesWeightNoEnergyPass(0.0);
    for(int k=0; k < 175; ++k){
        //// acceptance*xsec*number of photons from electron*expFactor
        if(debug){
            std::cout << "Xsec from the loop           : " << xsecFromMG[k] << std::endl;
            std::cout << "Acceptance from the loop     : " << weightAllPassEachBeam[k] << std::endl;
            std::cout << "N_i for this beam            : " << weightTotalALPEachBeam[k] << std::endl;
            std::cout << "Exp factor from the loop     : " << expFactorFullRange[k] << std::endl; 
        }
        sumXsecTimesWeight             += weightAllPassEachBeam[k]*xsecFromMG[k]*weightTotalALPEachBeam[k]; // *expFactorFullRange[k]*incomingBeamNormWeight
        sumXsecTimesWeightNoEnergyPass += weightNoEnergyPassEachBeam[k]*xsecFromMG[k]*weightTotalALPEachBeam[k]; // *expFactorFullRange[k]*incomingBeamNormWeight
    }
    
    double expectedNumber = (preFactor)*pow(10,-36)*sumXsecTimesWeight; /// converting the /cm2 to picobarn
    if(debug){
        std::cout << "... After the loop ..." << std::endl;
        std::cout << "The preFactor:                                  " << preFactor << std::endl;
        std::cout << "The prefactor times Ne:                         " << preFactor*1.5*pow(10,16) << std::endl;
        //std::cout << "Theory - xsec*luminosity*acceptance*expFactor : " << preFactor*pow(10,-36)*1.5*pow(10,16)*getXsecFromFile(mass, 5.0, "8")*getLambdaFactor(oneOverLambda)*expFactorFullRange[0]*weightAllPassEachBeam[0]*dwdeForEachBeam << std::endl;
        std::cout << "The average xsec:                               " << xsecFromMG[0] << " pb" << std::endl;
        std::cout << "The average acceptance:                         " << weightAllPassEachBeam[0] << std::endl;
        std::cout << "The weight value:                               " << weightTotalALPEachBeam[0] << std::endl;
        std::cout << "The average expFactor:                          " << expFactorFullRange[0] << std::endl;
        std::cout << "Code - xsec*luminosity*acceptance*expFactor:    " << expectedNumber << std::endl;
    }
    acceptanceFiles << mass << "\t" << oneOverLambda << "\t" << expectedNumber  << std::endl;
    acceptanceFiles.close();
    fOut->Write();
    fOut->Close();
    // Record end time
    auto finish = std::chrono::steady_clock::now();
    auto diff = finish - start;
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time : " << chrono::duration <double, milli> (diff).count()/1000.0 << " s" << endl;
   
}
