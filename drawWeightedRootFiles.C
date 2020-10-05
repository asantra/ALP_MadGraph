#define drawWeightedRootFiles_cxx
#include "drawWeightedRootFiles.h"
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
// 
using namespace std;

double openingAngle(TLorentzVector a, TLorentzVector b){
    double angle = a.Angle(b.Vect());
    return angle;
}

float dPhiCalc(float phiLead, float phiTrail){
  float dphi = fabs(phiLead - phiTrail);
  if(dphi > TMath::Pi()) dphi = TMath::Pi()*2. - dphi;
  return dphi;
}

float thetaCal(float eta){
    double theta = 2*atan(exp(-eta));
    return theta*180/TMath::Pi();
}

double getTau(float mALP){
      double tauA = 1.0; 
      if(mALP == 0.001)
           tauA = 1.0/(4.97359*pow(10,-18));
       else if(mALP == 0.01)
           tauA = 1.0/(4.97359*pow(10,-15));
       else if(mALP == 0.05)
           tauA = 1.0/(6.21699*pow(10,-13));
       else if(mALP == 0.5)
           tauA = 1.0/(6.21699*pow(10,-16)); // this is for lambda 1e6, not lambda 1e3
       else
           tauA = 1.0;
       
       return tauA;
}

void drawWeightedRootFiles::Loop(TString mass, TString version, TString cutPhoton)
{
//   In a ROOT session, you can do:
//      root> .L drawWeightedRootFiles.C
//      root> drawWeightedRootFiles t
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

    TFile *fOut                         = new TFile("PhotonAcceptanceFiles_ALPMass"+mass+"GeV"+photonNameFile+"_V"+version+".root", "RECREATE");
    fOut->cd();
    
    ///1D histogram
    TH1F *h_IncomingPhotonE             = new TH1F("h_IncomingPhotonE", "Incoming photon energy; E [GeV]; Events", 210, 0, 10.5);
    TH1F *h_IncomingPhotonPhi           = new TH1F("h_IncomingPhotonPhi", "Incoming photon #phi; #phi; Events", 64, -3.2, 3.2);
    TH1F *h_IncomingPhotonEta           = new TH1F("h_IncomingPhotonEta", "Incoming photon #eta; #eta; Events", 3000, -1500.0, 1500.0);
    TH1F *h_OutgoingPhotonE_NoCut       = new TH1F("h_OutgoingPhotonE_NoCut", "Outgoing photon energy, without cuts; E [GeV]; Events", 500, 0, 20);
    TH1F *h_OutgoingPhotonE             = new TH1F("h_OutgoingPhotonE", "Outgoing photon energy; E [GeV]; Events", 500, 0, 20);
    
    
    
    TH1F *h_rALP                        = new TH1F("h_rALP", "h_rALP",100,0,10);
    
    
    ///2D histograms
    
    
    /// all around the plane
    TH2F *h_TrackingPlane1m_XY = new TH2F("h_TrackingPlane1m_XY", "occupancy plot, 1 m from IP; X [m]; Y [m]", 1000,-50,50,1000,-50,50);
    
    Long64_t nentries = fChain->GetEntriesFast();
    double eachEntry  = 30000.0;
    map<double, vector<double> > weightInput;
    
    
    /// getting the xsec and nPhoton weights from a text file.
    string line;
    ifstream myfile("outputPhotonNumbers.txt");
    string dummyLine;
    std::getline(myfile, dummyLine); // go over the first line which is not numbers
    double energy, xsecv3, xsecv4, xsecv5, xsecv6, nPhoton;
    
    /// loop over the numbered lines in the text file
    if (myfile.is_open()){
        while (myfile >> energy >> xsecv3 >> xsecv4 >> xsecv5 >> xsecv6 >> nPhoton){
            vector<double> photonDetails;
            photonDetails.push_back(nPhoton);
            photonDetails.push_back(xsecv3);
            weightInput.insert(std::pair<double,vector<double> >(energy, photonDetails));
            //std::cout << energy << " " << xsecv3 << " " << xsecv4 << " " << xsecv5 << " " << xsecv6 << " " << nPhoton << std::endl;
        }
        myfile.close();
    }
    
    
    Long64_t nbytes = 0, nb = 0;
    double totalWeightAllPass(0.0), totalWeightEnergyPass(0.0), totalALP(0.0);
    int failEnergyCut(0);
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      //if(jentry>0)break;
      if(jentry%10000==0)std::cout << "processed: " << jentry << std::endl;
      // loop over particles
      double weightHistogram = 1.0;
      double incomingPhotonBeam = -999.0;
      ///select only 2.2 GeV beam
      for(int nParticle=0; nParticle < Particle_;++nParticle){
        if(Particle_Status[nParticle] == -1 && Particle_PID[nParticle] == 22)
            incomingPhotonBeam = Particle_E[nParticle];
      }
      
      /// need one specific photon beam
      if(std::abs(incomingPhotonBeam-2.2)>0.01)continue;
      std::cout << "Incoming photon beam " << incomingPhotonBeam << std::endl;
      
      int incomingPhotonOrder(-999);
      vector<int> outgoingPhotonOrder;
      
      for (int nParticle = 0; nParticle < Particle_; ++nParticle){
          ///select the incoming photon
          if(Particle_Status[nParticle] == -1 && Particle_PID[nParticle] == 22){
              //std::cout << "Incoming photon energy " << Particle_E[nParticle] << std::endl;
              
              /// selecting the weight for the incoming photon beam
              map<double, vector<double> >::iterator itr; 
              incomingPhotonOrder = nParticle;
//               for (itr = weightInput.begin(); itr != weightInput.end(); ++itr) {
//                  if(std::abs(itr->first - Particle_E[nParticle]) < 0.04 ){
//                     //cout << itr->first << " : " << Particle_E[nParticle] << endl;
//                      /// getting the weights, 30k is the total number of incoming photons
//                     weightHistogram = itr->second.at(0)/30000.0;
//                  }
//               }
              
              /// fill the incoming photon distributions
              h_IncomingPhotonE->Fill(Particle_E[nParticle], weightHistogram); 
              h_IncomingPhotonPhi->Fill(Particle_Phi[nParticle], weightHistogram);
              h_IncomingPhotonEta->Fill(Particle_Eta[nParticle], weightHistogram);
              
          }
          /// select the outgoing photon
          /// photon without any cuts
          if(Particle_Status[nParticle] == 1 && Particle_PID[nParticle] == 22)h_OutgoingPhotonE_NoCut->Fill(Particle_E[nParticle], weightHistogram);
          
          if(Particle_Status[nParticle] == 1 && Particle_PID[nParticle] == 22 && Particle_E[nParticle] > photonThreshold){
              h_OutgoingPhotonE->Fill(Particle_E[nParticle], weightHistogram);
              outgoingPhotonOrder.push_back(nParticle);
          }
       }
       /// only select events with two photons having the threshold cut, otherwise no need to go further
       if(outgoingPhotonOrder.size()<2)continue;
       
       ///  select two outgoing photons
       vector<TLorentzVector> photon;
       
       for(size_t i = 0; i < outgoingPhotonOrder.size(); ++i){
           TLorentzVector photonSample;
           photonSample.SetPtEtaPhiE(Particle_PT[outgoingPhotonOrder.at(i)], Particle_Eta[outgoingPhotonOrder.at(i)], Particle_Phi[outgoingPhotonOrder.at(i)], Particle_E[outgoingPhotonOrder.at(i)]);
           photon.push_back(photonSample);
       }
       /// get the properties of ALP by adding the two photons it decay into.
       TLorentzVector alp = photon.at(0) + photon.at(1);
       
       float alpEta   = alp.Eta();
       float alpTheta = thetaCal(alpEta);
       
       //// from here make the acceptance plots
       alpTheta = (TMath::Pi()/180.0)*alpTheta; // converting theta to radian
       
       
       
       //// add the parameters needed for ALP travel
       float mALP       = mass.Atof(); /// taking different ALP mass values
       float tauA       = getTau(mALP);
       float ctauA      = tauA*(0.1975)*pow(10,-15); // c is 1 in natural units, the additional factor is the conversion factor from GeV^{-1} to meters
       float pA         = alp.P(); // taking momentum from MG //sqrt(pow(Particle_E[incomingPhotonOrder],2) - pow(mALP,2));
       float LA         = ctauA*pA/mALP;
       std::cout << "The tauA: " << tauA << std::endl;
       std::cout << "The ALP momentum: " << pA << std::endl;
       std::cout << "The pA/mALP is " << pA/mALP << std::endl;
       std::cout << "ctauA : " << ctauA << std::endl; 
       std::cout << "The LA value for " << mALP << " GeV is " << LA << std::endl;
       float LS         = 0.5;
       float LD         = 5.5;
       double expFactor = exp(-(LS/LA)) - exp(-((LD+LS)/LA));
       
       TF1 *f1 = new TF1("f1","1/[0]*exp(-x/[0])",0.0,1.0);
       f1->SetParameter(0,LA);
       f1->SetParName(0,"decay length");
       double r = 0.0; // this is the r where ALP is decaying
       // draw the random number r 100 times
       double photonAllPass(0.0), photonEnergyPass(0.0), photonAll(0.0);
       float theta[2] = {-999.0,-999.0};
       theta[0] = TMath::Pi()/180.0*thetaCal(Particle_Eta[outgoingPhotonOrder.at(0)]);
       theta[1] = TMath::Pi()/180.0*thetaCal(Particle_Eta[outgoingPhotonOrder.at(1)]);
       double totalR = 0.0;
       int zAxisFail  = 0;
       int radiusFail = 0;
       for(int m=0;m<100;++m){
           double prob = f1->GetRandom();
           r = -LA*log(prob*LA);
           h_rALP->Fill(r);
           totalR += r;
           /// block to get occupancy plots after each meter.
           /// get x,y,z of the ALP
           float xALP   = r*sin(alpTheta)*cos(alp.Phi());
           float yALP   = r*sin(alpTheta)*sin(alp.Phi());
           float zALP   = r*cos(alp.Phi());
           
           bool zAxisPass   = (zALP >= LS && zALP <= (LS+LD));
           double x1m[2], y1m[2];
           x1m[0] = xALP + (LD+LS - zALP)*tan(theta[0])*cos(Particle_Phi[outgoingPhotonOrder.at(0)]); 
           x1m[1] = xALP + (LD+LS - zALP)*tan(theta[1])*cos(Particle_Phi[outgoingPhotonOrder.at(1)]);
           y1m[0] = yALP + (LD+LS - zALP)*tan(theta[0])*sin(Particle_Phi[outgoingPhotonOrder.at(0)]);
           y1m[1] = yALP + (LD+LS - zALP)*tan(theta[1])*sin(Particle_Phi[outgoingPhotonOrder.at(1)]);
           
           bool photon1RadiusPass = (sqrt(pow(x1m[0],2)+pow(y1m[0],2)) <= 1.0);
           bool photon2RadiusPass = (sqrt(pow(x1m[1],2)+pow(y1m[1],2)) <= 1.0);
           bool photon1EnergyPass = (Particle_E[outgoingPhotonOrder.at(0)] > 0.5);
           bool photon2EnergyPass = (Particle_E[outgoingPhotonOrder.at(1)] > 0.5);
           
           if(!zAxisPass)zAxisFail++;
           if(!(photon1RadiusPass && photon2RadiusPass))radiusFail++;
           
           if(zAxisPass && photon1RadiusPass && photon2RadiusPass && photon1EnergyPass && photon2EnergyPass)photonAllPass++;
           if(zAxisPass && photon1RadiusPass && photon2RadiusPass && (photon1EnergyPass || photon2EnergyPass))photonEnergyPass++;
        
           photonAll++;
           
       }
       std::cout << "From the 100 trial, failed trial for zAxis cut: " << zAxisFail << std::endl;
       std::cout << "From the 100 trial, failed trial for radius cut: " << radiusFail << std::endl;
       std::cout << "The average r value " << totalR/photonAll << std::endl;
       
       double weightAllPass    = photonAllPass/photonAll;
       double weightEnergyPass = photonEnergyPass/photonAll;
       
       totalWeightAllPass    += weightAllPass;
       totalWeightEnergyPass += weightEnergyPass;
       totalALP++;
       
       if((Particle_E[outgoingPhotonOrder.at(0)] < 0.5 || Particle_E[outgoingPhotonOrder.at(1)] < 0.5))failEnergyCut++;
       
       
    }/// loop over all events end
    
    double acceptanceAllPass    = totalWeightAllPass / totalALP;
    double acceptanceEnergyPass = totalWeightEnergyPass / totalALP;
    
    std::cout << "----------------------------" << std::endl;
    std::cout << "Total weight for all cuts: " << totalWeightAllPass << " and totalALP: " << totalALP << std::endl;
    std::cout << "The acceptance for all cuts: " << acceptanceAllPass << std::endl;
    std::cout << "Total weight for energy cuts: " << totalWeightEnergyPass << " and totalALP: " << totalALP << std::endl;
    std::cout << "The acceptance for energy cuts: " << acceptanceEnergyPass << std::endl;
    std::cout << "Failed the energy cut " << failEnergyCut << std::endl;
    
    fOut->Write();
    fOut->Close();
    // Record end time
    auto finish = std::chrono::steady_clock::now();
    auto diff = finish - start;
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time : " << chrono::duration <double, milli> (diff).count()/1000.0 << " s" << endl;
   
}
