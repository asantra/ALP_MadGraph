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
           tauA = 1.0/(6.21699*pow(10,-16));
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

    TFile *fOut                         = new TFile("combinedOutputFile_ALPMass"+mass+"GeV"+photonNameFile+"_V"+version+".root", "RECREATE");
    fOut->cd();
    
    ///1D histogram
    TH1F *h_IncomingPhotonE             = new TH1F("h_IncomingPhotonE", "Incoming photon energy; E [GeV]; Events", 210, 0, 10.5);
    TH1F *h_IncomingPhotonPhi           = new TH1F("h_IncomingPhotonPhi", "Incoming photon #phi; #phi; Events", 64, -3.2, 3.2);
    TH1F *h_IncomingPhotonEta           = new TH1F("h_IncomingPhotonEta", "Incoming photon #eta; #eta; Events", 3000, -1500.0, 1500.0);
    TH1F *h_OutgoingPhotonE_NoCut       = new TH1F("h_OutgoingPhotonE_NoCut", "Outgoing photon energy, without cuts; E [GeV]; Events", 500, 0, 20);
    TH1F *h_OutgoingPhotonE             = new TH1F("h_OutgoingPhotonE", "Outgoing photon energy; E [GeV]; Events", 500, 0, 20);
    
    
    TH1F *h_OutgoingPhotonPT            = new TH1F("h_OutgoingPhotonPT", "Outgoing photon Pt; p_{T} [GeV]; Events", 500, 0, 20);
    TH1F *h_OutgoingPhotonEta           = new TH1F("h_OutgoingPhotonEta", "Outgoing photon Eta; #eta; Events", 400, -10, 10);
    TH1F *h_OutgoingPhotonPhi           = new TH1F("h_OutgoingPhotonPhi", "Outgoing photon Phi; #phi (radian); Events", 32, -3.2, 3.2);
    
    TH1F *h_OutgoingPhotonOpeningAngle  = new TH1F("h_OutgoingPhotonOpeningAngle", "Opening angle of outgoing photons; angle (radian); Events", 128, 0.0, 3.2);
    TH1F *h_OutgoingPhotonTheta         = new TH1F("h_OutgoingPhotonTheta", "Outgoing photon theta; #theta (degrees); events", 190,0,190);
    TH1F *h_IncomingPhotonTheta         = new TH1F("h_IncomingPhotonTheta", "Incoming photon theta; #theta (degrees); events", 190,0,190);
    TH1F *h_ALPTheta                    = new TH1F("h_ALPTheta", "ALP polar angle; #theta [degrees]; events", 180,0,180);
    TH1F *h_ALPThetaZoomed              = new TH1F("h_ALPThetaZoomed", "ALP polar angle; #theta [degrees]; events", 100,0,10);
    TH1F *h_ALPEnergy                   = new TH1F("h_ALPEnergy", "ALP Energy; E (GeV); events", 100,0,20);
    TH1F *h_ALPMass                     = new TH1F("h_ALPMass", "ALP Mass; Mass (GeV); events", 100,0,1);
    TH1F *h_ALPMomentum                 = new TH1F("h_ALPMomentum", "ALP Momentum; P (GeV); events", 100,0,20);
    TH1F *h_ALPTransverseMomentum       = new TH1F("h_ALPTransverseMomentum", "ALP Transverse Momentum; p_{T} (GeV); events", 100,0,20);
    
    TH1F *h_rALP                        = new TH1F("h_rALP", "h_rALP",100,0,10);
    
    
    ///2D histograms
    TH2F *h_OutgoingPhotonE_Theta                         = new TH2F("h_OutgoingPhotonE_Theta", "Outgoing photon E vs Theta; #theta [degrees];E [GeV]", 180, 0, 180, 110, 0, 11.0);
    TH2F *h_OutgoingLeadingPhotonE_Theta                  = new TH2F("h_OutgoingLeadingPhotonE_Theta", "Outgoing leading photon #theta vs E; leading photon E [GeV]; #theta [degrees]", 110, 0, 11.0, 180, 0, 180);
    TH2F *h_OutgoingSubLeadingPhotonE_Theta               = new TH2F("h_OutgoingSubLeadingPhotonE_Theta", "Outgoing subleading photon #theta vs E; subleading photon E [GeV]; #theta [degrees]", 110, 0, 11.0, 180, 0, 180);
    TH2F *h_OutgoingPhoton_Theta1VsTheta2                 = new TH2F("h_OutgoingPhoton_Theta1VsTheta2", "Outgoing photon theta1 vs theta2; leading photon #theta [degrees]; subleading photon #theta [degrees]", 180, 0, 180, 180, 0, 180);
    TH2F *h_OutgoingPhoton_Theta1VsTheta2_WithALP         = new TH2F("h_OutgoingPhoton_Theta1VsTheta2_WithALP", "Outgoing photon theta1 vs theta2 from ALP; leading photon #theta [degrees]; subleading photon #theta [degrees]", 360, -180, 180, 360, -180, 180);
    TH2F *h_OutgoingPhoton_Theta1VsTheta2_WithCut         = new TH2F("h_OutgoingPhoton_Theta1VsTheta2_WithCut", "Outgoing photon theta1 vs theta2; leading photon #theta [degrees]; subleading photon #theta [degrees]", 500, 0, 5, 500, 0, 5);
    TH2F *h_OutgoingPhoton_diffAngleVsLeadingPhotonAngle  = new TH2F("h_OutgoingPhoton_diffAngleVsLeadingPhotonAngle", "Leading photon angle vs #Delta#theta; leading photon #theta [degrees]; #Delta#theta [degrees]", 180, 0, 180, 180, 0, 180);
    
    TH2F *h_OutgoingPhoton_Energy1VsEnergy2               = new TH2F("h_OutgoingPhoton_Energy1VsEnergy2", "Outgoing photon Energy1 vs Energy2; leading photon energy [GeV]; subleading photon energy [GeV]", 105, 0, 10.5, 105, 0, 10.5);
    TH2F *h_OutgoingPhoton_diffEnergyVsLeadingPhotonEnergy= new TH2F("h_OutgoingPhoton_diffEnergyVsLeadingPhotonEnergy", "Leading photon energy vs #DeltaE; leading photon energy [GeV]; #Delta E [GeV]", 105, 0, 10.5, 105, 0, 10.5);
    
    TH2F *h_OutgoingPhoton_diffALPEnergyVsLeadingPhotonEnergy  = new TH2F("h_OutgoingPhoton_diffALPEnergyVsLeadingPhotonEnergy", "Leading photon energy vs #DeltaE with ALP; leading photon energy [GeV]; #Delta E with ALP [GeV]", 105, 0, 10.5, 105, 0, 10.5);
    
    ///leading photon energy
    TH2F *h_OutgoingPhoton_LeadingPhotonThetaVsLeadingPhotonEnergyFractionALP  = new TH2F("h_OutgoingPhoton_LeadingPhotonThetaVsLeadingPhotonEnergyFractionALP", "Leading photon #theta vs leading photon energy (fraction of ALP); fraction of leading photon energy; leading photon #theta [degrees]", 60, 0.4, 1, 180, 0, 180);
    TH2F *h_OutgoingPhoton_LeadingPhotonThetaWrtALPVsLeadingPhotonEnergyFractionALP  = new TH2F("h_OutgoingPhoton_LeadingPhotonThetaWrtALPVsLeadingPhotonEnergyFractionALP", "Leading photon #theta wrt ALP vs leading photon energy (fraction of ALP); fraction of leading photon energy; leading photon #theta wrt ALP [degrees]", 60, 0.4, 1, 180, 0, 180);
    
    
    ///subleading photon energy
    TH2F *h_OutgoingPhoton_SubLeadingPhotonThetaVsSubLeadingPhotonEnergyFractionALP  = new TH2F("h_OutgoingPhoton_SubLeadingPhotonThetaVsSubLeadingPhotonEnergyFractionALP", "SubLeading photon #theta vs Subleading photon energy (fraction of ALP); fraction of subleading photon energy; subleading photon #theta [degrees]", 60, 0, 0.6, 180, 0, 180);
    TH2F *h_OutgoingPhoton_SubLeadingPhotonThetaWrtALPVsSubLeadingPhotonEnergyFractionALP  = new TH2F("h_OutgoingPhoton_SubLeadingPhotonThetaWrtALPVsSubLeadingPhotonEnergyFractionALP", "SubLeading photon #theta wrt ALP vs subleading photon energy (fraction of ALP); fraction of subleading photon energy; subleading photon #theta wrt ALP [degrees]", 60, 0, 0.6, 180, 0, 180);
    
    /// ALP parameters
    TH2F *h_ALPThetaVsEnergy   = new TH2F("h_ALPThetaVsEnergy", "ALP #theta vs energy; E (GeV); #theta [degrees]", 100,0,20, 180,0,180);
    
    /// all around the plane
    TH2F *h_TrackingPlane1m_XY = new TH2F("h_TrackingPlane1m_XY", "occupancy plot, 1 m from IP; X [m]; Y [m]", 1000,-50,50,1000,-50,50);
    TH2F *h_TrackingPlane2m_XY = new TH2F("h_TrackingPlane2m_XY", "occupancy plot, 2 m from IP; X [m]; Y [m]", 1000,-50,50,1000,-50,50);
    TH2F *h_TrackingPlane3m_XY = new TH2F("h_TrackingPlane3m_XY", "occupancy plot, 3 m from IP; X [m]; Y [m]", 1000,-50,50,1000,-50,50);
    TH2F *h_TrackingPlane4m_XY = new TH2F("h_TrackingPlane4m_XY", "occupancy plot, 4 m from IP; X [m]; Y [m]", 1000,-50,50,1000,-50,50);
    TH2F *h_TrackingPlane5m_XY = new TH2F("h_TrackingPlane5m_XY", "occupancy plot, 5 m from IP; X [m]; Y [m]", 1000,-50,50,1000,-50,50);
    TH2F *h_TrackingPlane6m_XY = new TH2F("h_TrackingPlane6m_XY", "occupancy plot, 6 m from IP; X [m]; Y [m]", 1000,-50,50,1000,-50,50);
    
    /// only within 1m radius
    TH2F *h_TrackingPlane1m_XY_1mRadius = new TH2F("h_TrackingPlane1m_XY_1mRadius", "occupancy plot, 1 m from IP, within 1m radius; X [m]; Y [m]", 300,-1.5,1.5,300,-1.5,1.5);
    TH2F *h_TrackingPlane2m_XY_1mRadius = new TH2F("h_TrackingPlane2m_XY_1mRadius", "occupancy plot, 2 m from IP, within 1m radius; X [m]; Y [m]", 300,-1.5,1.5,300,-1.5,1.5);
    TH2F *h_TrackingPlane3m_XY_1mRadius = new TH2F("h_TrackingPlane3m_XY_1mRadius", "occupancy plot, 3 m from IP, within 1m radius; X [m]; Y [m]", 300,-1.5,1.5,300,-1.5,1.5);
    TH2F *h_TrackingPlane4m_XY_1mRadius = new TH2F("h_TrackingPlane4m_XY_1mRadius", "occupancy plot, 4 m from IP, within 1m radius; X [m]; Y [m]", 300,-1.5,1.5,300,-1.5,1.5);
    TH2F *h_TrackingPlane5m_XY_1mRadius = new TH2F("h_TrackingPlane5m_XY_1mRadius", "occupancy plot, 5 m from IP, within 1m radius; X [m]; Y [m]", 300,-1.5,1.5,300,-1.5,1.5);
    TH2F *h_TrackingPlane6m_XY_1mRadius = new TH2F("h_TrackingPlane6m_XY_1mRadius", "occupancy plot, 6 m from IP, within 1m radius; X [m]; Y [m]", 300,-1.5,1.5,300,-1.5,1.5);
    
    
    /// only within 1m radius but outside 0.02m radius, zoomed
    TH2F *h_TrackingPlane1m_XY_1mRadiusZoomed = new TH2F("h_TrackingPlane1m_XY_1mRadiusZoomed", "occupancy plot, 1 m from IP, outside of 0.02 m and within 1m radius; X [m]; Y [m]", 50,-0.05,0.05,50,-0.05,0.05);
    TH2F *h_TrackingPlane2m_XY_1mRadiusZoomed = new TH2F("h_TrackingPlane2m_XY_1mRadiusZoomed", "occupancy plot, 2 m from IP, outside of 0.02 m and within 1m radius; X [m]; Y [m]", 50,-0.05,0.05,50,-0.05,0.05);
    TH2F *h_TrackingPlane3m_XY_1mRadiusZoomed = new TH2F("h_TrackingPlane3m_XY_1mRadiusZoomed", "occupancy plot, 3 m from IP, outside of 0.02 m and within 1m radius; X [m]; Y [m]", 50,-0.05,0.05,50,-0.05,0.05);
    TH2F *h_TrackingPlane4m_XY_1mRadiusZoomed = new TH2F("h_TrackingPlane4m_XY_1mRadiusZoomed", "occupancy plot, 4 m from IP, outside of 0.02 m and within 1m radius; X [m]; Y [m]", 50,-0.05,0.05,50,-0.05,0.05);
    TH2F *h_TrackingPlane5m_XY_1mRadiusZoomed = new TH2F("h_TrackingPlane5m_XY_1mRadiusZoomed", "occupancy plot, 5 m from IP, outside of 0.02 m and within 1m radius; X [m]; Y [m]", 50,-0.05,0.05,50,-0.05,0.05);
    TH2F *h_TrackingPlane6m_XY_1mRadiusZoomed = new TH2F("h_TrackingPlane6m_XY_1mRadiusZoomed", "occupancy plot, 6 m from IP, outside of 0.02 m and within 1m radius; X [m]; Y [m]", 50,-0.05,0.05,50,-0.05,0.05);
    
    
   
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
              float inTheta = thetaCal(Particle_Eta[nParticle]);
              h_IncomingPhotonTheta->Fill(inTheta, weightHistogram);
              
          }
          /// select the outgoing photon
          /// photon without any cuts
          if(Particle_Status[nParticle] == 1 && Particle_PID[nParticle] == 22)h_OutgoingPhotonE_NoCut->Fill(Particle_E[nParticle], weightHistogram);
          
          if(Particle_Status[nParticle] == 1 && Particle_PID[nParticle] == 22 && Particle_E[nParticle] > photonThreshold){
              h_OutgoingPhotonE->Fill(Particle_E[nParticle], weightHistogram);
              h_OutgoingPhotonPT->Fill(Particle_PT[nParticle], weightHistogram);
              h_OutgoingPhotonEta->Fill(Particle_Eta[nParticle], weightHistogram);
              h_OutgoingPhotonPhi->Fill(Particle_Phi[nParticle], weightHistogram);
              float outTheta = thetaCal(Particle_Eta[nParticle]);
              h_OutgoingPhotonTheta->Fill(outTheta, weightHistogram);
              h_OutgoingPhotonE_Theta->Fill(outTheta,Particle_E[nParticle]);
              outgoingPhotonOrder.push_back(nParticle);
          }
       }
       /// only select events with two photons having the threshold cut, otherwise no need to go further
       if(outgoingPhotonOrder.size()<2)continue;
       //// add the parameters needed for ALP travel
       float mALP       = mass.Atof(); /// taking different ALP mass values
       float tauA       = getTau(mALP);
       float ctauA      = tauA*(0.1975)*pow(10,-15); // c is 1 in natural units
       float pA         = sqrt(pow(Particle_E[incomingPhotonOrder],2) - pow(mALP,2));
       float LA         = ctauA*pA/mALP;
//        std::cout << "The tauA: " << tauA << std::endl;
//        std::cout << "The pA/mALP is " << pA/mALP << std::endl;
//        std::cout << "ctauA : " << ctauA << std::endl; 
       std::cout << "The LA value for " << mALP << " GeV is " << LA << std::endl;
       float LS         = 0.5;
       float LD         = 5.5;
       double expFactor = exp(-(LS/LA)) - exp(-((LD+LS)/LA));
       
       TF1 *f1 = new TF1("f1","1/[0]*exp(-x/[0])",0.0,10000.0);
       f1->SetParameter(0,LA);
       f1->SetParName(0,"decay length");
       double r = 0.0; // this is the r where ALP is decaying
       /// getting the random number for vertex position r
       do{
         r = f1->GetRandom();  
       }
       while((r < LS) || (r > (LS+LD)));
       
       //cout << "The r value after: " << r << endl;
       h_rALP->Fill(r);
       ///  select two outgoing photons
       vector<TLorentzVector> photon;
       
       for(size_t i = 0; i < outgoingPhotonOrder.size(); ++i){
           TLorentzVector photonSample;
           photonSample.SetPtEtaPhiE(Particle_PT[outgoingPhotonOrder.at(i)], Particle_Eta[outgoingPhotonOrder.at(i)], Particle_Phi[outgoingPhotonOrder.at(i)], Particle_E[outgoingPhotonOrder.at(i)]);
           photon.push_back(photonSample);
       }
       
       /// get the opneing angle between two photons
       double openAng = openingAngle(photon.at(0), photon.at(1));
       h_OutgoingPhotonOpeningAngle->Fill(openAng, weightHistogram);
       
       /// get the properties of ALP by adding the two photons it decay into.
       TLorentzVector alp = photon.at(0) + photon.at(1);
       
       float alpEta   = alp.Eta();
       float alpTheta = thetaCal(alpEta); 
       
       h_ALPEnergy->Fill(alp.E(), weightHistogram);
       h_ALPTheta->Fill(alpTheta, weightHistogram);
       h_ALPThetaZoomed->Fill(alpTheta, weightHistogram);
       h_ALPMomentum->Fill(alp.P(), weightHistogram);
       h_ALPTransverseMomentum->Fill(alp.Pt(), weightHistogram);
       h_ALPMass->Fill(alp.M(), weightHistogram);
       h_ALPThetaVsEnergy->Fill(alp.E(), alpTheta, weightHistogram);
       
       /// make some plots from the outgoing photons
       double thetaLeadingE     = -999.;
       double thetaSubLeadingE  = -999.;
       double energyLeadingE    = -999.;
       double energySubLeadingE = -999.;
       
       if(Particle_E[outgoingPhotonOrder.at(0)] > Particle_E[outgoingPhotonOrder.at(1)]){
           thetaLeadingE    = thetaCal(Particle_Eta[outgoingPhotonOrder.at(0)]);
           thetaSubLeadingE = thetaCal(Particle_Eta[outgoingPhotonOrder.at(1)]);
           energyLeadingE   = Particle_E[outgoingPhotonOrder.at(0)];
           energySubLeadingE= Particle_E[outgoingPhotonOrder.at(1)];
       }
       else{
           thetaLeadingE    = thetaCal(Particle_Eta[outgoingPhotonOrder.at(1)]);
           thetaSubLeadingE = thetaCal(Particle_Eta[outgoingPhotonOrder.at(0)]);
           energyLeadingE   =  Particle_E[outgoingPhotonOrder.at(1)];
           energySubLeadingE=  Particle_E[outgoingPhotonOrder.at(0)];
       }
       
       h_OutgoingPhoton_Theta1VsTheta2->Fill(thetaLeadingE, thetaSubLeadingE, weightHistogram);
       h_OutgoingLeadingPhotonE_Theta->Fill(energyLeadingE, thetaLeadingE, weightHistogram);
       h_OutgoingSubLeadingPhotonE_Theta->Fill(energySubLeadingE, thetaSubLeadingE, weightHistogram);
       
       h_OutgoingPhoton_Theta1VsTheta2_WithALP->Fill(thetaLeadingE - alpTheta, thetaSubLeadingE - alpTheta, weightHistogram);
       if(thetaLeadingE<=5 && thetaSubLeadingE<=5)h_OutgoingPhoton_Theta1VsTheta2_WithCut->Fill(thetaLeadingE, thetaSubLeadingE, weightHistogram);
       h_OutgoingPhoton_diffAngleVsLeadingPhotonAngle->Fill(thetaLeadingE, std::abs(thetaLeadingE-thetaSubLeadingE), weightHistogram);
       h_OutgoingPhoton_Energy1VsEnergy2->Fill(energyLeadingE, energySubLeadingE, weightHistogram);
       h_OutgoingPhoton_diffEnergyVsLeadingPhotonEnergy->Fill(energyLeadingE, std::abs(energyLeadingE - energySubLeadingE), weightHistogram);
       h_OutgoingPhoton_diffALPEnergyVsLeadingPhotonEnergy->Fill(energyLeadingE, std::abs(energyLeadingE - alp.E()), weightHistogram);
       
       h_OutgoingPhoton_LeadingPhotonThetaVsLeadingPhotonEnergyFractionALP->Fill(energyLeadingE/alp.E(), thetaLeadingE, weightHistogram);
       h_OutgoingPhoton_LeadingPhotonThetaWrtALPVsLeadingPhotonEnergyFractionALP->Fill(energyLeadingE/alp.E(), std::abs(thetaLeadingE - alpTheta), weightHistogram);
       
       h_OutgoingPhoton_SubLeadingPhotonThetaVsSubLeadingPhotonEnergyFractionALP->Fill(energySubLeadingE/alp.E(), thetaSubLeadingE, weightHistogram);
       h_OutgoingPhoton_SubLeadingPhotonThetaWrtALPVsSubLeadingPhotonEnergyFractionALP->Fill(energySubLeadingE/alp.E(), std::abs(thetaSubLeadingE - alpTheta), weightHistogram);
       
       
       
       //// from here make the acceptance plots
       alpTheta = (TMath::Pi()/180.0)*alpTheta; // converting theta to radian
       
       
       /// block to get occupancy plots after each meter.
       /// get x,y,z of the ALP
       float xALP   = r*sin(alpTheta)*cos(alp.Phi());
       float yALP   = r*sin(alpTheta)*sin(alp.Phi());
       float zALP   = r*cos(alp.Phi());
       
       float x1m[2] = {-999.0,-999.0};
       float x2m[2] = {-999.0,-999.0};
       float x3m[2] = {-999.0,-999.0};
       float x4m[2] = {-999.0,-999.0};
       float x5m[2] = {-999.0,-999.0};
       float x6m[2] = {-999.0,-999.0};
       
       float y1m[2] = {-999.0,-999.0};
       float y2m[2] = {-999.0,-999.0};
       float y3m[2] = {-999.0,-999.0};
       float y4m[2] = {-999.0,-999.0};
       float y5m[2] = {-999.0,-999.0};
       float y6m[2] = {-999.0,-999.0};
       
       float theta[2] = {-999.0,-999.0};
       
       theta[0] = TMath::Pi()/180.0*thetaCal(Particle_Eta[outgoingPhotonOrder.at(0)]);
       theta[1] = TMath::Pi()/180.0*thetaCal(Particle_Eta[outgoingPhotonOrder.at(1)]);
       
       /// when the detector is at 1 m distance
       if((1 - zALP)>0){
          x1m[0] = xALP + (1 - zALP)*tan(theta[0])*cos(Particle_Phi[outgoingPhotonOrder.at(0)]); 
          x1m[1] = xALP + (1 - zALP)*tan(theta[1])*cos(Particle_Phi[outgoingPhotonOrder.at(1)]);
          y1m[0] = yALP + (1 - zALP)*tan(theta[0])*sin(Particle_Phi[outgoingPhotonOrder.at(0)]);
          y1m[1] = yALP + (1 - zALP)*tan(theta[1])*sin(Particle_Phi[outgoingPhotonOrder.at(1)]);
          
          for(int l=0;l<2;++l){
              h_TrackingPlane1m_XY->Fill(x1m[l], y1m[l], weightHistogram);
          }
            
          if((sqrt(pow(x1m[0],2)+pow(y1m[0],2))<=1.0) && (sqrt(pow(x1m[1],2)+pow(y1m[1],2))<=1.0)){
              for(int l=0;l<2;++l){
                  h_TrackingPlane1m_XY_1mRadius->Fill(x1m[l], y1m[l], weightHistogram);
                  h_TrackingPlane1m_XY_1mRadiusZoomed->Fill(x1m[l], y1m[l],weightHistogram);
              }
          }
       }
       
       /// when the detector is at 2m distance
       
       if((2 - zALP)>0){
          x2m[0] = xALP + (2 - zALP)*tan(theta[0])*cos(Particle_Phi[outgoingPhotonOrder.at(0)]); 
          x2m[1] = xALP + (2 - zALP)*tan(theta[1])*cos(Particle_Phi[outgoingPhotonOrder.at(1)]);
          y2m[0] = yALP + (2 - zALP)*tan(theta[0])*sin(Particle_Phi[outgoingPhotonOrder.at(0)]);
          y2m[1] = yALP + (2 - zALP)*tan(theta[1])*sin(Particle_Phi[outgoingPhotonOrder.at(1)]);
          
          for(int l=0;l<2;++l){
              h_TrackingPlane2m_XY->Fill(x2m[l], y2m[l], weightHistogram);
          }
            
          if((sqrt(pow(x2m[0],2)+pow(y2m[0],2))<=1.0) && (sqrt(pow(x2m[1],2)+pow(y2m[1],2))<=1.0)){
              for(int l=0;l<2;++l){
                  h_TrackingPlane2m_XY_1mRadius->Fill(x2m[l], y2m[l], weightHistogram);
                  h_TrackingPlane2m_XY_1mRadiusZoomed->Fill(x2m[l], y2m[l],weightHistogram);
              }
          }
       }
       
       /// when the detector is at 3m distance
       
       if((3 - zALP)>0){
          x3m[0] = xALP + (3 - zALP)*tan(theta[0])*cos(Particle_Phi[outgoingPhotonOrder.at(0)]); 
          x3m[1] = xALP + (3 - zALP)*tan(theta[1])*cos(Particle_Phi[outgoingPhotonOrder.at(1)]);
          y3m[0] = yALP + (3 - zALP)*tan(theta[0])*sin(Particle_Phi[outgoingPhotonOrder.at(0)]);
          y3m[1] = yALP + (3 - zALP)*tan(theta[1])*sin(Particle_Phi[outgoingPhotonOrder.at(1)]);
          
          for(int l=0;l<2;++l){
              h_TrackingPlane3m_XY->Fill(x3m[l], y3m[l], weightHistogram);
          }
            
          if((sqrt(pow(x3m[0],2)+pow(y3m[0],2))<=1.0) && (sqrt(pow(x3m[1],2)+pow(y3m[1],2))<=1.0)){
              for(int l=0;l<2;++l){
                  h_TrackingPlane3m_XY_1mRadius->Fill(x3m[l], y3m[l], weightHistogram);
                  h_TrackingPlane3m_XY_1mRadiusZoomed->Fill(x3m[l], y3m[l],weightHistogram);
              }
          }
       }
       
       /// when the detector is at 4m distance
       
       if((4 - zALP)>0){
          x4m[0] = xALP + (4 - zALP)*tan(theta[0])*cos(Particle_Phi[outgoingPhotonOrder.at(0)]); 
          x4m[1] = xALP + (4 - zALP)*tan(theta[1])*cos(Particle_Phi[outgoingPhotonOrder.at(1)]);
          y4m[0] = yALP + (4 - zALP)*tan(theta[0])*sin(Particle_Phi[outgoingPhotonOrder.at(0)]);
          y4m[1] = yALP + (4 - zALP)*tan(theta[1])*sin(Particle_Phi[outgoingPhotonOrder.at(1)]);
          
          for(int l=0;l<2;++l){
              h_TrackingPlane4m_XY->Fill(x4m[l], y4m[l], weightHistogram);
          }
            
          if((sqrt(pow(x4m[0],2)+pow(y4m[0],2))<=1.0) && (sqrt(pow(x4m[1],2)+pow(y4m[1],2))<=1.0)){
              for(int l=0;l<2;++l){
                  h_TrackingPlane4m_XY_1mRadius->Fill(x4m[l], y4m[l], weightHistogram);
                  h_TrackingPlane4m_XY_1mRadiusZoomed->Fill(x4m[l], y4m[l],weightHistogram);
              }
          }
       }
       
       /// when the detector is at 5m distance
       
       if((5 - zALP)>0){
          x5m[0] = xALP + (5 - zALP)*tan(theta[0])*cos(Particle_Phi[outgoingPhotonOrder.at(0)]); 
          x5m[1] = xALP + (5 - zALP)*tan(theta[1])*cos(Particle_Phi[outgoingPhotonOrder.at(1)]);
          y5m[0] = yALP + (5 - zALP)*tan(theta[0])*sin(Particle_Phi[outgoingPhotonOrder.at(0)]);
          y5m[1] = yALP + (5 - zALP)*tan(theta[1])*sin(Particle_Phi[outgoingPhotonOrder.at(1)]);
          
          for(int l=0;l<2;++l){
              h_TrackingPlane5m_XY->Fill(x5m[l], y5m[l], weightHistogram);
          }
            
          if((sqrt(pow(x5m[0],2)+pow(y5m[0],2))<=1.0) && (sqrt(pow(x5m[1],2)+pow(y5m[1],2))<=1.0)){
              for(int l=0;l<2;++l){
                  h_TrackingPlane5m_XY_1mRadius->Fill(x5m[l], y5m[l], weightHistogram);
                  h_TrackingPlane5m_XY_1mRadiusZoomed->Fill(x5m[l], y5m[l],weightHistogram);
              }
          }
       }
       
       /// when the detector is at 6m distance
       
       if((6 - zALP)>0){
          x6m[0] = xALP + (6 - zALP)*tan(theta[0])*cos(Particle_Phi[outgoingPhotonOrder.at(0)]); 
          x6m[1] = xALP + (6 - zALP)*tan(theta[1])*cos(Particle_Phi[outgoingPhotonOrder.at(1)]);
          y6m[0] = yALP + (6 - zALP)*tan(theta[0])*sin(Particle_Phi[outgoingPhotonOrder.at(0)]);
          y6m[1] = yALP + (6 - zALP)*tan(theta[1])*sin(Particle_Phi[outgoingPhotonOrder.at(1)]);
          
          for(int l=0;l<2;++l){
              h_TrackingPlane6m_XY->Fill(x6m[l], y6m[l], weightHistogram);
          }
            
          if((sqrt(pow(x6m[0],2)+pow(y6m[0],2))<=1.0) && (sqrt(pow(x6m[1],2)+pow(y6m[1],2))<=1.0)){
              for(int l=0;l<2;++l){
                  h_TrackingPlane6m_XY_1mRadius->Fill(x6m[l], y6m[l], weightHistogram);
                  h_TrackingPlane6m_XY_1mRadiusZoomed->Fill(x6m[l], y6m[l],weightHistogram);
              }
          }
       }
       
    }/// loop over all events end
    
    
    fOut->Write();
    fOut->Close();
    // Record end time
    auto finish = std::chrono::steady_clock::now();
    auto diff = finish - start;
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time : " << chrono::duration <double, milli> (diff).count()/1000.0 << " s" << endl;
   
}
