/* This is the script to add many root files into one file according to the weight */
/* This is based on RDataFrame class of root */
/* Usage: inside root: .L combineRootFiles.C++; combineRootFiles()  */
/* written by Arka Santra, arka.santra@cern.ch */

#include <iostream>
#include "math.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TCollection.h"
#include "TKey.h"
#include "TClass.h"
#include "TMath.h"
#include <string>
#include <sys/stat.h>
#include "TChain.h"
#include "TTree.h"
#include <ROOT/RDataFrame.hxx>
#include <algorithm>
#include <time.h>
#include <chrono>  // for high_resolution_clock
#include <sstream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <limits>
#include <math.h>
#include <vector>
#include <tuple>

using namespace std;
using namespace ROOT;
using ints     = ROOT::VecOps::RVec<int>;
using floats   = ROOT::VecOps::RVec<float>;
using doubles  = ROOT::VecOps::RVec<double>;
bool debug = false;

int ProcessList(const std::string &fnamelist, std::vector<std::string> &flist);

void combineRootFiles(const char *fnlist = 0){
    auto start = std::chrono::steady_clock::now();
    string fileOut  =  "combinedRootFile_AllSamples.root"; //
    string treeOut  =  "outputTree";
    
    if (!fnlist){
       std::cout << "Usage: root -l combineRootFiles.C\'(\"file_with_list_of_mc_files\")\'\n";
       return;
    }
  
    std::string fnamelist(fnlist);
    std::vector<std::string>  flist;
    ProcessList(fnamelist, flist);  
    if (debug) {
       std::cout << "The following files will be processed:\n";
       std::for_each(flist.begin(), flist.end(), [](const std::string ss) {std::cout << ss << std::endl;});
    }
    
    TChain *LHEFtree = new TChain("LHEF");
    std::for_each(flist.begin(), flist.end(), [LHEFtree](const std::string ss) {LHEFtree->Add(ss.c_str());} );
    
    ROOT::RDataFrame d(*LHEFtree);
    
    auto count = d.Count();
    //# Determine the number of events to loop over
    unsigned long long rangeNumber = -1;
    
    rangeNumber = *count;
    std::cout << "The number of events processing: " << rangeNumber << std::endl;
    
    float photonEnergy = 2.2;
    int   nEvents      = 30000;
    
    auto CalcWeight = [](float photoE, int nEv){
        auto eventBeamWeight = photoE/nEv;
        return eventBeamWeight;
    };
    
    vector<string> colNames;
    colNames.push_back("photonE");
    colNames.push_back("nEvGenerated");
    colNames.push_back("correctWeights");
    colNames.push_back("event");
//     colNames.push_back("eventUniqueId");
//     colNames.push_back("eventfBits");
//     colNames.push_back("eventNumber");
//     colNames.push_back("eventNparticles");
//     colNames.push_back("eventProcessID");
//     colNames.push_back("eventWeight");
//     colNames.push_back("eventScalePDF");
//     colNames.push_back("eventCouplingQED");
//     colNames.push_back("eventCouplingQCD");
//     colNames.push_back("eventSize");
//     colNames.push_back("rwgtSize");
    colNames.push_back("rwgt");
    colNames.push_back("particle");
//     colNames.push_back("particlefUniqueID");
//     colNames.push_back("particlefBits");
//     colNames.push_back("particlePID");
//     colNames.push_back("particleStatus");
//     colNames.push_back("particleMother1");
//     colNames.push_back("particleMother2");
//     colNames.push_back("particleColorLine1");
//     colNames.push_back("particleColorLine2");
//     colNames.push_back("particlePx");
//     colNames.push_back("particlePy");
//     colNames.push_back("particlePz");
//     //colNames.push_back("particleE");
//     colNames.push_back("particleM");
//     colNames.push_back("particlePt");
//     colNames.push_back("particleEta");
//     colNames.push_back("particlePhi");
//     colNames.push_back("particleRapidity");
//     colNames.push_back("particleLifeTime");
//     colNames.push_back("particleSpin");
//     colNames.push_back("particleSize");
    
    auto dCut   = d.Define("photonE", [photonEnergy] { return photonEnergy; })
                   .Define("nEvGenerated", [nEvents] { return nEvents; })
                   .Define("correctWeights", CalcWeight, {"photonE", "nEvGenerated"})
                   .Define("event", "Event")
//                    .Define("eventUniqueId", "Event.fUniqueID")
//                    .Define("eventfBits", "Event.fBits")
//                    .Define("eventNumber", "Event.Number")
//                    .Define("eventNparticles", "Event.Nparticles")
//                    .Define("eventProcessID", "Event.ProcessID")
//                    .Define("eventWeight", "Event.Weight")
//                    .Define("eventScalePDF", "Event.ScalePDF")
//                    .Define("eventCouplingQED", "Event.CouplingQED")
//                    .Define("eventCouplingQCD", "Event.CouplingQCD")
//                    .Define("eventSize", "Event_size")
//                    .Define("rwgtSize", "Rwgt_size")
                   .Define("rwgt", "Rwgt")
                   .Define("particle", "Particle")
//                    .Define("particlefUniqueID", "Particle.fUniqueID")
//                    .Define("particlefBits", "Particle.fBits")
//                    .Define("particlePID","Particle.PID")
//                    .Define("particleStatus", "Particle.Status")
//                    .Define("particleMother1", "Particle.Mother1")
//                    .Define("particleMother2", "Particle.Mother2")
//                    .Define("particleColorLine1", "Particle.ColorLine1")
//                    .Define("particleColorLine2", "Particle.ColorLine2")
//                    .Define("particlePx", "Particle.Px")
//                    .Define("particlePy", "Particle.Py")
//                    .Define("particlePz", "Particle.Pz")
//                    //.Define("particleE", "Particle.E")
//                    //.Define("particleE", "Particle.E")
//                    .Define("particleM", "Particle.M")
//                    .Define("particlePt", "Particle.PT")
//                    .Define("particleEta", "Particle.Eta")
//                    .Define("particlePhi", "Particle.Phi")
//                    .Define("particleRapidity", "Particle.Rapidity")
//                    .Define("particleLifeTime", "Particle.LifeTime")
//                    .Define("particleSpin", "Particle.Spin")
//                    .Define("particleSize", "Particle_size")
                   .Snapshot(treeOut, fileOut, colNames);
                   
    // Record end time
    auto finish = std::chrono::steady_clock::now();
    auto diff = finish - start;
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time : " << chrono::duration <double, milli> (diff).count()/1000.0 << " s" << endl;
}


int ProcessList(const std::string &fnamelist, std::vector<std::string> &flist){
    std::fstream  fdata;
    fdata.open(fnamelist, std::ios::in);
    if (!fdata.is_open()) {
        throw std::runtime_error(std::string("Error reading data from the file ") + fnamelist);
    }
  
    unsigned long lid = 0;
    while (!fdata.eof()) {
        std::string  ffname;
        double fweight;
        fdata >> ffname;
        if (!fdata.fail()) { 
//         std::cout << "File name " << ffname << " is read from the list file" << std::endl;
           flist.push_back(ffname);
        }
        else if (fdata.eof()) { break; }
        else {
            std::cout << "ProcessList(..)  :  Error reading data from the file " << fnamelist 
                << ",  line: " << lid << ". Exit." << std::endl;
            fdata.close();          
            return -2;
        }
        ++lid;
    }
    fdata.close();

    return 0;
}
