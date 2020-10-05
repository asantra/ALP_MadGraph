import os, sys, time, glob
import math, argparse
import array
from ROOT import *

def main():
    
    parser   = argparse.ArgumentParser()
    parser.add_argument('-m', action="store", dest="massValue")
    parser.add_argument('-l', action="store", dest="fileValue", help="the name of the filenames")
    args = parser.parse_args()
    
    fOut   = open("outputPhotonNumbers.txt", "w")
    fOut.write("#energy\txsecv3 (m=0.5 GeV)\txsecv4 (m=0.001GeV)\txsecv5 (m=0.01GeV)\txsecv6 (m=0.05GeV)\tnPhoton\n")
    Ne     = 6*pow(10,9)
    tPulse = 35*pow(10,-15)
    rootOut = TFile("photonEnergy.root","RECREATE")
    rootOut.cd()
    h_InputPhotonEnergy  = TH1F("h_InputPhotonEnergy", "Input photon spectrum; E [GeV]; N",2000,0,0.0001)
    numberOfEnergyPoints = int((10.6 - float(args.massValue))/0.1)
    x = array.array('d')
    y = array.array('d')
    
    totalDwDE = 0
    
    with open(args.fileValue) as inFile:
        for lines in inFile.readlines():
            if "#" in lines: continue
            allWords    = lines.split()
            energy      = allWords[0]
            floatEnergy = float(energy)
            if(floatEnergy > 10.6):
                continue
            
            energyForXsec = round(floatEnergy - 0.001, 1)
            
            xsecValue3 = "0.0"
            xsecValue4 = "0.0"
            xsecValue5 = "0.0"
            xsecValue6 = "0.0"

            xsecFileName3  = "run_01_tag_1_banner_0.5GeV_beam"+str(energyForXsec)+"GeV_v3.txt"
            xsecFileName4  = "run_01_tag_1_banner_0.001GeV_beam"+str(energyForXsec)+"GeV_v1.txt"
            xsecFileName5  = "run_01_tag_1_banner_0.01GeV_beam"+str(energyForXsec)+"GeV_v1.txt"
            xsecFileName6  = "run_01_tag_1_banner_0.05GeV_beam"+str(energyForXsec)+"GeV_v1.txt"
                
            ### this is for 0.5 GeV of ALP mass, a special case as it does not have files for energy < 0.5 GeV
            if os.path.exists("/Users/arkasantra/MadGraphWorks/LHEFiles/"+xsecFileName3):
                fInXsec3       = open("/Users/arkasantra/MadGraphWorks/LHEFiles/"+xsecFileName3)
                for lineXsec3 in fInXsec3.readlines():
                    lineXsec3 = lineXsec3.rstrip()
                    if "#  Integrated weight (pb)  :" not in lineXsec3:
                        continue
                    xsecValue3 = lineXsec3.split(":")[1]
            else:
                xsecValue3 = "0.0"
                
                
            if os.path.exists("/Users/arkasantra/MadGraphWorks/LHEFiles/"+xsecFileName4):
                fInXsec4       = open("/Users/arkasantra/MadGraphWorks/LHEFiles/"+xsecFileName4)
                for lineXsec4 in fInXsec4.readlines():
                    lineXsec4 = lineXsec4.rstrip()
                    if "#  Integrated weight (pb)  :" not in lineXsec4:
                        continue
                    xsecValue4 = lineXsec4.split(":")[1]
            else:
                xsecValue4 = "0.0"
                
                
            if os.path.exists("/Users/arkasantra/MadGraphWorks/LHEFiles/"+xsecFileName5):
                fInXsec5       = open("/Users/arkasantra/MadGraphWorks/LHEFiles/"+xsecFileName5)
                for lineXsec5 in fInXsec5.readlines():
                    lineXsec5 = lineXsec5.rstrip()
                    if "#  Integrated weight (pb)  :" not in lineXsec5:
                        continue
                    xsecValue5 = lineXsec5.split(":")[1]
            else:
                xsecValue5 = "0.0"
                
                
                
            if os.path.exists("/Users/arkasantra/MadGraphWorks/LHEFiles/"+xsecFileName6):
                fInXsec6       = open("/Users/arkasantra/MadGraphWorks/LHEFiles/"+xsecFileName6)
                for lineXsec6 in fInXsec6.readlines():
                    lineXsec6 = lineXsec6.rstrip()
                    if "#  Integrated weight (pb)  :" not in lineXsec6:
                        continue
                    xsecValue6 = lineXsec6.split(":")[1]
            else:
                xsecValue6 = "0.0"
                
                
                
            dwde      = float(allWords[4])  ### for xi = 2.0
            nPhoton   = Ne*dwde*10
            
            print "energy: ", energy, " dwde: ", dwde
            totalDwDE += dwde*0.1
            
            
            
            h_InputPhotonEnergy.Fill(nPhoton)
            x.append(float(energyForXsec))
            y.append(nPhoton)
            
            fOut.write(energy+"\t"+xsecValue3+"\t"+xsecValue4+"\t"+xsecValue5+"\t"+xsecValue6+"\t"+str(nPhoton)+"\n")
            
    print "normalization of dwde: ", totalDwDE
            
    g_PhotonDistribution = TGraph(numberOfEnergyPoints, x, y)
    g_PhotonDistribution.Write("g_PhotonDistribution")
    rootOut.Write()
    rootOut.Close()
    fOut.close()


if __name__=="__main__":
  main()
        
