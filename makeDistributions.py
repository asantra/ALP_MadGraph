import os, sys
import time, argparse
from ROOT import *
sys.path.insert(0, '/Users/arkasantra/arka/include')
from Functions import *

    
def main():
    gROOT.SetBatch()
    parser   = argparse.ArgumentParser()
    parser.add_argument('-m', action="store", dest="massValue")
    parser.add_argument('-v', action="store", dest="versionValue")
    parser.add_argument('-c', action="store", dest="cutValue", default="0.0")
    args = parser.parse_args()
    fIn = TFile("combinedOutputFile_ALPMass"+args.massValue+"GeV_"+args.cutValue+"GeVPhotonCut_V"+args.versionValue+".root")
    
    if(args.cutValue == "0.5"):
        if(args.versionValue == "2"):
            directory = "plotsPhotonCutV2"
        elif(args.versionValue == "4"):
            directory = "plotsPhotonCutV4"
        else:
            directory = "plotsPhotonCut"
            
    else:
        if(args.versionValue == "2"):
            directory = "plotsV2"
        elif(args.versionValue == "4"):
            directory = "plotsV4"
        else:
            directory = "plots"
        
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    h_TrackingPlane_XY                       = []
    h_TrackingPlane_XY_1mRadius              = []
    h_TrackingPlane_XY_1mRadiusZoomed        = []
    
    for i in xrange(1,7):
        #/// all around the plane
        hName                         = fIn.Get("h_TrackingPlane"+str(i)+"m_XY")
        hNameXY_1mRadius              = fIn.Get("h_TrackingPlane"+str(i)+"m_XY_1mRadius")
        hNameXY_1mRadiusZoomed        = fIn.Get("h_TrackingPlane"+str(i)+"m_XY_1mRadiusZoomed")
        
        h_TrackingPlane_XY.append(hName)
        h_TrackingPlane_XY_1mRadius.append(hNameXY_1mRadius)
        h_TrackingPlane_XY_1mRadiusZoomed.append(hNameXY_1mRadiusZoomed)
        
        
    h_OutgoingPhotonE_NoCut = fIn.Get("h_OutgoingPhotonE_NoCut")
        
    #### make acceptance calculations
    
    fTextOut = open(directory+"/acceptanceNumbers_"+args.massValue+"GeV_V"+args.versionValue+".txt", "w")
    
    fTextOut.write("%%% mass value: "+args.massValue+" GeV, weighted numbers\n")
    fTextOut.write("distance from ALPDecay & all photons & in 1m & acceptance\n")
    
    
    
    for i in xrange(0,6):
        #### fill up the acceptance text file 
        allPhotonIntegral        = h_OutgoingPhotonE_NoCut.Integral(0,501)
        photonin1mIntegral       = h_TrackingPlane_XY_1mRadius[i].Integral(0,301,0,301)
        
        acceptance1m             = round(photonin1mIntegral/float(allPhotonIntegral),3)
        
        fTextOut.write(str(i+1)+" & "+str("{:.3e}".format(allPhotonIntegral))+" & "+str("{:.3e}".format(photonin1mIntegral))+" & "+str(acceptance1m)+"\\\ \n")
        
        
        LegendName   = ["occupancy plots"]
        xAxisName    = "X [m]"
        yAxisName    = "Y [m]"
        zAxisName    = "Z [weighted events]"
        xrange1down  = -50.0
        xrange1up    = 50.0
        yrange1down  = -50.0
        yrange1up    = 50.0
        zrange1down  = 1e4
        zrange1up    = 1e9
        logz         = True
        
        LatexName    = str(i+1)+"m from ALP decay point"
        LatexName2   = "ALP mass: "+args.massValue+" GeV"
        FirstTH1     = [h_TrackingPlane_XY[i]]
        Draw2DHists(FirstTH1, LegendName, xAxisName, yAxisName, zAxisName, xrange1down, xrange1up, yrange1down, yrange1up, zrange1down, zrange1up, directory+"/"+FirstTH1[0].GetName()+"_"+args.massValue+"GeV", logz, LatexName)
    
    
        ##### only around 3m radius
        xrange1down  = -1.5
        xrange1up    = 1.5
        yrange1down  = -1.5
        yrange1up    = 1.5
        zrange1down  = 1e4
        zrange1up    = 1e9
        LatexName    = str(i+1)+"m from ALP decay point, ALP mass: "+args.massValue+" GeV"
        
        LatexName2   = "within 1m radius"
        FirstTH1     = [h_TrackingPlane_XY_1mRadius[i]]
        Draw2DHists(FirstTH1, LegendName, xAxisName, yAxisName, zAxisName, xrange1down, xrange1up, yrange1down, yrange1up, zrange1down, zrange1up, directory+"/"+FirstTH1[0].GetName()+"_"+args.massValue+"GeV", logz, LatexName, LatexName2)
        
        
        #### within 3 m radius and outside 0.02 m radius, zoomed
        xrange1down  = -0.05
        xrange1up    = 0.05
        yrange1down  = -0.05
        yrange1up    = 0.05
        zrange1down  = 1e4
        zrange1up    = 1e8
        
        LatexName2   = ""
        FirstTH1     = [h_TrackingPlane_XY_1mRadiusZoomed[i]]
        Draw2DHists(FirstTH1, LegendName, xAxisName, yAxisName, zAxisName, xrange1down, xrange1up, yrange1down, yrange1up, zrange1down, zrange1up, directory+"/"+FirstTH1[0].GetName()+"_"+args.massValue+"GeV", logz, LatexName, LatexName2)
        
        
        
    ### more plots
    h_OutgoingPhotonPhi                                 = fIn.Get("h_OutgoingPhotonPhi")
    h_OutgoingPhoton_Theta1VsTheta2                     = fIn.Get("h_OutgoingPhoton_Theta1VsTheta2")
    h_OutgoingPhoton_Theta1VsTheta2_WithCut             = fIn.Get("h_OutgoingPhoton_Theta1VsTheta2_WithCut")
    h_OutgoingPhoton_diffEnergyVsLeadingPhotonEnergy    = fIn.Get("h_OutgoingPhoton_diffEnergyVsLeadingPhotonEnergy")
    h_ALPTheta                                          = fIn.Get("h_ALPTheta")
    h_ALPEnergy                                         = fIn.Get("h_ALPEnergy")
    h_ALPThetaZoomed                                    = fIn.Get("h_ALPThetaZoomed")
    h_ALPThetaVsEnergy                                  = fIn.Get("h_ALPThetaVsEnergy")
    h_OutgoingLeadingPhotonE_Theta                      = fIn.Get("h_OutgoingLeadingPhotonE_Theta")
    h_OutgoingSubLeadingPhotonE_Theta                   = fIn.Get("h_OutgoingSubLeadingPhotonE_Theta")
    
    ## new plots 
    h_OutgoingPhoton_diffALPEnergyVsLeadingPhotonEnergy  = fIn.Get("h_OutgoingPhoton_diffALPEnergyVsLeadingPhotonEnergy")
    h_OutgoingPhoton_LeadingPhotonThetaVsLeadingPhotonEnergyFractionALP  = fIn.Get("h_OutgoingPhoton_LeadingPhotonThetaVsLeadingPhotonEnergyFractionALP")
    h_OutgoingPhoton_LeadingPhotonThetaWrtALPVsLeadingPhotonEnergyFractionALP  = fIn.Get("h_OutgoingPhoton_LeadingPhotonThetaWrtALPVsLeadingPhotonEnergyFractionALP")
    h_OutgoingPhoton_SubLeadingPhotonThetaVsSubLeadingPhotonEnergyFractionALP  = fIn.Get("h_OutgoingPhoton_SubLeadingPhotonThetaVsSubLeadingPhotonEnergyFractionALP")
    h_OutgoingPhoton_SubLeadingPhotonThetaWrtALPVsSubLeadingPhotonEnergyFractionALP  = fIn.Get("h_OutgoingPhoton_SubLeadingPhotonThetaWrtALPVsSubLeadingPhotonEnergyFractionALP")
        
        
    ##### only around 3m radius
    LatexName2   = ""
    
    LatexName    = "ALP mass: "+args.massValue+" GeV"
    FirstTH1     = [h_OutgoingPhoton_Theta1VsTheta2]
    Draw2DHists(FirstTH1, LegendName, "#theta_{1} [degrees]", "#theta_{2} [degrees]", "weighted events", 0.0, 180.0, 0.0, 180.0, 1e4, 1e10, directory+"/"+FirstTH1[0].GetName()+"_"+args.massValue+"GeV", logz, LatexName, LatexName2)
    
    FirstTH1     = [h_OutgoingPhoton_Theta1VsTheta2_WithCut]
    Draw2DHists(FirstTH1, LegendName, "#theta_{1} [degrees]", "#theta_{2} [degrees]", "weighted events", 0.0, 15.0, 0.0, 15.0, 1e3, 1e8, directory+"/"+FirstTH1[0].GetName()+"_"+args.massValue+"GeV", logz, LatexName, LatexName2)
    
    FirstTH1     = [h_OutgoingPhoton_diffEnergyVsLeadingPhotonEnergy]
    Draw2DHists(FirstTH1, LegendName, "leading photon energy [GeV]", "#Delta E [GeV]", "weighted events", 0.0, 10.5, 0.0, 10.5, 1e4, 1e10, directory+"/"+FirstTH1[0].GetName()+"_"+args.massValue+"GeV", logz, LatexName, LatexName2)
    
    
    FirstTH1     = [h_OutgoingPhoton_diffALPEnergyVsLeadingPhotonEnergy]
    Draw2DHists(FirstTH1, LegendName, "leading photon energy [GeV]", "#Delta E with ALP [GeV]", "weighted events", 0.0, 10.5, 0.0, 10.5, 1e4, 1e10, directory+"/"+FirstTH1[0].GetName()+"_"+args.massValue+"GeV", logz, LatexName, LatexName2)
    
    FirstTH1     = [h_OutgoingPhoton_LeadingPhotonThetaVsLeadingPhotonEnergyFractionALP]
    Draw2DHists(FirstTH1, LegendName, "leading photon energy fraction with ALP", "Leading photon #theta [degrees]", "weighted events", 0.4, 1, -180.0, 180.0, 1e4, 1e10, directory+"/"+FirstTH1[0].GetName()+"_"+args.massValue+"GeV", logz, LatexName, LatexName2)
    
    FirstTH1     = [h_OutgoingPhoton_LeadingPhotonThetaWrtALPVsLeadingPhotonEnergyFractionALP]
    Draw2DHists(FirstTH1, LegendName, "leading photon energy fraction with ALP", "Leading photon #theta wrt ALP [degrees]", "weighted events", 0.4, 1, -180.0, 180.0, 1e4, 1e10, directory+"/"+FirstTH1[0].GetName()+"_"+args.massValue+"GeV", logz, LatexName, LatexName2)
    
    
    FirstTH1     = [h_OutgoingPhoton_SubLeadingPhotonThetaVsSubLeadingPhotonEnergyFractionALP]
    Draw2DHists(FirstTH1, LegendName, "subleading photon energy fraction with ALP", "subleading photon #theta [degrees]", "weighted events", 0.0, 0.6, -180.0, 180.0, 1e4, 1e10, directory+"/"+FirstTH1[0].GetName()+"_"+args.massValue+"GeV", logz, LatexName, LatexName2)
    
    FirstTH1     = [h_OutgoingPhoton_SubLeadingPhotonThetaWrtALPVsSubLeadingPhotonEnergyFractionALP]
    Draw2DHists(FirstTH1, LegendName, "subleading photon energy fraction with ALP", "subleading photon #theta wrt ALP [degrees]", "weighted events", 0.0, 0.6, -180.0, 180.0, 1e4, 1e10, directory+"/"+FirstTH1[0].GetName()+"_"+args.massValue+"GeV", logz, LatexName, LatexName2)
    
    FirstTH1     = [h_OutgoingLeadingPhotonE_Theta]
    Draw2DHists(FirstTH1, LegendName, "leading photon energy", "leading photon #theta [degrees]", "weighted events", 0.0, 10.5, 0.0, 180.0, 1e4, 1e10, directory+"/"+FirstTH1[0].GetName()+"_"+args.massValue+"GeV", logz, LatexName, LatexName2)
    
    FirstTH1     = [h_OutgoingSubLeadingPhotonE_Theta]
    Draw2DHists(FirstTH1, LegendName, "subleading photon energy", "subleading photon #theta [degrees]", "weighted events", 0.0, 10.5, 0.0, 180.0, 1e4, 1e10, directory+"/"+FirstTH1[0].GetName()+"_"+args.massValue+"GeV", logz, LatexName, LatexName2)
    
    FirstTH1     = [h_ALPThetaVsEnergy]
    Draw2DHists(FirstTH1, LegendName, "ALP energy", "ALP #theta [degrees]", "weighted events", 0.0, 10.5, 0.0, 180.0, 1e4, 1e10, directory+"/"+FirstTH1[0].GetName()+"_"+args.massValue+"GeV", logz, LatexName, LatexName2)
    
    
    
    ### 1d plots
    LatexName = "ALP mass "+args.massValue+" GeV"
    FirstTH1 = [h_OutgoingPhotonPhi]
    
    DrawHists(FirstTH1, ["outgoing photon"], [4], "#phi [radian]", "weighted events", -3.2, 3.2, 1e4, 1e12, directory+"/"+FirstTH1[0].GetName()+"_"+args.massValue+"GeV", 1, 1, False, True, LatexName)
    
    FirstTH1 = [h_ALPTheta]
    
    DrawHists(FirstTH1, ["ALP"], [4], "#theta [degrees]", "weighted events", 0.0, 180.0, 1e4, 1e12, directory+"/"+FirstTH1[0].GetName()+"_"+args.massValue+"GeV", 1, 1, False, True, LatexName)
    
    FirstTH1 = [h_ALPEnergy]
    
    DrawHists(FirstTH1, ["ALP"], [4], "Energy [GeV]", "weighted events", 0.0, 11.0, 1e4, 1e12, directory+"/"+FirstTH1[0].GetName()+"_"+args.massValue+"GeV", 1, 1, False, True, LatexName)
    
    FirstTH1 = [h_ALPThetaZoomed]
    
    DrawHists(FirstTH1, ["ALP"], [4], "#theta [degrees]", "weighted events", 0.0, 10.0, 1e6, 1e12, directory+"/"+FirstTH1[0].GetName()+"_"+args.massValue+"GeV", 1, 1, False, True, LatexName)
    
    
    fTextOut.write("\n\n%%% mass value: "+args.massValue+" GeV, unweighted numbers\n")
    fTextOut.write("distance from ALPDecay & all photons & in 1m & acceptance\n")
    for i in xrange(0,6):
        allPhotonEntries        = h_OutgoingPhotonE_NoCut.GetEntries()
        photonin1mEntries       = h_TrackingPlane_XY_1mRadius[i].GetEntries()
        acceptance1m            = round(photonin1mEntries/float(allPhotonEntries),3)
        
        fTextOut.write(str(i+1)+" & "+str("{:.3e}".format(allPhotonEntries))+" & "+str("{:.3e}".format(photonin1mEntries))+" & "+str(acceptance1m)+"\\\ \n")
        
        
    fTextOut.close()
        
        
if __name__=="__main__":
    start = time.time()
    main()
    print "The processing time: ", time.time() - start, " s"
