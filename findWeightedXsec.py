import os, sys, time, glob
import argparse


def main():
    
    parser   = argparse.ArgumentParser()
    parser.add_argument('-l', action="store", dest="fileValue", help="the name of the filenames")
    args = parser.parse_args()
    
    with open(args.fileValue) as inFile:
        totXsecM0p5     = 0
        totXsecM0p001   = 0
        totXsecM0p01    = 0
        totXsecM0p05    = 0
        
        totWeightM0p5   = 0
        totWeightM0p001 = 0
        totWeightM0p01  = 0
        totWeightM0p05  = 0
        
        
        for lines in inFile.readlines():
            if "#" in lines: continue
            lines        = lines.rstrip()
            eachWord     = lines.split()
            
            ### getting the xsec
            
            if("OLD" in args.fileValue):
                xsecM0p5     = float(eachWord[5])
                xsecM0p001   = float(eachWord[6])
                xsecM0p01    = float(eachWord[7])
                xsecM0p05    = float(eachWord[8])
                nPhotons     = float(eachWord[9])
            else:
                xsecM0p5     = float(eachWord[1])
                xsecM0p001   = float(eachWord[2])
                xsecM0p01    = float(eachWord[3])
                xsecM0p05    = float(eachWord[4])
                nPhotons     = float(eachWord[5])
            
            ## getting the weights
            weightM0p5   = nPhotons/30000.0
            weightM0p001 = nPhotons/30000.0
            weightM0p01  = nPhotons/30000.0
            weightM0p05  = nPhotons/30000.0
            
            ### M = 0.001 GeV
            totXsecM0p001   += xsecM0p001*weightM0p001
            if(xsecM0p001==0.0):
                totWeightM0p001   += 0
            else:
                totWeightM0p001   += weightM0p001
            
            ### M = 0.01 GeV
            totXsecM0p01    += xsecM0p01*weightM0p01
            if(xsecM0p01==0.0):
                totWeightM0p01   += 0
            else:
                totWeightM0p01   += weightM0p01
            
            ### M = 0.05 GeV
            totXsecM0p05    += xsecM0p05*weightM0p05
            if(xsecM0p05==0.0):
                totWeightM0p05   += 0
            else:
                totWeightM0p05   += weightM0p05
            
            ### M = 0.5 GeV
            totXsecM0p5     += xsecM0p5*weightM0p5
            if(xsecM0p5==0.0):
                totWeightM0p5   += 0
            else:
                totWeightM0p5   += weightM0p5
                
            #print "totXsecM0p5: ", totXsecM0p5, " totWeightM0p5: ", totWeightM0p5
            
            
        weightedXsecM0p001 = (totXsecM0p001 / totWeightM0p001)*0.7228041518  ### different weighting factor for mass 0.001 GeV
        weightedXsecM0p01  = (totXsecM0p01 / totWeightM0p01)*0.6897650549  ### weighting factor from mass 0.1 GeV
        weightedXsecM0p05  = (totXsecM0p05 / totWeightM0p05)*0.6897650549 ### weighting factor from mass 0.1 GeV
        weightedXsecM0p5   = (totXsecM0p5 / totWeightM0p5)*0.592532760524  ### weighting factor from mass 0.5 GeV
        
        print "totXsecM0p001/totWeightM0p001 : ", totXsecM0p001, " / ", totWeightM0p001, " = ", weightedXsecM0p001
        print "totXsecM0p01/totWeightM0p01   : ", totXsecM0p01, " / ", totWeightM0p01, " = ", weightedXsecM0p01
        print "totXsecM0p05/totWeightM0p05   : ", totXsecM0p05, " / ", totWeightM0p05, " = ", weightedXsecM0p05
        print "totXsecM0p5/totWeightM0p5     : ", totXsecM0p5, " / ", totWeightM0p5, " = ", weightedXsecM0p5
        
        
        print "xsec for ALP mass 0.001 GeV   : ", weightedXsecM0p001, " pb"
        print "xsec for ALP mass 0.01 GeV    : ", weightedXsecM0p01, " pb"
        print "xsec for ALP mass 0.05 GeV    : ", weightedXsecM0p05, " pb"
        print "xsec for ALP mass 0.5 GeV     : ", weightedXsecM0p5, " pb"
            
            
            
            
            
            
            
if __name__=="__main__":
    main()
