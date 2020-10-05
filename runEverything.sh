#! /bin/bash
#### make the histograms, version 2
bash runWeightedRootFiles.sh 0.001 2 0.5
bash runWeightedRootFiles.sh 0.001 2 0.0
bash runWeightedRootFiles.sh 0.01 2 0.5
bash runWeightedRootFiles.sh 0.01 2 0.0
bash runWeightedRootFiles.sh 0.05 2 0.5
bash runWeightedRootFiles.sh 0.05 2 0.0
bash runWeightedRootFiles.sh 0.5 2 0.5
bash runWeightedRootFiles.sh 0.5 2 0.0

#### make the histograms, version 4
bash runWeightedRootFiles.sh 0.001 4 0.5
bash runWeightedRootFiles.sh 0.001 4 0.0
bash runWeightedRootFiles.sh 0.01 4 0.5
bash runWeightedRootFiles.sh 0.01 4 0.0
bash runWeightedRootFiles.sh 0.05 4 0.5
bash runWeightedRootFiles.sh 0.05 4 0.0
bash runWeightedRootFiles.sh 0.5 4 0.5
bash runWeightedRootFiles.sh 0.5 4 0.0

### plot the histograms, version 2
python makeDistributions.py -m 0.001 -v 2 -c 0.5
python makeDistributions.py -m 0.001 -v 2 -c 0.0
python makeDistributions.py -m 0.01 -v 2 -c 0.5
python makeDistributions.py -m 0.01 -v 2 -c 0.0
python makeDistributions.py -m 0.05 -v 2 -c 0.5
python makeDistributions.py -m 0.05 -v 2 -c 0.0
python makeDistributions.py -m 0.5 -v 2 -c 0.5
python makeDistributions.py -m 0.5 -v 2 -c 0.0



### plot the histograms, version 4
python makeDistributions.py -m 0.001 -v 4 -c 0.5
python makeDistributions.py -m 0.001 -v 4 -c 0.0
python makeDistributions.py -m 0.01 -v 4 -c 0.5
python makeDistributions.py -m 0.01 -v 4 -c 0.0
python makeDistributions.py -m 0.05 -v 4 -c 0.5
python makeDistributions.py -m 0.05 -v 4 -c 0.0
python makeDistributions.py -m 0.5 -v 4 -c 0.5
python makeDistributions.py -m 0.5 -v 4 -c 0.0
