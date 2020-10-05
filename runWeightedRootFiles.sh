#! /bin/bash

massValue="$1"
versionValue="$2"
cutWanted="$3"
root -l -b << EOF
.L drawWeightedRootFiles.C++
std::cout << "---> adding the root files for mass " << ${massValue} << " and version " << ${versionValue} << " and cutValue wanted " << ${cutWanted} << std::endl;
TChain chain("LHEF")
chain.Add("/Users/arkasantra/MadGraphWorks/outputRootFiles/unweighted_events_mass${massValue}GeV_beam*GeV_v${versionValue}.root")
drawWeightedRootFiles t(&chain)
t.Loop("${massValue}", "${versionValue}","${cutWanted}")
.q
EOF
