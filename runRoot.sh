#! /bin/bash


massValue=$1
versionValue=$2
cutWanted=$3
oneOverLambda=$4

HOME="/srv01/agrp/arkas"
ln -s /usr/local/anaconda3/pkgs/zstd-1.3.7-h0b5b093_0/lib/libzstd.so.1.3.7 $HOME/libzstd.so.1
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/anaconda3/lib:$PWD:$HOME
echo "Installing root"
cd $HOME/arkaRoot/root/bin
source thisroot.sh
cd


root -l << EOF
gSystem->Load("drawGridRootFiles_C.so");
//.L drawGridRootFiles.C++
std::cout << "---> adding the root files for mass " << ${massValue} << " and version " << ${versionValue} << " and cutValue wanted " << ${cutWanted}  << " and oneOverLambda " << ${oneOverLambda} << std::endl;
std::cout << "---> adding these files /storage/agrp/arkas/ALPFiles/outputRootFiles/LeadDumpFiles/mass${massValue}GeV/unweighted_events_mass${massValue}GeV_beam*GeV_v2.root" << std::endl;
TChain chain("LHEF")
chain.Add("/storage/agrp/arkas/ALPFiles/outputRootFiles/LeadDumpFiles/mass${massValue}GeV/unweighted_events_mass${massValue}GeV_beam*GeV_v2.root") // always want v2 root files, because the only thing that changes between versions are the lambda and tau
drawGridRootFiles t(&chain)
t.Loop("${massValue}", "${versionValue}","${cutWanted}", "${oneOverLambda}")
.q
EOF
