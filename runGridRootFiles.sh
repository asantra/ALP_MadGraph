#! /bin/bash


lambdaValueList1="0.00000000316 0.00000000398 0.00000000501 0.00000000631 0.00000000794 0.00000001000 \
0.00000001259 0.00000001585"

lambdaValueList2="0.00000001995 0.00000002512 0.00000003162 \
0.00000003981 0.00000005012 0.00000006310 0.00000007943 0.00000010000"

lambdaValueList3="0.00000012589 0.00000015849 0.00000019953 0.00000025119 0.00000031623 \
0.00000039811 0.00000050119 0.00000063096"

lambdaValueList4="0.00000079433 0.00000100000 \ 
0.00000125893 0.00000158489 0.00000199526 0.00000251189 0.00000316228"

lambdaValueList5="0.00000398107 0.00000501187 0.00000630957 0.00000794328 0.00001000000 \
0.00001258925 0.00001584893"

lambdaValueList6="0.00001995262 0.00002511886 0.00003162278 \ 
0.00003981072 0.00005011872 0.00006309573 0.00007943282 0.00010000000"


lambdaValueList7="0.00012589254 0.00015848932 0.00019952623 0.00025118864 0.00031622777 \
0.00039810717 0.00050118723"

lambdaValueList8="0.00063095734 0.00079432823 0.00100000000 \
0.00125892541 0.00158489319 0.00199526231 0.00251188643 0.00316227766"

lambdaValueListSpecial="0.00002511886 0.00031622777 0.00039810717 0.00050118723 0.00063095734 0.00079432823 0.00100000000"

massValueList="0.01  0.0104713  0.0109648  0.0114815  0.0120226  0.0125893  \
0.0131826  0.0138038  0.0144544  0.0151356  0.0158489  0.0165959  \
0.017378  0.018197  0.0190546  0.0199526  0.020893  0.0218776  \
0.0229087  0.0239883  0.0251189  0.0263027  0.0275423  0.0288403  \
0.0301995  0.0316228  0.0331131  0.0346737  0.0363078  0.0380189  \
0.0398107  0.0416869  0.0436516  0.0457088  0.047863  0.0501187  \
0.0524807  0.0549541  0.057544  0.060256  0.0630957  0.0660693  \
0.0691831  0.0724436  0.0758578  0.0794328  0.0831764  0.0870964  \
0.0912011  0.0954993  0.1  0.104713  0.109648  0.114815  0.120226  \
0.125893  0.131826  0.138038  0.144544  0.151356  0.158489  0.165959  \
0.17378  0.18197  0.190546  0.199526  0.20893  0.218776  0.229087  \
0.239883  0.251189  0.263027  0.275423  0.288403  0.301995  0.316228  \
0.331131  0.346737  0.363078  0.380189  0.398107  0.416869  0.436516  \
0.457088  0.47863  0.501187  0.524807  0.549541  0.57544  0.60256  \
0.630957  0.660693  0.691831  0.724436  0.758578  0.794328  0.831764  \
0.870964  0.912011  0.954993  1."


massValueListSpecial="0.0380189"

#0.06294627058970834, 0.00012589254117941688, 13430.538545254938

versionValue="2"
cutWanted="0.0"
echo "symlinking libzstd"


DESTINATION="/storage/agrp/arkas/ALPFiles/findAcceptanceGridOutput"


echo "I am here:"$PWD
#echo "LD_LIBRARY_PATH is: "$LD_LIBRARY_PATH
echo "creating and copying the .so file"
root -l -b << EOF
.L drawGridRootFiles.C++
EOF
cp drawGridRootFiles_C.so $HOME
cp drawGridRootFiles_C_ACLiC_dict_rdict.pcm $HOME
cp drawGridRootFiles_C.d $HOME





if [[ $1 == "1" ]]
then
   counter=1
   lambdaList=${lambdaValueList1}
elif [[ $1 == "2" ]]
then
   counter=1000
   lambdaList=${lambdaValueList2}
elif [[ $1 == "3" ]]
then
   counter=2000
   lambdaList=${lambdaValueList3}
elif [[ $1 == "4" ]]
then
   counter=3000
   lambdaList=${lambdaValueList4}
elif [[ $1 == "5" ]]
then
   counter=4000
   lambdaList=${lambdaValueList5}
elif [[ $1 == "6" ]]
then
   counter=5000
   lambdaList=${lambdaValueList6}
elif [[ $1 == "7" ]]
then
   counter=6000
   lambdaList=${lambdaValueList7}
elif [[ $1 == "8" ]]
then
   counter=7000
   lambdaList=${lambdaValueList8}
else
   counter=8000
   lambdaList=${lambdaValueListSpecial}
fi


echo "The counter "$counter" for the lambdaList --> "${lambdaList}
echo "The grid output directory "${DESTINATION}"/Trial"${1}

runid=0
increase=0

### create the directory where the grid output will go
mkdir ${DESTINATION}"/Trial"${1}
### first working on Lambda value list 3, then 2 and then 4
for oneOverLambda in ${lambdaList}; do
    for massValue in ${massValueListSpecial}; do
        echo "oneOverLambda "${oneOverLambda}" and mass "${massValue}
        runid=$(( $increase + $counter ))
        mkdir -p ${DESTINATION}"/Trial"${1}"/run_"$runid
        PRESENTDIRECTORY=${PWD}
        echo "I am here: "$PRESENTDIRECTORY
        qsub -q N -N "run_${runid}" -l mem=8gb,vmem=10gb -o "${DESTINATION}/Trial"${1}"/run_"$runid -e "${DESTINATION}/Trial"${1}"/run_"$runid -v massValue=${massValue},versionValue=${versionValue},cutWanted=${cutWanted},oneOverLambda=${oneOverLambda} runRoot.sh 
        increase=$(( $increase + 1 ))
        break
        sleep 1s
    done
    break
done
