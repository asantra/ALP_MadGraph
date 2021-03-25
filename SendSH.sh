#! /bin/bash
# targetDirectory="/storage/agrp/arkas/ALPFiles/LocalAnalyzerDirectory"
targetDirectory="/storage/agrp/arkas/ALPFiles/GridAnalyzerDirectory"
sftp arkas@wipp-an1 << EOF
put drawGridRootFiles.* $targetDirectory/
put runGrid*.sh $targetDirectory/
put runRoot.sh $targetDirectory/
put outputPhotonNumbersUpdated.txt $targetDirectory/
put outputPhotonNumbers_photons_from_gbeam_JETI40_025fs_6500nm.txt $targetDirectory/
put outputPhotonNumbers_photons_from_gbeam_phase2_120fs_10000nm.txt $targetDirectory/
put outputPhotonNumbers_photons_from_gbeam_phase2_120fs_10000nm_tungsten.txt $targetDirectory/
EOF
