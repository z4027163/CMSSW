#!/bin/bash


mkdir -p /eos/uscms/store/user/wangz/ntuple/jobdir
mkdir -p /eos/uscms/store/user/wangz/ntuple/llbb/2mu2b/loose/iso02 #/amcatnloFXFX

workdir=${PWD}
echo "Running Mono-HtoZZto4Leptons Analysis with executables RunHZZ4LeptonsAnalysis"


savedir=`echo /eos/uscms/store/user/wangz/ntuple/llbb/2mu2b/SingleMuon`

echo "Path is" $PATH
echo "Search Path is" $CMSSW_SEARCH_PATH

echo "Working dir is $workdir"
echo "Saving dir is $savedir"


#./RunReference_llbb_single NO Spring16 >& ${workdir}/HZZ4LeptonsAnalysis_log

mv ${workdir}/HZZ4LeptonsAnalysis_log ${savedir}/.
mv ${workdir}/output*.root    ${savedir}/.
mv ${workdir}/output_bnn.txt ${savedir}/.
mv ${workdir}/output_bnn.root ${savedir}/.
mv ${workdir}/output_txt.txt ${savedir}/.



