#!/bin/bash


mkdir -p /eos/uscms/store/user/wangz/ntuple/jobdir
mkdir -p /eos/uscms/store/user/wangz/ntuple/4l2b #llbb/amcatnloFXFX

workdir=${PWD}
echo "Running Mono-HtoZZto4Leptons Analysis with executables RunHZZ4LeptonsAnalysis"


savedir=`echo /eos/uscms/store/user/wangz/ntuple/4l2b/` #llbb/amcatnloFXFX

echo "Path is" $PATH
echo "Search Path is" $CMSSW_SEARCH_PATH

echo "Working dir is $workdir"
echo "Saving dir is $savedir"


#./RunReference_4l2b 2016 NO >& ${workdir}/4l2b_log

mv ${workdir}/4l2b_log ${savedir}/.
mv ${workdir}/output*.root    ${savedir}/.
mv ${workdir}/output_bnn.txt ${savedir}/.
mv ${workdir}/output_bnn.root ${savedir}/.
mv ${workdir}/output_txt.txt ${savedir}/.



