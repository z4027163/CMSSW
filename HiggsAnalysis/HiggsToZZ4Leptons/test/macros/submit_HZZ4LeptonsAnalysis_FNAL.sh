#!/bin/bash


mkdir -p $_CONDOR_SCRATCH_DIR/histodir
mkdir -p /eos/uscms/store/user/ndefilip/jobdir
mkdir -p /eos/uscms/store/user/ndefilip/histodir

workdir=${PWD}
echo "Running HtoZZto4Leptons Analysis with executables RunHZZ4LeptonsAnalysis"
source /uscmst1/prod/sw/cms/bashrc prod
exedir=`echo /uscms/home/zwang4/nobackup/WORKSPCACE/tem/CMSSW_7_6_3/src/HiggsAnalysis/HiggsToZZ4Leptons/test/macros`
export LD_LIBRARY_PATH=${melalibdir}:$LD_LIBRARY_PATH
cd ${exedir}
eval `scramv1 runtime -sh`

melalibdir=${CMSSW_BASE}/lib/$SCRAM_ARCH/
export LD_LIBRARY_PATH=${melalibdir}:$LD_LIBRARY_PATH

if [ -d "$_CONDOR_SCRATCH_DIR/" ]; then
    workdir=`echo $_CONDOR_SCRATCH_DIR/`;
    cd ${workdir};
fi



# savedir=`echo ${exedir}/log`
savedir=`echo /eos/uscms/store/user/wangz/histodir`

echo "Working dir is $workdir"
echo "Executable dir is $exedir"
echo "Saving dir is $savedir"

${exedir}/RunReferenceAnalysis ${exedir}/sig_input_h150.txt 1 ${exedir}/bkg_input.txt 1 ${exedir}/data_input.txt 1 site year mc >& ${workdir}/HZZ4LeptonsAnalysis_log
mv ${workdir}/HZZ4LeptonsAnalysis_log /eos/uscms/store/user/ndefilip/jobdir/HZZ4LeptonsAnalysis_log
# rm -f  ${workdir}/HZZ4LeptonsAnalysis_log /eos/uscms/store/user/ndefilip/jobdir/HZZ4LeptonsAnalysis_log

mv ${workdir}/output.root    ${savedir}/.
mv ${workdir}/output_bnn.txt ${savedir}/.
mv ${workdir}/output_bnn.root ${savedir}/.
mv ${workdir}/output_txt.txt ${savedir}/.
mv ${workdir}/output_txt_vbf.txt ${savedir}/.


