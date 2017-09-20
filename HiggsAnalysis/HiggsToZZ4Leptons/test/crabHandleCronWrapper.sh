#!/bin/bash

#Takes care of necessary environment setup and such for cron to monitor crab tasks
#To setup: first argument should be SETUP, followed by CMSSW directory (wherever cmsenv should be run)
#Rest should just be directory of crab directories (will be a timestamp)a

#Example: bash crabHandleCronWrapper.sh SETUP ~/nobackup/CMSSW_X_Y_Z/src CRAB_DIR [crab_sadhjashdjk crab_asdhjkashdkjashd]

set -e 

if [ "$1" == "SETUP" ]; then
    shift
    if [[ "$@" =~ "%" ]]; then
        echo "No % signs allowed. crontab will lose its mind over those."
        exit 1
    fi
    if [ -d "$1" ]; then
        cat <(crontab -l 2>/dev/null) <(echo $(($(date +%M)+1)) '*' '*' '*' '*' 'bash' "`readlink -f "$BASH_SOURCE"`" `printf '%q ' "$(echo $(date +%s%N) "$@" | sha256sum | awk '{print$1}')" "$@"`) | crontab -
    else 
        exit 1
    fi
    exit 0
fi

removeCrontab() {
    crontab -l | grep -Fv "$SHA_SUM" | crontab -
}

SHA_SUM="$1"
shift
CMSSW_DIR="$1"
shift
#crontab -l
#removeCrontab "$@"
#crontab -l
#exit
. /cvmfs/cms.cern.ch/cmsset_default.sh
cd "$CMSSW_DIR"
eval `scramv1 runtime -sh`
voms-proxy-info #Check for proxy. set -e will kill the script if it is missing
. /cvmfs/cms.cern.ch/crab3/crab.sh
python crabHandle.py "$@" && removeCrontab "$@"
