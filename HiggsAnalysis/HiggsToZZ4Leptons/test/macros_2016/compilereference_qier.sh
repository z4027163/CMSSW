#!/bin/bash

exedir=${CMSSW_BASE}/src/HiggsAnalysis/HiggsToZZ4Leptons/test/macros_2016

melalibdir=${CMSSW_BASE}/lib/slc6_amd64_gcc530
melaincdir=${CMSSW_BASE}/src

export LD_LIBRARY_PATH=${melalibdir}:$LD_LIBRARY_PATH

    g++ -I $ROOTSYS/include -I $ROOFITSYS/include -I ${melaincdir} compilereference_single.C HZZ4LeptonsAnalysis.C `root-config --glibs` `root-config --libs` `root-config --cflags` -L $ROOFITSYS/lib  -lRooFit -lRooFitCore -L ${melalibdir} -o RunReferenceqier
#    g++ -I $ROOTSYS/include -I $ROOFITSYS/include -I ${melaincdir} compilereference_bkg.C HZZ4LeptonsAnalysis.C `root-config --glibs` `root-config --libs` `root-config --cflags` -L $ROOFITSYS/lib  -lRooFit -lRooFitCore -L ${melalibdir} -o RunReferenceqier_bkg
#    g++ -I $ROOTSYS/include -I $ROOFITSYS/include -I ${melaincdir} compilereference_signal.C HZZ4LeptonsAnalysis.C `root-config --glibs` `root-config --libs` `root-config --cflags` -L $ROOFITSYS/lib  -lRooFit -lRooFitCore -L ${melalibdir} -o RunReferenceqier_signal
#    g++ -I $ROOTSYS/include -I $ROOFITSYS/include -I ${melaincdir} compilereference_data.C HZZ4LeptonsAnalysis.C `root-config --glibs` `root-config --libs` `root-config --cflags` -L $ROOFITSYS/lib  -lRooFit -lRooFitCore -L ${melalibdir}  -o RunReferenceqier_data
