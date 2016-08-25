HiggsAnalysis-ZZMatrixElement
=============================

[Wiki page for the MEM packages](https://twiki.cern.ch/twiki/bin/viewauth/CMS/HZZ4lME)


Checking out on top of other CMSSW packages
-------------------------------------------

Assuming that you have set up a CMSSW area with packages from the cms-sw/cmssw git repository (`git cms-init`),  you have to set up this repository so that its .git folder does not interfere with the cms-sw/cmssw folder:

```
git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement
```

At this point, git commands apply to this repository if issued inside the ZZMatrixElement folder; or to the cmssw repository if issued outside (and in this case ZZMatrixElement is reported as an untracked folder by git status).


The package depends on:
```
## >>> CombinedLimit <<<
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git checkout -b from-V02-06-00 HiggsAnalysis-CombinedLimit-V02-06-00
cd -
## >>> Higgs_CS_and_Width <<<
git clone https://github.com/msnowball/HCSaW Higgs/Higgs_CS_and_Width
cd Higgs/Higgs_CS_and_Width
git filter-branch --subdirectory-filter Higgs_CS_and_Width
cd -
```

Mixing with packages from cvs
-------------------------------------------
Just set up the git part first, then checkout the cvs packages as usual. 
They will appear as untracked folders in git status.

You can move the ZZAnalysis folder into an existing area if needed, it will not break git providing that ZZAnalysis/.git/ is moved as well.


Committing and pushing
-------------------------------------------

```
git pull
[edit files]
git add [files to be added]
git commit -m ["commit message"] [files to be added]
git push origin master
```
