Limit_codes
===========

On CMSSW_6_1_1 you install the higgs combination tools according to slc5

https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideHiggsAnalysisCombinedLimit

As of 8 April 2015, it is necessary to make the following modification to the combination code.

edit `HiggsAnalysis/CombinedLimit/src/ToyMCSamplerOpt.cc` with:
   1. comment out line 156
   1. create a new line above 156 that reads: `nPA=nPA;`
	(scram treats warnings as errors, so this avoids the problem of nPA being unused.)

This modification forces the pseudo-asimov to be binned everywhere.
Likelihood evaluation will be correct, albeit slower.

This modification might be introduced as an option later in the slc6 version.

......................................................................
## Instructions

```bash
cd CMSSW_6_1_1/src/Limit_codes
cmsenv
source setup.sh
make

mkdir datacards/
mkdir workspaces/

# for 4 body fit, this runs a root script interactively
./fits_4body.sh

# for mgg and 2D fit, this runs a compiled program
./fits.sh

# for fitting all nonresonant samples either with Mgg or 2D, run
./fits_allnonres.sh
# the output is in limit_{dir_text}.txt
```
......................................................................
