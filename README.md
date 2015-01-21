Limit_codes
===========

On CMSSW_6_1_1 you install the higgs combination tools

https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideHiggsAnalysisCombinedLimit

......................................................................

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
......................................................................
