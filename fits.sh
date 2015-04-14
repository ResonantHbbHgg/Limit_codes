#!/bin/bash

doBlinding=1

version=44
limitdirs=("fitToMgg_nonresSearch_withKinFit" "fitToMgg_resSearch_withKinFit" "fitToMgg_resSearch_withRegKinFit" "fitTo2D_nonresSearch_withKinFit" "fitTo2D_resSearch_withRegKinFit" "fitTo2D_resSearch_withKinFit" "fitToFTR14001_nonresSearch_withKinFit")
doResLimits=("0" "1" "1" "0" "1" "1" "0")
do2DLimits=("0" "0" "0" "1" "1" "1" "1")

#If you only want to run on a subset of directories, edit this array with the appropriate indices.
runLimits=("3" "4")

for i in `echo ${runLimits[@]}`; do

    if [ ${doResLimits[$i]} == "0" ]; then
	masses=("0")
    else
	masses=("260" "270" "300" "350" "400") #the limit trees exist to do 450 and 500 too
    fi

    if [ ${do2DLimits[$i]} == "0" ]; then
	fitter="R2GGBBFitter_mgg_addhiggs.exe"
    else
	fitter="R2GGBBFitter_2D_addhiggs.exe"
    fi

    for imass in `echo ${masses[@]}`; do

	echo "Running limits for mass $imass on ${limitdirs[$i]}"

	outputdir="radlim_${limitdirs[$i]}/radlim${imass}"
	mkdir -p $outputdir

	ncat=2
	useSigTheoryUnc=0
	if [ "$imass" == "0" ]; then
	    ncat=4
	    useSigTheoryUnc=1
	fi

	./$fitter -v $version -n $ncat --sigMass $imass --analysisType ${limitdirs[$i]} --useSigTheoryUnc ${useSigTheoryUnc} >& ${outputdir}/log_radlim${imass}.txt

	mv workspaces/*.root $outputdir
	mv datacards/*.txt $outputdir

       #also colect the plots
	mv databkgoversig* $outputdir
	mv sigmodel* $outputdir
	mv higmodel* $outputdir

        ## create limits root files for each mass
	cd $outputdir
	
	if [ ${doBlinding} == 1 ]    
	then
	    combine -M Asymptotic --run blind hgg.mH125.0_8TeV.txt >> higgsCombineTest.Asymptotic.mH125.0.mR${imass}_higgs.txt
	    mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR${imass}_higgs.root
	    echo did with Higgs
	        
	    combine -M Asymptotic --run blind hgg.mH125.0_8TeVonecatnohiggs.txt >> higgsCombineTest.Asymptotic.mH125.0.mR${imass}_onecatnohiggs.txt
	    mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR${imass}_onecatnohiggs.root
	    echo did one categ no higgs

	    combine -M Asymptotic --run blind hgg.mH125.0_8TeV.txt -S 0 >> higgsCombineTest.Asymptotic.mH125.0.mR${imass}_nosyst_higgs.txt
	    mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR${imass}_nosyst_higgs.root
	    echo did with no syst with higgs

	    rm roostats*
	    cd ../..

	else
	    combine -M Asymptotic hgg.mH125.0_8TeV.txt >> higgsCombineTest.Asymptotic.mH125.0.mR${imass}_higgs.txt
	    mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR${imass}_higgs.root
	    echo did with Higgs
	        
	    combine -M Asymptotic hgg.mH125.0_8TeVonecatnohiggs.txt >> higgsCombineTest.Asymptotic.mH125.0.mR${imass}_onecatnohiggs.txt
	    mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR${imass}_onecatnohiggs.root
	    echo did one categ no higgs

	    combine -M Asymptotic hgg.mH125.0_8TeV.txt -S 0 >> higgsCombineTest.Asymptotic.mH125.0.mR${imass}_nosyst_higgs.txt
	    mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR${imass}_nosyst_higgs.root
	    echo did with no syst with higgs

	    rm roostats*
	    cd ../..

	fi
    done # mass
done #dirs


