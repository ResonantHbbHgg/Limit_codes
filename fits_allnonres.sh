#!/bin/bash

doBlinding=0

version=44
limitdirs=("fitToMgg_nonresSearch_withKinFit" "fitTo2D_nonresSearch_withKinFit")
do2DLimits=("0" "1")

#If you only want to run on a subset of directories, edit this array with the appropriate indices.
#It is recommended to only run on one of "0" or "1" at a time. There are lots of nonres samples.
runLimits=("1")

for i in `echo ${runLimits[@]}`; do

    limitOutputFile="limits_${limitdirs[$i]}.txt"
    if [ -f $limitOutputFile ]; then
	mv $limitOutputFile ${limitOutputFile}_old
    fi
    echo "c2 yt lambda_hhh -2_sigma -1_sigma median +1_sigma +2_sigma observed" > $limitOutputFile

    if [ ${do2DLimits[$i]} == "0" ]; then
	fitterScript=R2GGBBFitter_mgg_addhiggs.exe
    else
	fitterScript=R2GGBBFitter_2D_addhiggs.exe
    fi

    sampleList=(`cat nonres_samplelist.txt`)

    for isample in `echo ${sampleList[@]}`; do

	echo "Running limits for sample $isample on ${limitdirs[$i]}"

        #if [ "$isample" != "Lam_1d0_Yt_1d0_c2_0d0" ]; then
        #     continue
        #fi

	outputdir="radlim_${limitdirs[$i]}/${isample}"
	mkdir -p $outputdir

	if [ ${do2DLimits[$i]} == "0" ]; then
	    fitter="R2GGBBFitter_mgg_addhiggs.exe"
	else
	    fitter="R2GGBBFitter_2D_addhiggs.exe"
	fi

	./$fitter -v $version -n 4 --sigMass 0 --analysisType ${limitdirs[$i]} --nonresFile $isample --useSigTheoryUnc 0 >& ${outputdir}/log_radlim0.txt

	mv workspaces/*.root $outputdir
	mv datacards/*.txt $outputdir

       #also colect the plots
	mv databkgoversig* $outputdir
	mv sigmodel* $outputdir
	mv higmodel* $outputdir

        ## create limits root files for each mass
	cd $outputdir
	
	outputFile="higgsCombineTest.Asymptotic.mH125.0.mR0_higgs.txt"

	if [ ${doBlinding} == 1 ]    
	then
	    combine -M Asymptotic --run blind hgg.mH125.0_8TeV.txt >> $outputFile
	    mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR0_higgs.root
	    echo did with Higgs
	        
	    rm roostats*
	    rm log_radlim0.txt #this is a big file
	    cd ../..

	    ./outputNonresLimits.py $outputdir/$outputFile $isample $limitOutputFile

	else
	    combine -M Asymptotic hgg.mH125.0_8TeV.txt >> $outputFile
	    mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR0_higgs.root
	    echo did with Higgs
	        
	    rm roostats*
	    rm log_radlim0.txt
	    cd ../..

	    ./outputNonresLimits.py $outputdir/$outputFile $isample $limitOutputFile

	fi
    done # sample
done #dirs


