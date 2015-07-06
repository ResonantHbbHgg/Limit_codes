#!/bin/bash

doBlinding=0

declare -a radion=("350" "400" "450" "500" "550" "600" "650" "700" "800" "900" "1000" "1100" "1200" "1300" "1400")
declare -a winl=("300" "350" "400" "450" "500" "550" "600" "650" "750" "800" "900" "1000" "1100" "1200" "1300")
declare -a winu=("400" "450" "500" "550" "600" "650" "700" "750" "900" "1000" "1100" "1200" "1300" "1400" "1500")

mkdir radlim_Mggjj

for (( j = 1 ; j < 12 ; j++ )); do # for each mass 

    sed -i -r -e "s/m[0-9]+.rs/m${radion[$j]}.rs/g" R2GGBBFitter_mtot_range.cc
    sed -i -r -e "s/m[0-9]+_8TeV_m[0-9]+.root/m${radion[$j]}_8TeV_m${radion[$j]}.root/g" R2GGBBFitter_mtot_range.cc

    sed -i -r -e "s/mtot_sig_m0\[[0-9]+, [0-9]+, [0-9]+/mtot_sig_m0\[${radion[$j]}, ${winl[$j]}, ${winu[$j]}/g" models_mtot_range.rs
    sed -i -r -e "s/mtot_sig_m0_cat0\[[0-9]+, [0-9]+, [0-9]+/mtot_sig_m0_cat0\[${radion[$j]}, ${winl[$j]}, ${winu[$j]}/g" models_mtot_range.rs
    sed -i -r -e "s/mtot_sig_m0_cat1\[[0-9]+, [0-9]+, [0-9]+/mtot_sig_m0_cat1\[${radion[$j]}, ${winl[$j]}, ${winu[$j]}/g" models_mtot_range.rs

    # the legend
    sed -i -r -e "s/[0-9]+ GeV\"\);/${radion[$j]} GeV\"\);/g" R2GGBBFitter_mtot_range.cc
    sed -i -r -e "s/mtot > [0-9]+ \|\| mtot < [0-9]+/mtot > ${winu[$j]} \|\| mtot < ${winl[$j]}/g" R2GGBBFitter_mtot_range.cc

    echo MR ${radion[$j]}
    # the mass to fit around
    sed -i -r -e "s/runfits\([0-9]+/runfits\(${radion[$j]}/g" runfits.C

    outputdir="radlim_Mggjj/radlim${radion[$j]}"
    mkdir $outputdir
    root -l -b -q runfits.C >> ${outputdir}/log_radlim${radion[$j]}.txt
    mv workspaces/hgghbb.* ${outputdir}
    mv datacards/* ${outputdir}
    #also colect the plots
    mv databkgoversig*  ${outputdir}
    mv sigmodel*  ${outputdir}

    # create limits root files for each mass
    cd ${outputdir}/
    cp ../../models_mtot_range_m${radion[$j]}.rs .

    if [ ${doBlinding} == 1 ]    
    then
	combine -M Asymptotic --run blind hgghbb.mH${radion[$j]}.0_8TeVonecat.txt >> higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}onecat.txt
	mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mR${radion[$j]}_onecat.root
	echo did with rep 2btag only
	       
	combine -M Asymptotic --run blind hgghbb.mH${radion[$j]}.0_8TeVrep.txt -S 0 >> higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}_nosyst.txt
	mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mR${radion[$j]}_nosyst.root
	echo with rep no syst	

	combine -M Asymptotic --run blind hgghbb.mH${radion[$j]}.0_8TeVrep.txt >> higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}_reparametrization.txt
	mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}.root
	echo did with rep 
	cd ../..
    else
	combine -M Asymptotic hgghbb.mH${radion[$j]}.0_8TeVonecat.txt >> higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}onecat.txt
	mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mR${radion[$j]}_onecat.root
	echo did with rep 2btag only
	
	combine -M Asymptotic hgghbb.mH${radion[$j]}.0_8TeVrep.txt -S 0 >> higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}_nosyst.txt
	mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mR${radion[$j]}_nosyst.root
	echo with rep no syst	

	combine -M Asymptotic hgghbb.mH${radion[$j]}.0_8TeVrep.txt >> higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}_reparametrization.txt
	mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}.root
	echo did with rep 
	cd ../..
    fi   
done # mass

#  root -l -q Brazilianflag.cc
#sed -i -r -e "s/mR300.Xcut/Mcut_cutbased/g" Brazilianflag.cc #adapt this name to the cut
#root -l Brazilianflag.cc

