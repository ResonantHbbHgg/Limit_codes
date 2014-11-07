#!/bin/bash

#change windows for cuts
declare -a radion=("260" "270" "300" "350" "400" "450" "500"  )
declare -a winl0=("200" "250" "250" "300" "350" "400" "450" )
declare -a winu0=("280" "290" "350" "400" "450" "500" "550" )
#
declare -a winl1=("200" "250" "250" "300" "350" "400" "450" )
declare -a winu1=("280" "290" "350" "400" "450" "500" "550" )
# R2GGBBFitter_mgg_addhiggs.cc
# models_mtot_exp.rs

for (( i = 4 ; i < 5 ; i++ )); do # for each working point
  mkdir radlim_CSV_WP$i
#   for (( j = 0 ; j < ${#radion[@]} ; j++ )); do # for each mass
#  for (( j = 1 ; j < 7; j++ )); do # for each mass
#   for (( j = 0 ; j < 1 ; j++ )); do # for each mass
   for (( j = 1 ; j < 2 ; j++ )); do # for each mass
        # the signal file
        sed -i -r -e "s/WP[0-9]/WP$i/g" R2GGBBFitter_mgg_addhiggs.cc
        sed -i -r -e "s/m[0-9]+.root/m${radion[$j]}.root/g" R2GGBBFitter_mgg_addhiggs.cc 
        sed -i -r -e "s/Radion_m[0-9]+_8TeV_m[0-9]+.root/Radion_m${radion[$j]}_8TeV_m${radion[$j]}.root/g" R2GGBBFitter_mgg_addhiggs.cc
        # the ggH file
        #sed -i -r -e "s/WP[0-9]/WP$i/g" R2GGBBFitter_mgg_addhigg_ggH.cc
        #sed -i -r -e "s/m[0-9]+.root/m${radion[$j]}.root/g" R2GGBBFitter_mgg_addhigg_ggH.cc 
        # the ttH file
        #sed -i -r -e "s/WP[0-9]/WP$i/g" R2GGBBFitter_mgg_addhigg_ttH.cc
        #sed -i -r -e "s/m[0-9]+.root/m${radion[$j]}.root/g" R2GGBBFitter_mgg_addhigg_ttH.cc 
        # the signal file
        #sed -i -r -e "s/WP[0-9]/WP$i/g" R2GGBBFitter_mgg_addhigg_vbf.cc
        #sed -i -r -e "s/m[0-9]+.root/m${radion[$j]}.root/g" R2GGBBFitter_mgg_addhigg_vbf.cc 
        # the signal file
        #sed -i -r -e "s/WP[0-9]/WP$i/g" R2GGBBFitter_mgg_addhigg_vH.cc
        #sed -i -r -e "s/m[0-9]+.root/m${radion[$j]}.root/g" R2GGBBFitter_mgg_addhigg_vH.cc 
        # the legend
        sed -i -r -e "s/[0-9]+ GeV\"\);/${radion[$j]} GeV\"\);/g" R2GGBBFitter_mgg_addhiggs.cc
        # name of files
        #sed -i -r -e "s/mH[0-9]+/mH${radion[$j]}/g" R2GGBBFitter_mgg_addhiggs.cc
        #sed -i -r -e "s/m[0-9]+\_regression/${radion[$j]}\_regression/g" R2GGBBFitter_mgg_addhiggs.cc
        #sed -i -r -e "s/m[0-9]+\_/m${radion[$j]}\_/g" R2GGBBFitter_mgg_addhiggs.cc
        #sed -i -r -e "s/m[0-9]+\_/m${radion[$j]}\_/g" models_mtot_exp.rs
        # the window mtot > 550 || mtot < 450
        #sed -i -r -e "s/mtot > [0-9]+ \|\| mtot < [0-9]+/mtot > ${winu[$j]} \|\| mtot < ${winl[$j]}/g" R2GGBBFitter_mgg_addhiggs.cc
        echo WP$i MR ${radion[$j]}
        # the cuts
        # TString cut0 = "&& mtot > 455 && mtot < 550 "; //"&& 1>0";//
        # TString cut1 = "&& mtot > 455 && mtot < 550 "; // "&& 1>0";//
        # TString cutj0 = "&& mjj > 90 && mjj < 170 "; //"&& 1>0";//
        # TString cutj1 = "&& mjj > 100 && mjj < 160 "; // "&& 1>0";//
        # the signal parameters (median to fit)
#        sed -i -r -e "s/TString cut0 = \"\&\& mtot > [0-9]+ \&\& mtot < [0-9]+/TString cut0 = \"\&\& mtot > ${winl0[$j]} \&\& mtot < ${winu0[$j]}/g" R2GGBBFitter_mgg_addhiggs.cc
#        sed -i -r -e "s/TString cut1 = \"\&\& mtot > [0-9]+ \&\& mtot < [0-9]+/TString cut1 = \"\&\& mtot > ${winl1[$j]} \&\& mtot < ${winu1[$j]}/g" R2GGBBFitter_mgg_addhiggs.cc
        #
        mkdir radlim_CSV_WP$i/radlim${radion[$j]}_CSV/
        root -l -q runmgg.C >> radlim_CSV_WP$i/radlim${radion[$j]}_CSV/log_radlim${radion[$j]}.txt
        # overwrite the higgs WS
        #root -l -q runmgg_ggH.C >> radlim_CSV_WP$i/radlim${radion[$j]}_CSV/log_radlim${radion[$j]}_ggH.txt
        #root -l -q runmgg_ttH.C >> radlim_CSV_WP$i/radlim${radion[$j]}_CSV/log_radlim${radion[$j]}_ttH.txt
        #root -l -q runmgg_vbf.C >> radlim_CSV_WP$i/radlim${radion[$j]}_CSV/log_radlim${radion[$j]}_vbf.txt
        #root -l -q runmgg_vH.C >> radlim_CSV_WP$i/radlim${radion[$j]}_CSV/log_radlim${radion[$j]}_vH.txt
	mv workspaces/*.root radlim_CSV_WP$i/radlim${radion[$j]}_CSV
	mv datacards/*.txt radlim_CSV_WP$i/radlim${radion[$j]}_CSV

        #also colect the plots
        mv databkgoversig* radlim_CSV_WP$i/radlim${radion[$j]}_CSV/
        mv sigmodel* radlim_CSV_WP$i/radlim${radion[$j]}_CSV/
        mv higmodel* radlim_CSV_WP$i/radlim${radion[$j]}_CSV/
        mv remenber.txt radlim${radion[$j]}_CSV
        ## create limits root files for each mass
        cd radlim_CSV_WP$i/radlim${radion[$j]}_CSV/
	#
        #combine -M Asymptotic hgg.mH125.6_8TeVrep.txt >> higgsCombineTest.Asymptotic.mH125.6.mR${radion[$j]}.txt
        #mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}.root
	#echo did with reparametrization
	#
        combine -M Asymptotic hgg.mH125.0_8TeV.txt >> higgsCombineTest.Asymptotic.mH125.6.mR${radion[$j]}_higgs.txt
        mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}_higgs.root
	echo did with Higgs
	#
        #combine -M Asymptotic hgg.mH125.6_8TeVlnu.txt >> higgsCombineTest.Asymptotic.mH125.6.mR${radion[$j]}lnu.txt
        #mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}lnu.root
	#echo did with lnu
	#
        #combine -M Asymptotic hgg.mH125.6_8TeVrep.txt  -S 0 >> higgsCombineTest.Asymptotic.mH125.6.mR${radion[$j]}_nosyst.txt
        #mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}_nosyst.root
	#echo did no syst
	#
        #combine -M Asymptotic hgg.mH125.6_8TeVonecat.txt >> higgsCombineTest.Asymptotic.mH125.6.mR${radion[$j]}_onecat.txt
        #mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}_onecat.root
	#echo did one categ with higgs
	#
        combine -M Asymptotic hgg.mH125.0_8TeVonecatnohiggs.txt >> higgsCombineTest.Asymptotic.mH125.6.mR${radion[$j]}_onecatnohiggs.txt
        mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}_onecatnohiggs.root
	echo did one categ no higgs
	#
        combine -M Asymptotic hgg.mH125.0_8TeV.txt -S 0 >> higgsCombineTest.Asymptotic.mH125.6.mR${radion[$j]}_nosyst_higgs.txt
        mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR${radion[$j]}_nosyst_higgs.root
	echo did with no syst with higgs
        rm roostats*
        cd ../..
  done # mass
done # wp

