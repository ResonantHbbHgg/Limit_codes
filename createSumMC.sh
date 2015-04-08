#!/bin/bash

version=v38
basedir=/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees
#dir=$basedir/$version/${version}_fitToFTR14001_nonresSearch_withKinFit
#dir=$basedir/$version/${version}_fitTo2D_nonresSearch_withKinFit
dir=$basedir/$version/${version}_fitToMgg_nonresSearch_withKinFit


hadd -f sumMC_m0.root \
$dir/DYJetsToLL_m0.root \
$dir/LNuGG_FSR_8TeV_m0.root \
$dir/LNuGG_ISR_8TeV_m0.root \
$dir/TTGJets_8TeV_m0.root \
$dir/ZGToLLG_8TeV_m0.root \
$dir/bbh_m125_8TeV_m0.root \
$dir/diphojet_sherpa_8TeV_m0.root \
$dir/ggh_m125_powheg_8TeV_m0.root \
$dir/gjet_20_8TeV_pf_m0.root \
$dir/gjet_40_8TeV_pf_m0.root \
$dir/qcd_30_8TeV_ff_m0.root \
$dir/qcd_30_8TeV_pf_m0.root \
$dir/qcd_40_8TeV_ff_m0.root \
$dir/qcd_40_8TeV_pf_m0.root \
$dir/tGG_8TeV_m0.root \
$dir/ttGG_8TeV_m0.root \
$dir/tth_m125_8TeV_m0.root \
$dir/vbf_m125_8TeV_m0.root \
$dir/wzh_m125_8TeV_wh_m0.root \
$dir/wzh_m125_8TeV_zh_m0.root
