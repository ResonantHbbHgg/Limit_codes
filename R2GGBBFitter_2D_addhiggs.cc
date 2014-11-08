/** \macro H2GGFitter.cc
*
* $Id$
*
* Software developed for the CMS Detector at LHC
*
*
* Template Serguei Ganjour - CEA/IRFU/SPP, Saclay
*
*
* Macro is implementing the unbinned maximum-likelihood model for
* the Higgs to gamma gamma analysis. PDF model and RooDataSets
* are stored in the workspace which is feeded to HiggsAnalysis/CombinedLimit tools:
*
*/
// this one is for mgg fit
using namespace RooFit;
using namespace RooStats ;
const Int_t NCAT = 2;
float integral[NCAT];
bool addHiggs=true;
void AddSigData(RooWorkspace*, Float_t);
void AddHigData(RooWorkspace*, Float_t,int);
void AddBkgData(RooWorkspace*);
void SigModelFit(RooWorkspace*, Float_t);
void MakePlots(RooWorkspace*, Float_t);
void MakePlotsHiggs(RooWorkspace* w, Float_t Mass);
void MakeSigWS(RooWorkspace* w, const char* filename);
void MakeHigWS(RooWorkspace* w, const char* filename,int);
void MakeBkgWS(RooWorkspace* w, const char* filename);//, 
         // const char* filenameh0, const char* filenameh1, const char* filenameh2, const char* filenameh4);
void MakeDataCard(RooWorkspace* w, const char* filename, const char* filename1,
                  const char* filename2, const char*, const char*, const char*);
void MakeDataCardREP(RooWorkspace* w, const char* filename, const char* filename1);
void MakeDataCardonecat(RooWorkspace* w, const char* filename, const char* filename1, 
                         const char* filename2, const char*, const char*, const char*);
void MakeDataCardLNU(RooWorkspace* w, const char* filename, const char* filename1);
void MakeDataCardonecatnohiggs(RooWorkspace* w, const char* filename, const char* filename1, const char* filename2, const char*, const char*, const char*);
// MakeDataCardonecat(w, fileBaseName, fileBkgName);
// MakeDataCardREP(w, fileBaseName, fileBkgName);
//void MakeDataCardnohiggs(RooWorkspace* w, const char* filename, const char* filename1, const char* filename2);
//void MakeDataCardonecat(RooWorkspace* w, const char* filename, const char* filename1, const char* filename2);
void SetParamNames(RooWorkspace*);
void SetConstantParams(const RooArgSet* params);

RooFitResult* fitresult[NCAT]; // container for the fit results
RooFitResult* BkgModelFitBernstein(RooWorkspace*, Bool_t);

RooArgSet* defineVariables()
{
  RooRealVar* mgg = new RooRealVar("mgg","M(#gamma#gamma)",100,180,"GeV");
  //RooRealVar* mtot = new RooRealVar("mtot","M(#gamma#gammajj)",200,1600,"GeV");
  RooRealVar* mjj = new RooRealVar("mjj","M(jj)",60,180,"GeV");
  RooRealVar* evWeight = new RooRealVar("evWeight","HqT x PUwei",0,100,"");
  RooCategory* cut_based_ct = new RooCategory("cut_based_ct","event category 4") ;
  //
  cut_based_ct->defineType("cat4_0",0);
  cut_based_ct->defineType("cat4_1",1);
  //
  RooArgSet* ntplVars = new RooArgSet(*mgg, *mjj, *cut_based_ct, *evWeight);
  ntplVars->add(*mgg);
  //ntplVars->add(*mtot);
  ntplVars->add(*mjj);
  ntplVars->add(*cut_based_ct);
  return ntplVars;
}

void runfits(const Float_t mass=120, Int_t mode=1, Bool_t dobands = false)
{
  style();
  TString fileBaseName(TString::Format("hgg.mH%.1f_8TeV", mass));
  TString fileHiggsNameggh(TString::Format("hgg.hig.mH%.1f_8TeV.ggh", mass));
  TString fileHiggsNametth(TString::Format("hgg.hig.mH%.1f_8TeV.tth", mass));
  TString fileHiggsNamevbf(TString::Format("hgg.hig.mH%.1f_8TeV.vbf", mass));
  TString fileHiggsNamevh(TString::Format("hgg.hig.mH%.1f_8TeV.vh", mass));
  TString fileBkgName(TString::Format("hgg.inputbkg_8TeV", mass));
  TString card_name("models_test.rs"); // put the model parameters here!
  HLFactory hlf("HLFactory", card_name, false);

  RooFitResult* fitresults;
  bool cutbased=true;
  // the minitree to be addeed
  //
  TString hhiggsggh = "/afs/cern.ch/work/h/hebda/HggHbb/CMSSW_6_1_1/src/BJetRegression/v30_withRunLumiEvent_fitToMgg_noKinFit/ggh_m125_powheg_8TeV_m300.root";
  TString hhiggstth = "/afs/cern.ch/work/h/hebda/HggHbb/CMSSW_6_1_1/src/BJetRegression/v30_withRunLumiEvent_fitToMgg_noKinFit/tth_m125_8TeV_m300.root";
  TString hhiggsvbf = "/afs/cern.ch/work/h/hebda/HggHbb/CMSSW_6_1_1/src/BJetRegression/v30_withRunLumiEvent_fitToMgg_noKinFit/vbf_m125_8TeV_m300.root";
  TString hhiggsvh = "/afs/cern.ch/work/h/hebda/HggHbb/CMSSW_6_1_1/src/BJetRegression/v30_withRunLumiEvent_fitToMgg_noKinFit/wzh_m125_8TeV_zh_m300.root";///afs/cern.ch/work/a/acarvalh/CMSSW_6_1_1/src/code/Limit_codes/MiniTrees/v28/wzh_m125_8TeV_m300.root";
  //

  //ftr analysis
  //TString ssignal = "/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v36/v36_fitToFTR14001_nonresSearch_noKinFit/ggHH_8TeV_m0.root";
  //TString ddata = "/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v36/v36_fitToFTR14001_nonresSearch_noKinFit/DataCS_m0.root";
  //TString ddata = "/afs/cern.ch/work/h/hebda/HggHbb/CMSSW_6_1_1/src/fitter_2d/v36_FTR/sumMC_m0.root";//sum of all MC bkgs, res and non-res

  //baseline analysis v20
  TString ssignal = "/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v36/v36_fitTo2D_nonresSearch_withKinFit/ggHH_8TeV_m0.root";
  //TString ddata = "/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v36/v36_fitTo2D_nonresSearch_withKinFit/DataCS_m0.root";
  TString ddata = "/afs/cern.ch/work/h/hebda/HggHbb/CMSSW_6_1_1/src/fitter_2d/v36_BaseSM/sumMC_m0.root";//sum of all MC bkgs, res and non-res

  //Radion m300
  //TString ssignal = "/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v36/v36_fitTo2D_resSearch_withRegKinFit/Radion_m300_8TeV_m300.root";
  //TString ddata = "/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v36/v36_fitTo2D_resSearch_withRegKinFit/DataCS_m300.root";
  //TString ddata = "/afs/cern.ch/work/h/hebda/HggHbb/CMSSW_6_1_1/src/fitter_2d/v36_m300/sumMC_m300.root";//sum of all MC bkgs, res and non-res


  //TString ssignal = "/afs/cern.ch/work/h/hebda/HggHbb/CMSSW_6_1_1/src/Selection/v37_fitToMgg_noKinFit_noRegression_2D/Radion_m300_8TeV_m300.root ";
  //TString ddata = "/afs/cern.ch/work/h/hebda/HggHbb/CMSSW_6_1_1/src/Selection/v37_fitToMgg_noKinFit_noRegression_2D/Data_m300.root";
  //TString ssignal = "/afs/cern.ch/work/h/hebda/HggHbb/CMSSW_6_1_1/src/Selection/v37_fitToMgg_noKinFit_noRegression_2D/MSSM_m260_8TeV_m300.root";
  //TString ddata = "/afs/cern.ch/work/h/hebda/HggHbb/CMSSW_6_1_1/src/Selection/v37_fitToMgg_noKinFit_noRegression_2D/Data_m300.root";
  //TString ssignal = "/afs/cern.ch/work/h/hebda/HggHbb/CMSSW_6_1_1/src/Selection/v37_fitToMgg_noKinFit_noRegression_2D/ggHH_8TeV_m300.root";
  //TString ddata = "/afs/cern.ch/work/h/hebda/HggHbb/CMSSW_6_1_1/src/Selection/v37_fitToMgg_noKinFit_noRegression_2D/Data_m300.root";
  //
  // TString hhiggs = "MiniTrees/OlivierOc13/v16_base_mgg_0_massCutVersion0/02013-11-05-Radion_m300_8TeV_nm_m300.root";
  // TString ssignal = "MiniTrees/OlivierOc13/v16_base_mgg_0_massCutVersion0/02013-11-05-Radion_m300_8TeV_nm_m300.root";
  // TString ddata = "MiniTrees/OlivierOc13/v16_base_mgg_0_massCutVersion0/02013-11-05-Data_m300.root";
   //
  // TString hhiggs = "MiniTrees/OlivierOc13/v15_regkin_mgg_0_massCutVersion0/02013-10-30-Radion_m300_8TeV_nm_m300.root";
  // TString ssignal = "MiniTrees/OlivierOc13/v15_regkin_mgg_0_massCutVersion0/02013-10-30-Radion_m300_8TeV_nm_m300.root";
  // TString ddata = "MiniTrees/OlivierOc13/v15_regkin_mgg_0_massCutVersion0/02013-10-30-Data_m300.root";
  //
  cout<<"Signal: "<<ssignal<<endl;
  cout<<"Data: "<<ddata<<endl;
  //
  RooWorkspace* w = hlf.GetWs();
  AddSigData(w, mass,ssignal);
  cout<<"SIGNAL ADDED"<<endl;
  SigModelFit(w, mass); // constructing signal pdf
  MakeSigWS(w, fileBaseName);
  MakePlots(w, mass);
  cout<<" did signal WS's"<<endl;
  //

  cout<<"Higgs: "<<hhiggsggh<<endl;
  AddHigData(w, mass,hhiggsggh,0);
  HigModelFit(w, mass,0); // constructing higgs pdf
  MakeHigWS(w, fileHiggsNameggh,0);
  //
  cout<<"Higgs: "<<hhiggstth<<endl;
  AddHigData(w, mass,hhiggstth,1);
  HigModelFit(w, mass,1); // constructing higgs pdf
  MakeHigWS(w, fileHiggsNametth,1);
  //
  cout<<"Higgs: "<<hhiggsvbf<<endl;
  AddHigData(w, mass,hhiggsvbf,2);
  HigModelFit(w, mass,2); // constructing higgs pdf
  MakeHigWS(w, fileHiggsNamevbf,2);
  //
  cout<<"Higgs: "<<hhiggsvh<<endl;
  AddHigData(w, mass,hhiggsvh,3);
  HigModelFit(w, mass,3); // constructing higgs pdf
  MakeHigWS(w, fileHiggsNamevh,3);
  cout<<"HIGGS ADDED"<<endl;
  //MakePlotsHiggs(w, mass);
  //

  AddBkgData(w,ddata);
  w->Print("v");
  cout<<"BKG ADDED"<<endl;
  bool dobands=true;
  fitresults = BkgModelFitBernstein(w, dobands); // this is berestein 3
  MakeBkgWS(w, fileBkgName);
  // construct the models to fit
  //
  MakeDataCardonecat(w, fileBaseName, fileBkgName, fileHiggsNameggh, fileHiggsNametth, fileHiggsNamevbf, fileHiggsNamevh);
  MakeDataCardREP(w, fileBaseName, fileBkgName);
  MakeDataCardLNU(w, fileBaseName, fileBkgName);
  MakeDataCardonecatnohiggs(w, fileBaseName, fileBkgName, fileHiggsNameggh, fileHiggsNametth, fileHiggsNamevbf, fileHiggsNamevh);
  MakeDataCard(w, fileBaseName, fileBkgName, fileHiggsNameggh, fileHiggsNametth, fileHiggsNamevbf, fileHiggsNamevh);
  // MakeDataCardonecat(w, fileBaseName, fileBkgName, fileHiggsName);//MakeDataCardnohiggs 
  cout<< "here"<<endl;

  return;
} // close runfits
////////////////////////////////////////////////////////////////////
// we add the data to the workspace in categories
void AddSigData(RooWorkspace* w, Float_t mass, TString signalfile) {
  cout << "================= Add Signal==============================" << endl;
  const Int_t ncat = NCAT;
  Float_t MASS(mass);
  // Luminosity:
  Float_t Lum = 19785.0; // pb-1
  RooRealVar lumi("lumi","lumi",Lum);
  w->import(lumi);
  RooArgSet* ntplVars = defineVariables();
  TFile sigFile(signalfile);
  TTree* sigTree = (TTree*) sigFile.Get("TCVARS");
  // common preselection cut
  TString mainCut("1");
  RooDataSet sigScaled(
"sigScaled",
"dataset",
sigTree,
*ntplVars,
mainCut,
"evWeight");
  cout << "======================================================================" << endl;
  RooDataSet* sigToFit[ncat];
  TString cut0 = " && 1>0";// "&& mtot > 955 && mtot < 1150 "; //
  TString cut1 = " && 1>0";//"&& mtot > 955 && mtot < 1150 "; // "&& 1>0";//
  //
  TString cutj0 = " && 1>0";//"&& mjj_wokinfit > 90 && mjj_wokinfit < 170 "; //"&& 1>0";//
  TString cutj1 = " && 1>0";//"&& mjj_wokinfit > 100 && mjj_wokinfit < 160 "; // "&& 1>0";//
  //
  // we take only mtot to fit to the workspace, we include the cuts
  sigToFit[0] = (RooDataSet*) sigScaled.reduce(
RooArgList(*w->var("mgg"),*w->var("mjj")),
mainCut+TString::Format(" && cut_based_ct==%d ",0)+cut0+cutj0);
  w->import(*sigToFit[0],Rename(TString::Format("Sig_cat%d",0)));
    //
  sigToFit[1] = (RooDataSet*) sigScaled.reduce(
RooArgList(*w->var("mgg"),*w->var("mjj")),
mainCut+TString::Format(" && cut_based_ct==%d ",1)+cut1+cutj1);
  w->import(*sigToFit[1],Rename(TString::Format("Sig_cat%d",1)));
  // Create full signal data set without categorization
  RooDataSet* sigToFitAll = (RooDataSet*) sigScaled.reduce(RooArgList(*w->var("mgg"),*w->var("mjj")),mainCut);
  cout << "======================================================================" << endl;
  //w->import(*sigToFitAll,Rename("Sig"));
  // here we print the number of entries on the different categories
  cout << "========= the number of entries on the different categories ==========" << endl;
  cout << "---- one channel: " << sigScaled.sumEntries() << endl;
  for (int c = 0; c < ncat; ++c) {
    Float_t nExpEvt = sigToFit[c]->sumEntries();
    integral[c] = sigToFit[c]->sumEntries();
    cout << TString::Format("nEvt exp. cat%d : ",c) << nExpEvt
<< TString::Format(" eff x Acc cat%d : ",c)
<< "%"
<< endl;
  } // close ncat
  cout << "======================================================================" << endl;
  sigScaled.Print("v");
  return;
} // end add signal function
///////////////////////////////////////////////////////////////////////////////////
// we add the data to the workspace in categories
void AddBkgData(RooWorkspace* w, TString datafile) {
  const Int_t ncat = NCAT;
  // common preselection cut
  TString mainCut("1");
  RooArgSet* ntplVars = defineVariables();
  RooRealVar weightVar("weightVar","",1,0,1000);
  weightVar.setVal(1.);
  // no common preselection cut applied yet;
  TFile dataFile(datafile);
  TTree* dataTree = (TTree*) dataFile.Get("TCVARS");
  RooDataSet Data("Data","dataset",dataTree,*ntplVars,"","evWeight");
  // evweight is 1 anyway...
  RooDataSet* dataToFit[ncat];
  RooDataSet* dataToPlot[ncat];
  TString cut0 = "&& 1>0";//"&& mtot > 955 && mtot < 1150 "; //"&& 1>0";//
  TString cut1 = "&& 1>0";//"&& mtot > 955 && mtot < 1150 "; //"&& 1>0";//
  //
  TString cutj0 = "&& 1>0";//"&& mjj_wokinfit > 90 && mjj_wokinfit < 170 "; //"&& 1>0";//
  TString cutj1 = "&& 1>0";//"&& mjj_wokinfit > 100 && mjj_wokinfit < 160 "; // "&& 1>0";//
  //
  cout<<" HERE TAKE DATASET"<<endl;

dataToFit[0] = (RooDataSet*) Data.reduce(
RooArgList(*w->var("mgg"),*w->var("mjj")),
mainCut+TString::Format(" && cut_based_ct==%d",0)+cut0+cutj0);
    dataToPlot[0] = (RooDataSet*) Data.reduce(
RooArgList(*w->var("mgg"),*w->var("mjj")),
mainCut+TString::Format(" && cut_based_ct==%d",0)
//+TString::Format(" && (mgg > 130 || mgg < 120)")// blind
+cut0+cutj0); 
   
dataToFit[1] = (RooDataSet*) Data.reduce(
RooArgList(*w->var("mgg"),*w->var("mjj")),
mainCut+TString::Format(" && cut_based_ct==%d",1)+cut1);
    dataToPlot[1] = (RooDataSet*) Data.reduce(
RooArgList(*w->var("mgg"),*w->var("mjj")),
mainCut+TString::Format(" && cut_based_ct==%d",1)
//+TString::Format(" && (mgg > 130 || mgg < 120)") // blind
+cut1); 

  for (int c = 0; c < ncat; ++c) {
    w->import(*dataToFit[c],Rename(TString::Format("Data_cat%d",c)));
    w->import(*dataToPlot[c],Rename(TString::Format("Dataplot_cat%d",c)));
  }
  // Create full data set without categorization
  RooDataSet* data = (RooDataSet*) Data.reduce(RooArgList(*w->var("mgg"),*w->var("mjj")),mainCut);
  w->import(*data, Rename("Data"));
  data->Print("v");
  return;
} // close add data ..
/////////////////////////////////////////////////////////////////////
// we make the fit model
void SigModelFit(RooWorkspace* w, Float_t mass) {
  const Int_t ncat = NCAT;
  Float_t MASS(mass);
  //******************************************//
  // Fit signal with model pdfs
  //******************************************//
  // four categories to fit
  RooDataSet* sigToFit[ncat];
  RooAbsPdf* mggSig[ncat];
  RooAbsPdf*  mjjSig[9];
  // fit range
  Float_t minSigMggFit(115),maxSigMggFit(135);
  Float_t minSigMjjFit(60), maxSigMjjFit(180);
  for (int c = 0; c < ncat; ++c) {
    // import sig and data from workspace
    sigToFit[c] = (RooDataSet*) w->data(TString::Format("Sig_cat%d",c));
    mggSig[c] = (RooAbsPdf*) w->pdf(TString::Format("mggSig_cat%d",c));
      ((RooRealVar*) w->var(TString::Format("mgg_sig_m0_cat%d",c)))->setVal(MASS);
    cout << "OK up to now..." <<MASS<< endl;
    // Fit model as M(x|y) to D(x,y)
    mggSig[c]->fitTo(*sigToFit[c],Range(minSigMggFit,maxSigMggFit),SumW2Error(kTRUE));
    cout << "old = " << ((RooRealVar*) w->var(TString::Format("mgg_sig_m0_cat%d",c)))->getVal() << endl;
    double mPeak = ((RooRealVar*) w->var(TString::Format("mgg_sig_m0_cat%d",c)))->getVal()+0.0; // shift the peak
    ((RooRealVar*) w->var(TString::Format("mgg_sig_m0_cat%d",c)))->setVal(mPeak); // shift the peak
    cout << "mPeak = " << mPeak << endl;
    cout << "new mPeak position = " << ((RooRealVar*) w->var(TString::Format("mgg_sig_m0_cat%d",c)))->getVal() << endl;
    // IMPORTANT: fix all pdf parameters to constant, why?
    w->defineSet(TString::Format("SigMggParam_cat%d",c),
        RooArgSet(
	*w->var(TString::Format("mgg_sig_m0_cat%d",c)),
	*w->var(TString::Format("mgg_sig_sigma_cat%d",c)),
	*w->var(TString::Format("mgg_sig_alpha_cat%d",c)),
	*w->var(TString::Format("mgg_sig_n_cat%d",c)),
	*w->var(TString::Format("mgg_sig_gsigma_cat%d",c)),
	*w->var(TString::Format("mgg_sig_frac_cat%d",c))) );
    SetConstantParams(w->set(TString::Format("SigMggParam_cat%d",c)));
    mjjSig[c] = (RooAbsPdf*) w->pdf(TString::Format("mjjSig_cat%d",c));
    cout << "OK up to now..." <<MASS<< endl;
    // Fit model as M(x|y) to D(x,y)
    mjjSig[c]->fitTo(*sigToFit[c],Range(minSigMjjFit,maxSigMjjFit),SumW2Error(kTRUE));
    // IMPORTANT: fix all pdf parameters to constant, why?
    w->defineSet(TString::Format("SigMjjParam_cat%d",c),
        RooArgSet(
	*w->var(TString::Format("mjj_sig_m0_cat%d",c)),
	*w->var(TString::Format("mjj_sig_sigma_cat%d",c)),
	*w->var(TString::Format("mjj_sig_alpha_cat%d",c)),
	*w->var(TString::Format("mjj_sig_n_cat%d",c)),
	*w->var(TString::Format("mjj_sig_gsigma_cat%d",c)),
	*w->var(TString::Format("mjj_sig_frac_cat%d",c))) );
    SetConstantParams(w->set(TString::Format("SigMjjParam_cat%d",c)));
    RooProdPdf SigPdf(TString::Format("SigPdf_cat%d",c),"",RooArgSet(*mggSig[c], *mjjSig[c]));    
    SetConstantParams(w->set(TString::Format("SigMggParam_cat%d",c)));
    SetConstantParams(w->set(TString::Format("SigMjjParam_cat%d",c)));
    w->import(SigPdf);
  } // close for ncat
} // close signal model fit
/////////////////////////////////////////
// we make the higgs model
void HigModelFit(RooWorkspace* w, Float_t mass, int higgschannel) {
  const Int_t ncat = NCAT;
  Float_t MASS(mass);
  // four categories to fit
  RooDataSet* higToFit[ncat];
  RooAbsPdf* mggHig[ncat];
  RooAbsPdf* mjjHig[ncat];
  // fit range
  Float_t minSigFit(115),maxSigFit(135);
  for (int c = 0; c < ncat; ++c) {
    // import sig and data from workspace
    higToFit[c] = (RooDataSet*) w->data(TString::Format("Hig_%d_cat%d",higgschannel,c));
    mggHig[c] = (RooAbsPdf*) w->pdf(TString::Format("mggHig_%d_cat%d",higgschannel,c));
    ((RooRealVar*) w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)))->setVal(MASS);
    cout << "OK up to now..." <<MASS<< endl;
    cout << "old = " << ((RooRealVar*) w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)))->getVal() << endl;
    double mPeak = ((RooRealVar*) w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)))->getVal()+0.0; // shift the peak
    ((RooRealVar*) w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)))->setVal(mPeak); // shift the peak
    cout << "mPeak = " << mPeak << endl;
    cout << "new mPeak position = " << ((RooRealVar*) w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)))->getVal() << endl;
    // IMPORTANT: fix all pdf parameters to constant
    w->defineSet(TString::Format("HigPdfParam_%d_cat%d",higgschannel,c),
        RooArgSet(
	*w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)),
	*w->var(TString::Format("mgg_hig_sigma_%d_cat%d",higgschannel,c)),
	*w->var(TString::Format("mgg_hig_alpha_%d_cat%d",higgschannel,c)),
	*w->var(TString::Format("mgg_hig_n_%d_cat%d",higgschannel,c)),
	*w->var(TString::Format("mgg_hig_gsigma_%d_cat%d",higgschannel,c)),
	*w->var(TString::Format("mgg_hig_frac_%d_cat%d",higgschannel,c))) );
    SetConstantParams(w->set(TString::Format("HigPdfParam_%d_cat%d",higgschannel,c)));

  } // close for ncat
  cout<<" ---- here fittted higs  -----"<< endl;
} // close higgs model fit
////////////////////////////////////////////////////////////
// BKG model berestein 3
RooFitResult* BkgModelFitBernstein(RooWorkspace* w, Bool_t dobands) {
  const Int_t ncat = NCAT;
  std::vector<TString> catdesc;
  catdesc.push_back("High purity");//2 btag");
  catdesc.push_back("Medium purity");//1 btag");
  catdesc.push_back("cat 2");
  catdesc.push_back("cat 3");
  //******************************************//
  // Fit background with model pdfs
  //******************************************//
  // retrieve pdfs and datasets from workspace to fit with pdf models
  RooDataSet* data[ncat];
  RooDataSet* dataplot[ncat]; // the data
  RooBernstein* mggBkg[ncat];// the polinomial of 4* order
  RooBernstein* mjjBkg[ncat];// the polinomial of 4* order
   RooGaussian* Higgs[ncat]; // the higgs to sum
  RooPlot* plotmggBkg[ncat];
  RooPlot* plotmjjBkg[ncat];
  RooDataSet* sigToFit0[ncat];
  RooDataSet* sigToFit1[ncat];
  RooDataSet* sigToFit2[ncat];
  RooDataSet* sigToFit3[ncat];
  RooAbsPdf* mggSig[ncat];
  Float_t minMassFitMgg(100),maxMassFitMgg(180);
  Float_t minMassFitMjj(60), maxMassFitMjj(180);
  // Fit data with background pdf for data limit
  RooRealVar* mgg = w->var("mgg");
  mgg->setUnit("GeV");
  mgg->setRange("MassRange",minMassFitMgg,maxMassFitMgg);
  RooRealVar* mjj = w->var("mjj");
  mjj->setUnit("GeV");
  mjj->setRange("MassRange",minMassFitMjj,maxMassFitMjj);
  //
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  //
  for (int c = 0; c < ncat; ++c) { // to each category
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    cout << "!!!!!!!!!!!!!" << endl;
    ////////////////////////////////////
    // these are the parameters for the bkg polinomial
    // one slope by category - float from -10 > 10
    // the parameters are squared
    RooFormulaVar *mgg_p1mod = new RooFormulaVar(TString::Format("mgg_p1mod_cat%d",c),"","@0*@0",*w->var(TString::Format("mgg_bkg_8TeV_slope1_cat%d",c)));
    //RooFormulaVar *mgg_p2mod = new RooFormulaVar(TString::Format("mgg_p2mod_cat%d",c),"","@0*@0",*w->var(TString::Format("mgg_bkg_8TeV_slope2_cat%d",c)));
    //RooFormulaVar *mgg_p3mod = new RooFormulaVar(TString::Format("mgg_p3mod_cat%d",c),"","@0*@0",*w->var(TString::Format("mgg_bkg_8TeV_slope3_cat%d",c)));
    //RooFormulaVar *mgg_p4mod = new RooFormulaVar(TString::Format("mgg_p4mod_cat%d",c),"","@0*@0",*w->var(TString::Format("mgg_bkg_8TeV_slope4_cat%d",c)));
    //RooFormulaVar *mgg_p5mod = new RooFormulaVar(TString::Format("mgg_p5mod_cat%d",c),"","@0*@0",*w->var(TString::Format("mgg_bkg_8TeV_slope5_cat%d",c)));
    RooFormulaVar *mjj_p1mod = new RooFormulaVar(TString::Format("mjj_p1mod_cat%d",c),"","@0*@0",*w->var(TString::Format("mjj_bkg_8TeV_slope1_cat%d",c)));
    RooFormulaVar *mjj_p2mod = new RooFormulaVar(TString::Format("mjj_p2mod_cat%d",c),"","@0*@0",*w->var(TString::Format("mjj_bkg_8TeV_slope2_cat%d",c)));
    RooFormulaVar *mjj_p3mod = new RooFormulaVar(TString::Format("mjj_p3mod_cat%d",c),"","@0*@0",*w->var(TString::Format("mjj_bkg_8TeV_slope3_cat%d",c)));
    RooFormulaVar *mjj_p4mod = new RooFormulaVar(TString::Format("mjj_p4mod_cat%d",c),"","@0*@0",*w->var(TString::Format("mjj_bkg_8TeV_slope4_cat%d",c)));
    RooFormulaVar *mjj_p5mod = new RooFormulaVar(TString::Format("mjj_p5mod_cat%d",c),"","@0*@0",*w->var(TString::Format("mjj_bkg_8TeV_slope5_cat%d",c)));

    RooAbsPdf* mggBkgTmp0 = new RooGenericPdf( // if exp function
                TString::Format("mggBkgTmp0_%d",c),
                "1./pow(@0,@1)",
                RooArgList(*mgg, *mgg_p1mod));
    RooAbsPdf* mjjBkgTmp0 = 0;
    if(c == 0)
      //mjjBkgTmp0 = new  RooBernstein(TString::Format("mjjBkgTmp0_cat%d",c), "", *mjj, RooArgList(RooConst(1.0),
      //											 *mjj_p1mod, *mjj_p2mod, *mjj_p3mod, *mjj_p4mod, *mjj_p5mod));
      mjjBkgTmp0 = new RooGenericPdf( TString::Format("mjjBkgTmp0_cat%d",c), "1./pow(@0,@1)", RooArgList(*mjj,*mjj_p1mod));
    if(c == 1)
      mjjBkgTmp0 = new RooGenericPdf( TString::Format("mjjBkgTmp0_cat%d",c), "1./pow(@0,@1)", RooArgList(*mjj,*mjj_p1mod));
    if(c == 2)
      mjjBkgTmp0 = new  RooBernstein(TString::Format("mjjBkgTmp0_cat%d",c), "", *mjj, RooArgList(RooConst(1.0),
												  *mjj_p1mod, *mjj_p2mod, *mjj_p3mod, *mjj_p4mod, *mjj_p5mod));

    RooProdPdf     BkgPdfTmp(TString::Format("BkgPdfTmp%d",c), "Background Pdf", RooArgList(*mggBkgTmp0, *mjjBkgTmp0));
    w->factory(TString::Format("bkg_8TeV_norm_cat%d[800.0,0.0,100000]",c));
    RooExtendPdf BkgPdf(TString::Format("BkgPdf_cat%d",c),"", BkgPdfTmp,*w->var(TString::Format("bkg_8TeV_norm_cat%d",c)));
    fitresult[c] = BkgPdf.fitTo(*data[c], Strategy(1),Minos(kFALSE), Range("MassRange"),SumW2Error(kTRUE), Save(kTRUE));
    w->import(BkgPdf);

   //************************************************//
   // Plot mgg background fit results per categories
   //************************************************//
   TCanvas* ctmp = new TCanvas("ctmp","mgg Background Categories",0,0,500,500);
   Int_t nBinsMass(80);
   plotmggBkg[c] = mgg->frame(nBinsMass);
    cout<<" here 1"<<endl;
   dataplot[c] = (RooDataSet*) w->data(TString::Format("Dataplot_cat%d",c));
   cout<<" here 1"<<endl;
   data[c]->plotOn(plotmggBkg[c],LineColor(kWhite),MarkerColor(kWhite)); //
   mggBkgTmp0.plotOn(
plotmggBkg[c],
LineColor(kBlue),
Range("fitrange"),NormRange("fitrange"));
    dataplot[c]->plotOn(plotmggBkg[c]);

    cout << "!!!!!!!!!!!!!!!!!" << endl;
    cout << "!!!!!!!!!!!!!!!!!" << endl; // now we fit the gaussian on signal
    //plotmggBkg[c]->SetMinimum(0.01); // no error bar in bins with zero events
      if(c==0)plotmggBkg[c]->SetMinimum(0.005); // no error bar in bins with zero events
      if(c==1)plotmggBkg[c]->SetMinimum(0.001); // no error bar in bins with zero events
    plotmggBkg[c]->Draw();
    //plotmggBkg[c]->SetTitle("CMS preliminary 19.7/fb");
    //////////////////////////////////////////////////////////////////
  TPaveText *pt = new TPaveText(0.2,0.93,0.9,0.99, "brNDC");
  //   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   //   pt->SetShadowColor(kWhite);
   pt->AddText("            CMS Preliminary                     L = 19.7 fb^{-1}    #sqrt{s} = 8 TeV   ");
   pt->SetTextSize(0.035);
   pt->Draw();
    ////////////////////////////////////////////////////////////////////
    //plotmggBkg[c]->SetMinimum(0.01); // no error bar in bins with zero events
    plotmggBkg[c]->SetMaximum(1.40*plotmggBkg[c]->GetMaximum());
    plotmggBkg[c]->GetXaxis()->SetTitle("M_{#gamma#gamma} (GeV)");
    //double test = sigToFit[c]->sumEntries();
    //cout<<"number of events on dataset "<<test<<endl;
  if (dobands) {
      RooAbsPdf *cpdf; cpdf = mggBkgTmp0;
      TGraphAsymmErrors *onesigma = new TGraphAsymmErrors();
      TGraphAsymmErrors *twosigma = new TGraphAsymmErrors();
      RooRealVar *nlim = new RooRealVar(TString::Format("nlim%d",c),"",0.0,0.0,10.0);
      nlim->removeRange();
      RooCurve *nomcurve = dynamic_cast<RooCurve*>(plotmggBkg[c]->getObject(1));
      for (int i=1; i<(plotmggBkg[c]->GetXaxis()->GetNbins()+1); ++i) {
        double lowedge = plotmggBkg[c]->GetXaxis()->GetBinLowEdge(i);
        double upedge = plotmggBkg[c]->GetXaxis()->GetBinUpEdge(i);
        double center = plotmggBkg[c]->GetXaxis()->GetBinCenter(i);
        double nombkg = nomcurve->interpolate(center);
        nlim->setVal(nombkg);
        mgg->setRange("errRange",lowedge,upedge);
        RooAbsPdf *epdf = 0;
        epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
        RooAbsReal *nll = epdf->createNLL(*(data[c]),Extended());
        RooMinimizer minim(*nll);
        minim.setStrategy(0);
        double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
        double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
        minim.migrad();
        minim.minos(*nlim);
        // printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
        onesigma->SetPoint(i-1,center,nombkg);
        onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
        minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2));
        // the 0.5 is because qmu is -2*NLL
        // eventually if cl = 0.95 this is the usual 1.92!
        minim.migrad();
        minim.minos(*nlim);
        twosigma->SetPoint(i-1,center,nombkg);
        twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
        delete nll;
        delete epdf;
         } // close for bin
      mgg->setRange("errRange",minMassFitMgg,maxMassFitMgg);
      twosigma->SetLineColor(kYellow);
      twosigma->SetFillColor(kYellow);
      twosigma->SetMarkerColor(kYellow);
      twosigma->Draw("L3 SAME");
      onesigma->SetLineColor(kGreen);
      onesigma->SetFillColor(kGreen);
      onesigma->SetMarkerColor(kGreen);
      onesigma->Draw("L3 SAME");
      plotmggBkg[c]->Draw("SAME");
    } else plotmggBkg[c]->Draw("SAME"); // close dobands
   //plotmggBkg[c]->getObject(1)->Draw("SAME");
   //plotmggBkg[c]->getObject(2)->Draw("P SAME");
/*   ////////////////////////////////////////////////////////// plot higgs
	sigToFit0[c] = (RooDataSet*) w->data(TString::Format("Hig_0_cat%d",c));
	double norm0; norm0 = 1.0*sigToFit0[c]->sumEntries(); //
        //norm0 = 0.0000001;
	mggSig0[c] = (RooAbsPdf*) w->pdf(TString::Format("mggHig_0_cat%d",c));
	// we are not constructing signal pdf, this is constructed on sig to fit function...
	mggSig0[c] ->plotOn(
	plotmggBkg[c],
	Normalization(norm0,RooAbsPdf::NumEvent),
	DrawOption("F"),
	LineColor(kRed),FillStyle(1001),FillColor(19));
	mggSig0[c]->plotOn(
	plotmggBkg[c],
	Normalization(norm0,RooAbsPdf::NumEvent),LineColor(kRed),LineStyle(1));
	//
	sigToFit1[c] = (RooDataSet*) w->data(TString::Format("Hig_1_cat%d",c));
	double norm1 = 1.0*sigToFit1[c]->sumEntries(); //
	mggSig1[c] = (RooAbsPdf*) w->pdf(TString::Format("mggHig_1_cat%d",c));
	// we are not constructing signal pdf, this is constructed on sig to fit function...
	mggSig1[c] ->plotOn(
	plotmggBkg[c],
	Normalization(norm1,RooAbsPdf::NumEvent),
	DrawOption("F"),
	LineColor(kGreen),FillStyle(1001),FillColor(19));
	mggSig1[c]->plotOn(
	plotmggBkg[c],
	Normalization(norm1,RooAbsPdf::NumEvent),LineColor(kGreen),LineStyle(1));
	//
	sigToFit2[c] = (RooDataSet*) w->data(TString::Format("Hig_2_cat%d",c));
        double norm2;
	//if(sigToFit2[c]->sumEntries()>0) 
        norm2 = 1.0*sigToFit2[c]->sumEntries(); //else 
        //norm2 = 0.0000000000001; //
	mggSig2[c] = (RooAbsPdf*) w->pdf(TString::Format("mggHig_2_cat%d",c));
	// we are not constructing signal pdf, this is constructed on sig to fit function...
	mggSig2[c] ->plotOn(
	plotmggBkg[c],
	Normalization(norm2,RooAbsPdf::NumEvent),
	DrawOption("F"),
	LineColor(kMagenta),FillStyle(1001),FillColor(19));
	mggSig2[c]->plotOn(
	plotmggBkg[c],
	Normalization(norm2,RooAbsPdf::NumEvent),LineColor(kMagenta),LineStyle(1));
	//
	sigToFit3[c] = (RooDataSet*) w->data(TString::Format("Hig_3_cat%d",c));
	double norm3 = 1.0*sigToFit1[c]->sumEntries(); //
	mggSig3[c] = (RooAbsPdf*) w->pdf(TString::Format("mggHig_3_cat%d",c));
	// we are not constructing signal pdf, this is constructed on sig to fit function...
	mggSig3[c] ->plotOn(
	plotmggBkg[c],
	Normalization(norm3,RooAbsPdf::NumEvent),
	DrawOption("F"),
	LineColor(kCyan),FillStyle(1001),FillColor(19));
	mggSig3[c]->plotOn(
	plotmggBkg[c],
	Normalization(norm3,RooAbsPdf::NumEvent),LineColor(kCyan),LineStyle(1));
   //////////////////////////////////////////////////////////
*/
   plotmggBkg[c]->Draw("SAME"); 
      if(c==0)plotmggBkg[c]->SetMinimum(0.005); // no error bar in bins with zero events
      if(c==1)plotmggBkg[c]->SetMinimum(0.01); // no error bar in bins with zero events
      if(c==0)plotmggBkg[c]->SetMaximum(5.3); // 260 && 270 mass-points
      //if(c==0)plotmggBkg[c]->SetMaximum(6.3); // 300 && 350 mass-points
      if(c==1)plotmggBkg[c]->SetMaximum(20); // 
     // plotmggBkg[c]->SetMinimum(0.005); // no error bar in bins with zero events
    //plotmggBkg[c]->SetLogy(0);
    plotmggBkg[c]->GetYaxis()->SetTitleSize(0.06);
    plotmggBkg[c]->GetYaxis()->SetTitleOffset(1.4);
    cout << "!!!!!!!!!!!!!!!!!" << endl;
    TLegend *legmcH = new TLegend(0.40,0.72,0.62,0.9);
    TLegend *legmc = new TLegend(0.64,0.69,0.94,0.89);
    legmc->AddEntry(plotmggBkg[c]->getObject(2),"Data ","LPE"); // not...
    legmc->AddEntry(plotmggBkg[c]->getObject(1),"Fit","L");
    if(dobands)legmc->AddEntry(twosigma,"Fit #pm 2 #sigma","F"); // not...
    if(dobands)legmc->AddEntry(onesigma,"Fit #pm 1 #sigma","F");
    legmcH->AddEntry(plotmggBkg[c]->getObject(3),"ggH ","LPE"); // not...
    legmcH->AddEntry(plotmggBkg[c]->getObject(5),"ttH ","LPE"); // not...
    legmcH->AddEntry(plotmggBkg[c]->getObject(7),"VBF ","LPE"); // not...
    legmcH->AddEntry(plotmggBkg[c]->getObject(9),"VH ","LPE"); // not...
    legmc->SetHeader(" m_{X} = SM GeV");
    legmcH->SetHeader(" Higgs");
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();
//    legmcH->Draw();
    TLatex *lat2 = new TLatex(minMassFitMgg+2.3,0.8*plotmggBkg[c]->GetMaximum(),catdesc.at(c));
    lat2->Draw();
    TLatex *lat3 = new TLatex(minMassFitMgg+2.3,0.9*plotmggBkg[c]->GetMaximum(),"X #rightarrow HH #rightarrow #gamma#gammab#bar{b}");
    //lat3->SetTextSize(0.035);
    lat3->Draw();

    ctmp->SaveAs(TString::Format("databkgoversigMgg_cat%d.pdf",c));
    ctmp->SaveAs(TString::Format("databkgoversigMgg_cat%d.png",c));
    ctmp->SaveAs(TString::Format("databkgoversigMgg_cat%d.root",c));

      if(c==0)plotmggBkg[c]->SetMaximum(100); // no error bar in bins with zero events
      if(c==1)plotmggBkg[c]->SetMaximum(1000); // no error bar in bins with zero events
  //  ctmp->SetLogy(1);
  //  ctmp->SaveAs(TString::Format("databkgoversig_cat%d_log.pdf",c));
  //  ctmp->SaveAs(TString::Format("databkgoversig_cat%d_log.png",c));
  //  ctmp->SaveAs(TString::Format("databkgoversig_cat%d_log.root",c));
  // ctmp->SaveAs(TString::Format("databkgoversig_cat%d.C",c));

   //************************************************//
   // Plot mjj background fit results per categories
   //************************************************//
   ctmp = new TCanvas("ctmp2","mjj Background Categories",0,0,500,500);
   nBinsMass = 60;
   plotmjjBkg[c] = mjj->frame(nBinsMass);
   dataplot[c] = (RooDataSet*) w->data(TString::Format("Dataplot_cat%d",c));
   data[c]->plotOn(plotmjjBkg[c],LineColor(kWhite),MarkerColor(kWhite)); //
   mjjBkgTmp0.plotOn(plotmjjBkg[c],LineColor(kBlue),Range("fitrange"),NormRange("fitrange"));
    dataplot[c]->plotOn(plotmjjBkg[c]);

      if(c==0)plotmjjBkg[c]->SetMinimum(0.005); // no error bar in bins with zero events
      if(c==1)plotmjjBkg[c]->SetMinimum(0.001); // no error bar in bins with zero events
    plotmjjBkg[c]->Draw();

    pt = new TPaveText(0.2,0.93,0.9,0.99, "brNDC");
    //   pt->SetName("title");
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    //   pt->SetShadowColor(kWhite);
    pt->AddText("            CMS Preliminary                     L = 19.7 fb^{-1}    #sqrt{s} = 8 TeV   ");
    pt->SetTextSize(0.035);
    pt->Draw();
    plotmjjBkg[c]->SetMaximum(1.40*plotmjjBkg[c]->GetMaximum());
    plotmjjBkg[c]->GetXaxis()->SetTitle("M_{jj} (GeV)");
    //double test = sigToFit[c]->sumEntries();
    //cout<<"number of events on dataset "<<test<<endl;
  if (dobands) {
      RooAbsPdf *cpdf; cpdf = mjjBkgTmp0;
      TGraphAsymmErrors *onesigma = new TGraphAsymmErrors();
      TGraphAsymmErrors *twosigma = new TGraphAsymmErrors();
      RooRealVar *nlim = new RooRealVar(TString::Format("nlim%d",c),"",0.0,0.0,10.0);
      nlim->removeRange();
      RooCurve *nomcurve = dynamic_cast<RooCurve*>(plotmjjBkg[c]->getObject(1));
      for (int i=1; i<(plotmjjBkg[c]->GetXaxis()->GetNbins()+1); ++i) {
        double lowedge = plotmjjBkg[c]->GetXaxis()->GetBinLowEdge(i);
        double upedge = plotmjjBkg[c]->GetXaxis()->GetBinUpEdge(i);
        double center = plotmjjBkg[c]->GetXaxis()->GetBinCenter(i);
        double nombkg = nomcurve->interpolate(center);
        nlim->setVal(nombkg);
        mjj->setRange("errRange",lowedge,upedge);
        RooAbsPdf *epdf = 0;
        epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
        RooAbsReal *nll = epdf->createNLL(*(data[c]),Extended());
        RooMinimizer minim(*nll);
        minim.setStrategy(0);
        double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
        double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
        minim.migrad();
        minim.minos(*nlim);
        // printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
        onesigma->SetPoint(i-1,center,nombkg);
        onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
        minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2));
        // the 0.5 is because qmu is -2*NLL
        // eventually if cl = 0.95 this is the usual 1.92!
        minim.migrad();
        minim.minos(*nlim);
        twosigma->SetPoint(i-1,center,nombkg);
        twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
        delete nll;
        delete epdf;
      } // close for bin
      mjj->setRange("errRange",minMassFitMgg,maxMassFitMgg);
      twosigma->SetLineColor(kYellow);
      twosigma->SetFillColor(kYellow);
      twosigma->SetMarkerColor(kYellow);
      twosigma->Draw("L3 SAME");
      onesigma->SetLineColor(kGreen);
      onesigma->SetFillColor(kGreen);
      onesigma->SetMarkerColor(kGreen);
      onesigma->Draw("L3 SAME");
      plotmjjBkg[c]->Draw("SAME");
  } else plotmjjBkg[c]->Draw("SAME"); // close dobands
   plotmjjBkg[c]->Draw("SAME"); 
      if(c==0)plotmjjBkg[c]->SetMinimum(0.005); // no error bar in bins with zero events
      if(c==1)plotmjjBkg[c]->SetMinimum(0.01); // no error bar in bins with zero events
      if(c==0)plotmjjBkg[c]->SetMaximum(5.3); // 260 && 270 mass-points
      //if(c==0)plotmggBkg[c]->SetMaximum(6.3); // 300 && 350 mass-points
      if(c==1)plotmjjBkg[c]->SetMaximum(20); // 
     // plotmggBkg[c]->SetMinimum(0.005); // no error bar in bins with zero events
    //plotmggBkg[c]->SetLogy(0);
    plotmjjBkg[c]->GetYaxis()->SetTitleSize(0.06);
    plotmjjBkg[c]->GetYaxis()->SetTitleOffset(1.4);
    cout << "!!!!!!!!!!!!!!!!!" << endl;
    legmcH = new TLegend(0.40,0.72,0.62,0.9);
    legmc = new TLegend(0.64,0.69,0.94,0.89);
    legmc->AddEntry(plotmjjBkg[c]->getObject(2),"Data ","LPE"); // not...
    legmc->AddEntry(plotmjjBkg[c]->getObject(1),"Fit","L");
    if(dobands)legmc->AddEntry(twosigma,"Fit #pm 2 #sigma","F"); // not...
    if(dobands)legmc->AddEntry(onesigma,"Fit #pm 1 #sigma","F");
    legmcH->AddEntry(plotmjjBkg[c]->getObject(3),"ggH ","LPE"); // not...
    legmcH->AddEntry(plotmjjBkg[c]->getObject(5),"ttH ","LPE"); // not...
    legmcH->AddEntry(plotmjjBkg[c]->getObject(7),"VBF ","LPE"); // not...
    legmcH->AddEntry(plotmjjBkg[c]->getObject(9),"VH ","LPE"); // not...
    legmc->SetHeader(" m_{X} = SM GeV");
    legmcH->SetHeader(" Higgs");
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();
//    legmcH->Draw();
    TLatex *lat2 = new TLatex(minMassFitMjj+4.6,0.8*plotmjjBkg[c]->GetMaximum(),catdesc.at(c));
    lat2->Draw();
    TLatex *lat3 = new TLatex(minMassFitMjj+4.6,0.9*plotmjjBkg[c]->GetMaximum(),"X #rightarrow HH #rightarrow #gamma#gammab#bar{b}");
    //lat3->SetTextSize(0.035);
    lat3->Draw();

    ctmp->SaveAs(TString::Format("databkgoversigMjj_cat%d.pdf",c));
    ctmp->SaveAs(TString::Format("databkgoversigMjj_cat%d.png",c));
    ctmp->SaveAs(TString::Format("databkgoversigMjj_cat%d.root",c));

      if(c==0)plotmjjBkg[c]->SetMaximum(100); // no error bar in bins with zero events
      if(c==1)plotmjjBkg[c]->SetMaximum(1000); // no error bar in bins with zero event

  } // close to each category

  RooBernstein *mggBkgAll = new RooBernstein("mggBkgAll", "", *mgg,RooArgList(RooConst(1.0),*w->var("mgg_bkg_8TeV_slope1"),*w->var("mgg_bkg_8TeV_slope2"),*w->var("mgg_bkg_8TeV_slope3")));
  RooBernstein *mjjBkgAll = new RooBernstein("mjjBkgAll", "", *mjj,RooArgList(RooConst(1.0),*w->var("mjj_bkg_8TeV_slope1"),*w->var("mjj_bkg_8TeV_slope2"),
							  *w->var("mjj_bkg_8TeV_slope3"),*w->var("mjj_bkg_8TeV_slope4"),*w->var("mjj_bkg_8TeV_slope5")));
  RooProdPdf     BkgPdfAll("BkgPdfAll", "Background Pdf", *mggBkgAll, *mjjBkgAll);
  RooFitResult* fitresults = BkgPdfAll.fitTo( // save results to workspace
*w->data("Data"),
Range("MassRange"),
SumW2Error(kTRUE), Save(kTRUE));
  fitresults->Print();
  w->import(BkgPdfAll);
  return fitresults;
} // close berestein 3
///////////////////////////////////////////////////////////////
void MakeSigWS(RooWorkspace* w, const char* fileBaseName) {
  TString wsDir = "workspaces/";
  const Int_t ncat = NCAT;
  //**********************************************************************//
  // Write pdfs and datasets into the workspace before to save
  // for statistical tests.
  //**********************************************************************//
  RooAbsPdf* SigPdf[ncat];
  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");
  for (int c = 0; c < ncat; ++c) {
    SigPdf[c] = (RooAbsPdf*) w->pdf(TString::Format("SigPdf_cat%d",c));
    wAll->import(*w->pdf(TString::Format("SigPdf_cat%d",c)));
  }
  // (2) Systematics on energy scale and resolution
  // 1,1,1 statistical to be treated on the datacard
  wAll->factory("CMS_hgg_sig_m0_absShift[1,1,1]");
  wAll->factory("prod::CMS_hgg_sig_m0_cat0(mgg_sig_m0_cat0, CMS_hgg_sig_m0_absShift)");
  wAll->factory("prod::CMS_hgg_sig_m0_cat1(mgg_sig_m0_cat1, CMS_hgg_sig_m0_absShift)");
  wAll->factory("CMS_hbb_sig_m0_absShift[1,1,1]");
  wAll->factory("prod::CMS_hbb_sig_m0_cat0(mjj_sig_m0_cat0, CMS_hbb_sig_m0_absShift)");
  wAll->factory("prod::CMS_hbb_sig_m0_cat1(mjj_sig_m0_cat1, CMS_hbb_sig_m0_absShift)");
  // (3) Systematics on resolution
  wAll->factory("CMS_hgg_sig_sigmaScale[1,1,1]");
  wAll->factory("prod::CMS_hgg_sig_sigma_cat0(mgg_sig_sigma_cat0, CMS_hgg_sig_sigmaScale)");
  wAll->factory("prod::CMS_hgg_sig_sigma_cat1(mgg_sig_sigma_cat1, CMS_hgg_sig_sigmaScale)");
  wAll->factory("prod::CMS_hgg_sig_gsigma_cat0(mgg_sig_gsigma_cat0, CMS_hgg_sig_sigmaScale)");
  wAll->factory("prod::CMS_hgg_sig_gsigma_cat1(mgg_sig_gsigma_cat1, CMS_hgg_sig_sigmaScale)");
  wAll->factory("CMS_hbb_sig_sigmaScale[1,1,1]");
  wAll->factory("prod::CMS_hbb_sig_sigma_cat0(mjj_sig_sigma_cat0, CMS_hbb_sig_sigmaScale)");
  wAll->factory("prod::CMS_hbb_sig_sigma_cat1(mjj_sig_sigma_cat1, CMS_hbb_sig_sigmaScale)");
  wAll->factory("prod::CMS_hbb_sig_gsigma_cat0(mjj_sig_gsigma_cat0, CMS_hbb_sig_sigmaScale)");
  wAll->factory("prod::CMS_hbb_sig_gsigma_cat1(mjj_sig_gsigma_cat1, CMS_hbb_sig_sigmaScale)");
  // save the other parameters
/* for (int c = 0; c < ncat; ++c) {
wAll->factory(
TString::Format("CMS_hgg_sig_alpha_cat%d[%g,0.5,5]",
c, wAll->var(TString::Format("mgg_sig_alpha_cat%d",c))->getVal()));
wAll->factory(
TString::Format("CMS_hgg_sig_n_cat%d[%g,0.5,20]",
c, wAll->var(TString::Format("mgg_sig_n_cat%d",c))->getVal()));
wAll->factory(
TString::Format("CMS_hgg_sig_frac_cat%d[%g,0.0,1.0]",
c, wAll->var(TString::Format("mgg_sig_frac_cat%d",c))->getVal()));
}
*/
  // (4) do reparametrization of signal
  for (int c = 0; c < ncat; ++c){
    wAll->factory(TString::Format("EDIT::CMS_sig_cat%d(SigPdf_cat%d,",c,c) + 
		  TString::Format(" mgg_sig_m0_cat%d=CMS_hgg_sig_m0_cat%d, ", c,c) +
		  TString::Format(" mgg_sig_sigma_cat%d=CMS_hgg_sig_sigma_cat%d, ", c,c) + 
		  TString::Format(" mgg_sig_gsigma_cat%d=CMS_hgg_sig_gsigma_cat%d, ", c,c) + 
		  TString::Format(" mjj_sig_m0_cat%d=CMS_hbb_sig_m0_cat%d, ", c,c) +
		  TString::Format(" mjj_sig_sigma_cat%d=CMS_hbb_sig_sigma_cat%d, ", c,c) + 
		  TString::Format(" mjj_sig_gsigma_cat%d=CMS_hbb_sig_gsigma_cat%d)", c,c)  );
  }
  TString filename(wsDir+TString(fileBaseName)+".inputsig.root");
  wAll->writeToFile(filename);
  cout << "Write signal workspace in: " << filename << " file" << endl;
  return;
} // close make signal WP
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void MakeBkgWS(RooWorkspace* w, const char* fileBaseName) {
  TString wsDir = "workspaces/";
  const Int_t ncat = NCAT;

  //**********************************************************************//
  // Write pdfs and datasets into the workspace before to save to a file
  // for statistical tests.
  //**********************************************************************//
  RooDataSet* data[ncat];
  RooProdPdf* BkgPdf[ncat];
  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");
  for (int c = 0; c < ncat; ++c) {
    cout<<"here "<<endl;
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    BkgPdf[c] = (RooProdPdf*) w->pdf(TString::Format("BkgPdf_cat%d",c));
    wAll->import(*data[c], Rename(TString::Format("data_obs_cat%d",c)));
    wAll->import(*w->pdf(TString::Format("BkgPdf_cat%d",c)));
    wAll->factory(TString::Format("CMS_bkg_8TeV_cat%d_norm[%g,0.0,100000.0]",c, wAll->var(TString::Format("bkg_8TeV_norm_cat%d",c))->getVal()));
    wAll->factory(TString::Format("CMS_hgg_bkg_8TeV_slope1_cat%d[%g,-10,10]",c, wAll->var(TString::Format("mgg_bkg_8TeV_slope1_cat%d",c))->getVal()));
    wAll->factory(TString::Format("CMS_hbb_bkg_8TeV_slope1_cat%d[%g,-10,10]",c, wAll->var(TString::Format("mjj_bkg_8TeV_slope1_cat%d",c))->getVal()));
    //if(c==1) wAll->factory(TString::Format("CMS_hbb_bkg_8TeV_slope2_cat%d[%g,-10,10]",c, wAll->var(TString::Format("mjj_bkg_8TeV_slope2_cat%d",c))->getVal()));
    //wAll->factory(TString::Format("CMS_hbb_bkg_8TeV_slope3_cat%d[%g,-10,10]",c, wAll->var(TString::Format("mjj_bkg_8TeV_slope3_cat%d",c))->getVal()));
    //wAll->factory(TString::Format("CMS_hbb_bkg_8TeV_slope4_cat%d[%g,-10,10]",c, wAll->var(TString::Format("mjj_bkg_8TeV_slope4_cat%d",c))->getVal()));
    //wAll->factory(TString::Format("CMS_hbb_bkg_8TeV_slope5_cat%d[%g,-10,10]",c, wAll->var(TString::Format("mjj_bkg_8TeV_slope5_cat%d",c))->getVal()));
  } // close ncat
  // (2) do reparametrization of background
  for (int c = 0; c < ncat; ++c){
      wAll->factory(TString::Format("EDIT::CMS_bkg_8TeV_cat%d(BkgPdf_cat%d, ",c,c) +
		    TString::Format("bkg_8TeV_norm_cat%d=CMS_bkg_8TeV_cat%d_norm, ", c,c)+
		    TString::Format("mgg_bkg_8TeV_slope1_cat%d=CMS_hgg_bkg_8TeV_slope1_cat%d, ", c,c)+
		    TString::Format("mjj_bkg_8TeV_slope1_cat%d=CMS_hbb_bkg_8TeV_slope1_cat%d)", c,c));
      //if more parameters in Mjj fit
      //wAll->factory(TString::Format("EDIT::CMS_bkg_8TeV_cat%d(BkgPdf_cat%d, ",c,c) +
      //	    TString::Format("bkg_8TeV_norm_cat%d=CMS_bkg_8TeV_cat%d_norm, ", c,c)+
      //	    TString::Format("mgg_bkg_8TeV_slope1_cat%d=CMS_hgg_bkg_8TeV_slope1_cat%d, ", c,c)+
      //	    TString::Format("mjj_bkg_8TeV_slope1_cat%d=CMS_hbb_bkg_8TeV_slope1_cat%d, ", c,c)+
      //	    TString::Format("mjj_bkg_8TeV_slope2_cat%d=CMS_hbb_bkg_8TeV_slope2_cat%d)", c,c));

 } // close for cat

  TString filename(wsDir+TString(fileBaseName)+".root");
  wAll->writeToFile(filename);
  cout << "Write background workspace in: " << filename << " file" << endl;
  std::cout << "observation ";
  for (int c = 0; c < ncat; ++c) {
    std::cout << " " << wAll->data(TString::Format("data_obs_cat%d",c))->sumEntries();
  }
  std::cout << std::endl;
  return;
} // close make BKG workspace
////////////////////////////////////////////////////////////////////
void SetConstantParams(const RooArgSet* params) {
  // set constant parameters for signal fit, ... NO IDEA !!!!
  TIterator* iter(params->createIterator());
  for (TObject *a = iter->Next(); a != 0; a = iter->Next()) {
    RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
    if (rrv) { rrv->setConstant(true); std::cout << " " << rrv->GetName(); }
  }
} // close set const parameters
////////////////////////////////////////////////////////////////////////
void MakePlots(RooWorkspace* w, Float_t Mass) {
  const Int_t ncat = NCAT;
  std::vector<TString> catdesc;
  catdesc.push_back("High purity ");//2 btag");
  catdesc.push_back("Medium purity ");//1 btag");
  catdesc.push_back("cat 2");
  catdesc.push_back("cat 3");
  // retrieve data sets from the workspace
  // RooDataSet* dataAll = (RooDataSet*) w->data("Data");
  //RooDataSet* signalAll = (RooDataSet*) w->data("Sig");
  //RooDataSet* higgsAll = (RooDataSet*) w->data("Hig");
  // blinded dataset
  // RooDataSet* data[ncat];
  RooDataSet* sigToFit[ncat];
  RooAbsPdf* mggGaussSig[ncat];
  RooAbsPdf* mggCBSig[ncat];
  RooAbsPdf* mggSig[ncat];
  RooAbsPdf* mjjGaussSig[ncat];
  RooAbsPdf* mjjCBSig[ncat];
  RooAbsPdf* mjjSig[ncat];
  //
  RooAbsPdf* mggBkg[ncat];
  RooAbsPdf* mjjBkg[ncat];
  for (int c = 0; c < ncat; ++c) {
  // data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    sigToFit[c] = (RooDataSet*) w->data(TString::Format("Sig_cat%d",c));
    mggGaussSig[c] = (RooAbsPdf*) w->pdf(TString::Format("mggGaussSig_cat%d",c));
    mggCBSig[c] = (RooAbsPdf*) w->pdf(TString::Format("mggCBSig_cat%d",c));
    mggSig[c] = (RooAbsPdf*) w->pdf(TString::Format("mggSig_cat%d",c));
    mggBkg[c] = (RooAbsPdf*) w->pdf(TString::Format("mggBkg_cat%d",c));
    mjjGaussSig[c] = (RooAbsPdf*) w->pdf(TString::Format("mjjGaussSig_cat%d",c));
    mjjCBSig[c] = (RooAbsPdf*) w->pdf(TString::Format("mjjCBSig_cat%d",c));
    mjjSig[c] = (RooAbsPdf*) w->pdf(TString::Format("mjjSig_cat%d",c));
    mjjBkg[c] = (RooAbsPdf*) w->pdf(TString::Format("mjjBkg_cat%d",c));
  } // close categories
  RooRealVar* mgg = w->var("mgg");
  mgg->setUnit("GeV");
  RooAbsPdf* mggGaussSigAll = w->pdf("mggGaussSig");
  RooAbsPdf* mggCBSigAll = w->pdf("mggCBSig");
  RooAbsPdf* mggSigAll = w->pdf("mggSig");
  RooRealVar* mjj = w->var("mjj");
  mjj->setUnit("GeV");
  RooAbsPdf* mjjGaussSigAll = w->pdf("mjjGaussSig");
  RooAbsPdf* mjjCBSigAll = w->pdf("mjjCBSig");
  RooAbsPdf* mjjSigAll = w->pdf("mjjSig");
  //RooAbsPdf* mggBkgAll = w->pdf("mggBkg_cat1");
  //
  //****************************//
  // Plot mgg Fit results
  //****************************//
  // Set P.D.F. parameter names
  // WARNING: Do not use it if Workspaces are created
  // SetParamNames(w);
  Float_t minSigFitMgg(120),maxSigFitMgg(130);
  Float_t minSigFitMjj(60), maxSigFitMjj(180);
  Float_t MASS(Mass);
  Int_t nBinsMass(20); // just need to plot
  //RooPlot* plotmggAll = mgg->frame(Range(minSigFit,maxSigFit),Bins(nBinsMass));
  //signalAll->plotOn(plotmggAll);
  gStyle->SetOptTitle(0);
  TCanvas* c1 = new TCanvas("c1","mgg",0,0,500,500);
  c1->cd(1);
  //********************************************//
  // Plot Signal Categories
  //****************************//
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  RooPlot* plotmgg[ncat];
  for (int c = 0; c < ncat; ++c) {
    plotmgg[c] = mgg->frame(Range(minSigFitMgg,maxSigFitMgg),Bins(nBinsMass));
    sigToFit[c]->plotOn(plotmgg[c],LineColor(kWhite),MarkerColor(kWhite),Rescale(1./integral[c]));
    mggSig[c] ->plotOn(plotmgg[c],Rescale(1./integral[c]));
    double chi2n = plotmgg[c]->chiSquare(0) ;
    cout << "------------------------- Experimentakl chi2 = " << chi2n << endl;
    mggSig[c] ->plotOn(
plotmgg[c],
Components(TString::Format("mggGaussSig_cat%d",c)),
LineStyle(kDashed),LineColor(kGreen),Rescale(1./integral[c]));
    mggSig[c] ->plotOn(
plotmgg[c],
Components(TString::Format("mggCBSig_cat%d",c)),
LineStyle(kDashed),LineColor(kRed),Rescale(1./integral[c]));
    //mggSig[c] ->paramOn(plotmgg[c]);
    sigToFit[c] ->plotOn(plotmgg[c],MarkerStyle(25),Rescale(1./integral[c]),RooFit::XErrorSize(0));
    /*TH1F* dataHist = (TH1F*)sigToFit[c] ->createHistogram("mgg");
    dataHist->Scale(1./dataHist->Integral());
    for(ii = 1; ii <= dataHist->GetNbinsX(); ii++)
        std::cout << "bin - " << ii << " - " << dataHist->GetBinContent(ii)/dataHist->Integral() << "\n" << std::endl;
   std::cout << "DataHist - " << c << " - Maximum = " << dataHist->GetMaximum()/integral[c] << std::endl;*/
// TCanvas* dummy = new TCanvas("dummy", "dummy",0, 0, 400, 400);
    TH1F *hist = new TH1F("hist", "hist", 400, minSigFitMgg, maxSigFitMgg);
    //plotmgg[c]->SetTitle("CMS preliminary 19.7/fb ");
    plotmgg[c]->SetMinimum(0.0);
    //plotmgg[c]->SetMaximum(1.*plotmgg[c]->GetMaximum());
    plotmgg[c]->GetXaxis()->SetTitle("M_{#gamma#gamma} (GeV)");
    std::string Ytitle = std::string(plotmgg[c]->GetYaxis()->GetTitle());
    Ytitle.replace(0,6,"Fraction of events");
    plotmgg[c]->GetYaxis()->SetTitle(Ytitle.c_str());
    plotmgg[c]->GetYaxis()->SetTitleSize(0.06);
    plotmgg[c]->GetYaxis()->SetTitleOffset(1.4);
    TCanvas* ctmp = new TCanvas("ctmp","Background Categories",0,0,500,500);
    plotmgg[c]->Draw();
    //plotmgg[c]->Draw("SAME");
    TLegend *legmc = new TLegend(0.65,0.6,0.95,0.99);
    legmc->AddEntry(plotmgg[c]->getObject(4),"Simulation","PE");
    legmc->AddEntry(plotmgg[c]->getObject(1),"Parametric Model","L");
    legmc->AddEntry(plotmgg[c]->getObject(2),"Gaussian ","L");
    legmc->AddEntry(plotmgg[c]->getObject(3),"Crystal Ball","L");
    legmc->SetHeader(" ");
    legmc->SetLineColor(kWhite);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();
    /////////////////////////////////////////////////////
  TPaveText *pt = new TPaveText(0.1,0.93,0.7,0.99, "brNDC");
   //pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetShadowColor(kWhite);
   pt->SetLineColor(kWhite);
   pt->AddText("CMS Preliminary Simulation ");
   pt->SetTextSize(0.035);
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetShadowColor(kWhite);
   pt->SetLineColor(kWhite);
   pt->Draw();
    ////////////////////////////////////////////////  
    // float effS = effSigma(hist);
    TLatex *lat = new TLatex(minSigFitMgg+0.5,0.7*plotmgg[c]->GetMaximum(),"m_{X} = SM GeV");
    lat->Draw();
    TLatex *lat2 = new TLatex(minSigFitMgg+0.5,0.8*plotmgg[c]->GetMaximum(),catdesc.at(c));
    lat2->Draw();
    TLatex *lat3 = new TLatex(minSigFitMgg+0.5,0.9*plotmgg[c]->GetMaximum(),"X#rightarrow HH #rightarrow #gamma#gammab#bar{b}");
    lat3->Draw();
    ///////
    char myChi2buffer[50];
    sprintf(myChi2buffer,"#chi^{2}/ndof = %f",chi2n);
    TLatex* latex = new TLatex(0.52, 0.7, myChi2buffer);
    latex -> SetNDC();
    latex -> SetTextFont(42);
    latex -> SetTextSize(0.04);
    //latex -> Draw("same");
    ctmp->SaveAs(TString::Format("sigmodelMgg_cat%d.pdf",c));
    ctmp->SaveAs(TString::Format("sigmodelMgg_cat%d.png",c));
    ctmp->SaveAs(TString::Format("sigmodelMgg_cat%d.root",c));
  }
  
  c1 = new TCanvas("c1","mjj",0,0,500,500);
  c1->cd(1);
  //********************************************//
  // Plot Signal Categories
  //****************************//
  text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  RooPlot* plotmjj[ncat];
  for (int c = 0; c < ncat; ++c) {
    plotmjj[c] = mjj->frame(Range(minSigFitMjj,maxSigFitMjj),Bins(nBinsMass));
    sigToFit[c]->plotOn(plotmjj[c],LineColor(kWhite),MarkerColor(kWhite),Rescale(1./integral[c]));
    mjjSig[c] ->plotOn(plotmjj[c],Rescale(1./integral[c]));
    double chi2n = plotmjj[c]->chiSquare(0) ;
    cout << "------------------------- Experimentakl chi2 = " << chi2n << endl;
    mjjSig[c] ->plotOn(
plotmjj[c],
Components(TString::Format("mjjGaussSig_cat%d",c)),
LineStyle(kDashed),LineColor(kGreen),Rescale(1./integral[c]));
    mjjSig[c] ->plotOn(
plotmjj[c],
Components(TString::Format("mjjCBSig_cat%d",c)),
LineStyle(kDashed),LineColor(kRed),Rescale(1./integral[c]));
    //mjjSig[c] ->paramOn(plotmjj[c]);
    sigToFit[c] ->plotOn(plotmjj[c],MarkerStyle(25),Rescale(1./integral[c]),RooFit::XErrorSize(0));
    /*TH1F* dataHist = (TH1F*)sigToFit[c] ->createHistogram("mjj");
    dataHist->Scale(1./dataHist->Integral());
    for(ii = 1; ii <= dataHist->GetNbinsX(); ii++)
        std::cout << "bin - " << ii << " - " << dataHist->GetBinContent(ii)/dataHist->Integral() << "\n" << std::endl;
   std::cout << "DataHist - " << c << " - Maximum = " << dataHist->GetMaximum()/integral[c] << std::endl;*/
// TCanvas* dummy = new TCanvas("dummy", "dummy",0, 0, 400, 400);
    TH1F *hist = new TH1F("hist", "hist", 400, minSigFitMjj, maxSigFitMjj);
    //plotmjj[c]->SetTitle("CMS preliminary 19.7/fb ");
    plotmjj[c]->SetMinimum(0.0);
    //plotmjj[c]->SetMaximum(1.*plotmjj[c]->GetMaximum());
    plotmjj[c]->GetXaxis()->SetTitle("M_{jj} (GeV)");
    std::string Ytitle = std::string(plotmjj[c]->GetYaxis()->GetTitle());
    Ytitle.replace(0,6,"Fraction of events");
    plotmjj[c]->GetYaxis()->SetTitle(Ytitle.c_str());
    plotmjj[c]->GetYaxis()->SetTitleSize(0.06);
    plotmjj[c]->GetYaxis()->SetTitleOffset(1.4);
    TCanvas* ctmp = new TCanvas("ctmp","Background Categories",0,0,500,500);
    plotmjj[c]->Draw();
    //plotmjj[c]->Draw("SAME");
    TLegend *legmc = new TLegend(0.65,0.6,0.95,0.99);
    legmc->AddEntry(plotmjj[c]->getObject(4),"Simulation","PE");
    legmc->AddEntry(plotmjj[c]->getObject(1),"Parametric Model","L");
    legmc->AddEntry(plotmjj[c]->getObject(2),"Gaussian ","L");
    legmc->AddEntry(plotmjj[c]->getObject(3),"Crystal Ball","L");
    legmc->SetHeader(" ");
    legmc->SetLineColor(kWhite);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();
    /////////////////////////////////////////////////////
  TPaveText *pt = new TPaveText(0.1,0.93,0.7,0.99, "brNDC");
   //pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetShadowColor(kWhite);
   pt->SetLineColor(kWhite);
   pt->AddText("CMS Preliminary Simulation ");
   pt->SetTextSize(0.035);
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetShadowColor(kWhite);
   pt->SetLineColor(kWhite);
   pt->Draw();
    ////////////////////////////////////////////////  
    // float effS = effSigma(hist);
    TLatex *lat = new TLatex(minSigFitMjj+10.0,0.7*plotmjj[c]->GetMaximum(),"m_{X} = SM GeV");
    lat->Draw();
    TLatex *lat2 = new TLatex(minSigFitMjj+10.0,0.8*plotmjj[c]->GetMaximum(),catdesc.at(c));
    lat2->Draw();
    TLatex *lat3 = new TLatex(minSigFitMjj+10.0,0.9*plotmjj[c]->GetMaximum(),"X#rightarrow HH #rightarrow #gamma#gammab#bar{b}");
    lat3->Draw();
    ///////
    char myChi2buffer[50];
    sprintf(myChi2buffer,"#chi^{2}/ndof = %f",chi2n);
    TLatex* latex = new TLatex(0.52, 0.7, myChi2buffer);
    latex -> SetNDC();
    latex -> SetTextFont(42);
    latex -> SetTextSize(0.04);
    //latex -> Draw("same");
    ctmp->SaveAs(TString::Format("sigmodelMjj_cat%d.pdf",c));
    ctmp->SaveAs(TString::Format("sigmodelMjj_cat%d.png",c));
    ctmp->SaveAs(TString::Format("sigmodelMjj_cat%d.root",c));
  } // close categories
  return;
} // close makeplots signal
////////////////////////////////////////////////////////////////////////
void MakePlotsHiggs(RooWorkspace* w, Float_t Mass) {
  const Int_t ncat = NCAT;
  std::vector<TString> catdesc;
  catdesc.push_back(" 2 btag");
  catdesc.push_back(" 1 btag");
  catdesc.push_back("cat 2");
  catdesc.push_back("cat 3");
  // retrieve data sets from the workspace
  // RooDataSet* dataAll = (RooDataSet*) w->data("Data");
  //RooDataSet* signalAll = (RooDataSet*) w->data("Sig");
  //RooDataSet* higgsAll = (RooDataSet*) w->data("Hig");
  // blinded dataset
  // RooDataSet* data[ncat];
    TString component[4] = {"ggH","ttH","VBF","VH"}; 
for (int d = 0; d < 4; ++d){
  RooDataSet* sigToFit[ncat];
  RooAbsPdf* mggGaussSig[ncat];
  RooAbsPdf* mggCBSig[ncat];
  RooAbsPdf* mggSig[ncat];
  //
  RooAbsPdf* mggBkg[ncat];
  for (int c = 0; c < ncat; ++c) {
  // data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    sigToFit[c] = (RooDataSet*) w->data(TString::Format("Hig_%d_cat%d",d,c));
    mggGaussSig[c] = (RooAbsPdf*) w->pdf(TString::Format("mggGaussHig_%d_cat%d",d,c));
    mggCBSig[c] = (RooAbsPdf*) w->pdf(TString::Format("mggCBHig_%d_cat%d",d,c));
    mggSig[c] = (RooAbsPdf*) w->pdf(TString::Format("mggHig_%d_cat%d",d,c));
    mggBkg[c] = (RooAbsPdf*) w->pdf(TString::Format("mggBkg_%d_cat%d",d,c));
  } // close categories
  RooRealVar* mgg = w->var("mgg");
  mgg->setUnit("GeV");
  //RooAbsPdf* mggBkgAll = w->pdf("mggBkg_cat1");
  //
  //****************************//
  // Plot mgg Fit results
  //****************************//
  // Set P.D.F. parameter names
  // WARNING: Do not use it if Workspaces are created
  // SetParamNames(w);
  Float_t minSigFit(120),maxSigFit(130);
  Float_t MASS(Mass);
  Int_t nBinsMass(20); // just need to plot
  //RooPlot* plotmggAll = mgg->frame(Range(minSigFit,maxSigFit),Bins(nBinsMass));
  //higgsAll->plotOn(plotmggAll);
  gStyle->SetOptTitle(0);
  TCanvas* c1 = new TCanvas("c1","mgg",0,0,500,500);
  c1->cd(1);
  //********************************************//
  // Plot Signal Categories
  //****************************//
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  RooPlot* plotmgg[ncat];
  for (int c = 0; c < ncat; ++c) {
    plotmgg[c] = mgg->frame(Range(minSigFit,maxSigFit),Bins(nBinsMass));
    plotmgg[c]->getAttMarker()->setMarkerType(25);
    sigToFit[c]->plotOn(plotmgg[c],LineColor(kWhite),MarkerColor(kWhite));
    mggSig[c] ->plotOn(plotmgg[c]);
    double chi2n = plotmgg[c]->chiSquare(0) ;
    cout << "------------------------- Experimentakl chi2 = " << chi2n << endl;
    mggSig[c] ->plotOn(
	plotmgg[c],
	Components(TString::Format("mggGaussHig_%d_cat%d",d,c)),
	LineStyle(kDashed),LineColor(kGreen));
    mggSig[c] ->plotOn(
	plotmgg[c],
	Components(TString::Format("mggCBHig_%d_cat%d",d,c)),
	LineStyle(kDashed),LineColor(kRed));
    //mggSig[c] ->paramOn(plotmgg[c]);
    sigToFit[c] ->plotOn(plotmgg[c]);
// TCanvas* dummy = new TCanvas("dummy", "dummy",0, 0, 400, 400);
    TH1F *hist = new TH1F("hist", "hist", 400, minSigFit, maxSigFit);
    plotmgg[c]->SetTitle("CMS preliminary 19.7/fb ");
    plotmgg[c]->SetMinimum(0.0);
    plotmgg[c]->SetMaximum(1.40*plotmgg[c]->GetMaximum());
    plotmgg[c]->GetXaxis()->SetTitle("M_{#gamma#gamma} (GeV)");
    TCanvas* ctmp = new TCanvas("ctmp","Background Categories",0,0,500,500);
    plotmgg[c]->Draw();
    plotmgg[c]->Draw("SAME");
    TLegend *legmc = new TLegend(0.62,0.75,0.99,0.99);

    legmc->AddEntry(plotmgg[c]->getObject(5),component[d],"LPE");
    legmc->AddEntry(plotmgg[c]->getObject(1),"Parametric Model","L");
    legmc->AddEntry(plotmgg[c]->getObject(2),"Gaussian","L");
    legmc->AddEntry(plotmgg[c]->getObject(3),"Crystal Ball ","L");
    legmc->SetHeader(" ");
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();
    // float effS = effSigma(hist);
    TLatex *lat = new TLatex(minSigFit+0.5,0.85*plotmgg[c]->GetMaximum()," M_{X} = SM GeV");
    lat->Draw();
    TLatex *lat2 = new TLatex(minSigFit+1.5,0.75*plotmgg[c]->GetMaximum(),catdesc.at(c));
    lat2->Draw();
    ///////
    char myChi2buffer[50];
    sprintf(myChi2buffer,"#chi^{2}/ndof = %f",chi2n);
    TLatex* latex = new TLatex(0.52, 0.7, myChi2buffer);
    latex -> SetNDC();
    latex -> SetTextFont(42);
    latex -> SetTextSize(0.04);
    //latex -> Draw("same");
    ctmp->SaveAs(TString::Format("higmodel_%d_cat%d.pdf",d,c));
    ctmp->SaveAs(TString::Format("higmodel_%d_cat%d.png",d,c));
    //ctmp->SaveAs(TString::Format("sigmodel_cat%d.C",c));
  } // close categories
 } // close to higgs component 
    return;
} // close makeplots signal
////////////////////////////////////////////////////////////////////
// we add the higgs to the workspace in categories
void AddHigData(RooWorkspace* w, Float_t mass, TString signalfile, int higgschannel) {
  const Int_t ncat = NCAT;
  Float_t MASS(mass);
  RooArgSet* ntplVars = defineVariables();
  TFile higFile(signalfile);
  TTree* higTree = (TTree*) higFile.Get("TCVARS");
  // common preselection cut
  TString mainCut("1");
  // one channel with right weights
  //RooDataSet higScaled1(
  RooDataSet higScaled(
  "higScaled",
  "dataset",
  higTree, // all variables of RooArgList
  *ntplVars,
  mainCut,
  "evWeight");
  /*
  RooRealVar *evWeight = (RooRealVar*) (*ntplVars)["evWeight"] ;
  RooRealVar *k = new RooRealVar("k", "k", 0.0006424);
  RooFormulaVar *nw = new RooFormulaVar("nw", "nw", "@1", RooArgSet(*evWeight, *k));
  higScaled1.addColumn(*nw);
  RooArgSet *ntplVars1 = higScaled1.get();
  RooDataSet *higScaled = new RooDataSet("higScaled", "dataset",higTree, *ntplVars1,"", "nw");
  */
  //
  RooDataSet* higToFit[ncat];
  TString cut0 = "&& 1>0";//"&& mtot > 955 && mtot < 1150 "; //
  TString cut1 = "&& 1>0";//"&& mtot > 955 && mtot < 1150 "; // "&& 1>0";//
  //
  TString cutj0 = "&& 1>0";//"&& mjj_wokinfit > 90 && mjj_wokinfit < 160 "; //"&& 1>0";//
  TString cutj1 = "&& 1>0";//"&& mjj_wokinfit > 100 && mjj_wokinfit < 170 "; // "&& 1>0";//
  //
  // we take only mtot to fit to the workspace, we include the cuts
  higToFit[0] = (RooDataSet*) higScaled.reduce(
    RooArgList(*w->var("mgg"),*w->var("mjj")),
     mainCut+TString::Format(" && cut_based_ct==%d ",0)+cut0+cutj0);
  w->import(*higToFit[0],Rename(TString::Format("Hig_%d_cat%d",higgschannel,0)));
  higToFit[1] = (RooDataSet*) higScaled.reduce(
     RooArgList(*w->var("mgg"),*w->var("mjj")),
     mainCut+TString::Format(" && cut_based_ct==%d ",1)+cut1+cutj1);
  w->import(*higToFit[1],Rename(TString::Format("Hig_%d_cat%d",higgschannel,1))); // Create full signal data set without categorization
  // here we print the number of entries on the different categories
  cout << "========= the number of entries on the different categories ==========" << endl;
  cout << "---- one channel: " << higScaled.sumEntries() << endl;
  for (int c = 0; c < ncat; ++c) {
    Float_t nExpEvt = higToFit[c]->sumEntries();
    cout << TString::Format("nEvt exp. cat%d : ",c) << nExpEvt
<< TString::Format(" eff x Acc cat%d : ",c)
<< "%"
<< endl;
  }
  cout << "======================================================================" << endl;
  higScaled.Print("v");
  return;
} // end add higgs function
///////////////////////////////////////////////////////////////
void MakeHigWS(RooWorkspace* w, const char* fileHiggsName,int higgschannel) {
  TString wsDir = "workspaces/";
  const Int_t ncat = NCAT;
  //**********************************************************************//
  // Write pdfs and datasets into the workspace before to save to a file
  // for statistical tests.
  //**********************************************************************//
  RooAbsPdf* mggHigPdf[ncat];
  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");
  for (int c = 0; c < ncat; ++c) {
    mggHigPdf[c] = (RooAbsPdf*) w->pdf(TString::Format("mggHig_%d_cat%d",higgschannel,c));
    wAll->import(*w->pdf(TString::Format("mggHig_%d_cat%d",higgschannel,c)));
  }
  // (2) Systematics on energy scale and resolution
  // 1,1,1 statistical to be treated on the datacard
  wAll->factory("CMS_hgg_sig_m0_absShift[1,1,1]");//TString::Format("CMS_hgg_hig_%d_m0_absShift[1,1,1]",higgschannel));
  wAll->factory(TString::Format("prod::CMS_hgg_hig_m0_%d_cat0(mgg_hig_m0_%d_cat0, CMS_hgg_sig_m0_absShift)",higgschannel,higgschannel));
  wAll->factory(TString::Format("prod::CMS_hgg_hig_m0_%d_cat1(mgg_hig_m0_%d_cat1, CMS_hgg_sig_m0_absShift)",higgschannel,higgschannel));
  // (3) Systematics on resolution
  wAll->factory("CMS_hgg_sig_sigmaScale[1,1,1]");//TString::Format("CMS_hgg_hig_%d_sigmaScale[1,1,1]",higgschannel));
  wAll->factory(TString::Format("prod::CMS_hgg_hig_sigma_%d_cat0(mgg_hig_sigma_%d_cat0, CMS_hgg_sig_sigmaScale)",higgschannel,higgschannel));

  wAll->factory(TString::Format("prod::CMS_hgg_hig_sigma_%d_cat1(mgg_hig_sigma_%d_cat1, CMS_hgg_sig_sigmaScale)",higgschannel,higgschannel));
  wAll->factory(TString::Format("prod::CMS_hgg_hig_gsigma_%d_cat0(mgg_hig_gsigma_%d_cat0, CMS_hgg_sig_sigmaScale)",higgschannel,higgschannel));
  wAll->factory(TString::Format("prod::CMS_hgg_hig_gsigma_%d_cat1(mgg_hig_gsigma_%d_cat1, CMS_hgg_sig_sigmaScale)",higgschannel,higgschannel));
  // save the other parameters
  /* for (int c = 0; c < ncat; ++c) {
wAll->factory(
TString::Format("CMS_hgg_sig_alpha_cat%d[%g,0.5,5]",
c, wAll->var(TString::Format("mgg_sig_alpha_cat%d",c))->getVal()));
wAll->factory(
TString::Format("CMS_hgg_sig_n_cat%d[%g,0.5,20]",
c, wAll->var(TString::Format("mgg_sig_n_cat%d",c))->getVal()));
wAll->factory(
TString::Format("CMS_hgg_sig_frac_cat%d[%g,0.0,1.0]",
c, wAll->var(TString::Format("mgg_sig_frac_cat%d",c))->getVal()));
}
*/
  // (4) do reparametrization of signal
  for (int c = 0; c < ncat; ++c) wAll->factory(
TString::Format("EDIT::CMS_hgg_hig_%d_cat%d(mggHig_%d_cat%d,",higgschannel,c,higgschannel,c) +
TString::Format(" mgg_hig_m0_%d_cat%d=CMS_hgg_hig_m0_%d_cat%d, ",higgschannel, c,higgschannel,c) +
TString::Format(" mgg_hig_sigma_%d_cat%d=CMS_hgg_hig_sigma_%d_cat%d, ",higgschannel, c,higgschannel,c) +
TString::Format(" mgg_hig_gsigma_%d_cat%d=CMS_hgg_hig_gsigma_%d_cat%d)",higgschannel, c,higgschannel,c)
  );
  TString filename(wsDir+TString(fileHiggsName)+".inputsig.root");
  wAll->writeToFile(filename);
  cout << "Write signal workspace in: " << filename << " file" << endl;
  return;
} // close make higgs WP
///////////////////////////////////////////////////////////
// declare histos or what -> NOT USED
Double_t effSigma(TH1 *hist) {
  TAxis *xaxis = hist->GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    std::cout << "effsigma: Not a valid histo. nbins = " << nb << std::endl;
    return 0.;
  }
  Double_t bwid = xaxis->GetBinWidth(1);
  if(bwid == 0) {
    std::cout << "effsigma: Not a valid histo. bwid = " << bwid << std::endl;
    return 0.;
  }
  Double_t xmax = xaxis->GetXmax();
  Double_t xmin = xaxis->GetXmin();
  Double_t ave = hist->GetMean();
  Double_t rms = hist->GetRMS();
  Double_t total=0.;
  for(Int_t i=0; i<nb+2; i++) {
    total+=hist->GetBinContent(i);
  }
  if(total < 100.) {
    std::cout << "effsigma: Too few entries " << total << std::endl;
    return 0.;
  }
  Int_t ierr=0;
  Int_t ismin=999;
  Double_t rlim=0.683*total;
  Int_t nrms=rms/(bwid); // Set scan size to +/- rms
  if(nrms > nb/10) nrms=nb/10; // Could be tuned...
  Double_t widmin=9999999.;
  for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
    Int_t ibm=(ave-xmin)/bwid+1+iscan;
    Double_t x=(ibm-0.5)*bwid+xmin;
    Double_t xj=x;
    Double_t xk=x;
    Int_t jbm=ibm;
    Int_t kbm=ibm;
    Double_t bin=hist->GetBinContent(ibm);
    total=bin;
    for(Int_t j=1;j<nb;j++){
      if(jbm < nb) {jbm++; xj+=bwid; bin=hist->GetBinContent(jbm); total+=bin; if(total > rlim) break;} else ierr=1;
      if(kbm > 0) {kbm--; xk-=bwid; bin=hist->GetBinContent(kbm); total+=bin; if(total > rlim) break; } else ierr=1;
    }
    Double_t dxf=(total-rlim)*bwid/bin;
    Double_t wid=(xj-xk+bwid-dxf)*0.5;
    if(wid < widmin) { widmin=wid; ismin=iscan; }
  } // Scan window centre
  if(ismin == nrms || ismin == -nrms) ierr=3;
  if(ierr != 0) std::cout << "effsigma: Error of type " << ierr << std::endl;
  return widmin;
} // close effSigma
//////////////////////////////////////////////////
// with higgs
void MakeDataCard(RooWorkspace* w, const char* fileBaseName, const char* fileBkgName , const char* fileHiggsNameggh, const char* fileHiggsNametth, const char* fileHiggsNamevbf, const char* fileHiggsNamevh) {
  TString cardDir = "datacards/";
  const Int_t ncat = NCAT;
  RooDataSet* data[ncat];
  RooDataSet* sigToFit[ncat];
  RooDataSet* higToFitggh[ncat];
  RooDataSet* higToFittth[ncat];
  RooDataSet* higToFitvbf[ncat];
  RooDataSet* higToFitvh[ncat];
  for (int c = 0; c < ncat; ++c) {
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    sigToFit[c] = (RooDataSet*) w->data(TString::Format("Sig_cat%d",c));
    //
    higToFitggh[c] = (RooDataSet*) w->data(TString::Format("Hig_%d_cat%d",0,c));
    higToFittth[c] = (RooDataSet*) w->data(TString::Format("Hig_%d_cat%d",1,c));
    higToFitvbf[c] = (RooDataSet*) w->data(TString::Format("Hig_%d_cat%d",2,c));
    higToFitvh[c] = (RooDataSet*) w->data(TString::Format("Hig_%d_cat%d",3,c));
  } // close cat
  ////////////////////////////////////////////////////////////////////////////////////
  //RooRealVar* lumi = w->var("lumi");
  cout << "======== Expected Events Number =====================" << endl;
  cout << ".........Measured Data for L = " << "19785" << " pb-1 ............................" << endl;
  cout << "#Events data: " << w->data("Data")->sumEntries() << endl;
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events data cat%d: ",c) << data[c]->sumEntries() << endl;
  }
  cout << ".........Expected Signal for L = " << "19785" << " pb-1 ............................" << endl;
  cout << "#Events Signal: " << w->data("Data")->sumEntries() << endl;
  Float_t siglikeErr[ncat];
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events Signal cat%d: ",c) << sigToFit[c]->sumEntries() << endl;
    siglikeErr[c]=0.6*sigToFit[c]->sumEntries();
  }
  cout << "====================================================" << endl;
  TString filename(cardDir+TString(fileBaseName)+".txt");
  ofstream outFile(filename);

// outFile << "#CMS-HGG DataCard for Unbinned Limit Setting, " << lumi->getVal() << " pb-1 " << endl;
  outFile << "#Run with: combine -d hgg.mH260.0.shapes-Unbinned.txt -U -m 130 -H ProfileLikelihood -M MarkovChainMC --rMin=0 --rMax=20.0 -b 3000 -i 50000 --optimizeSim=1 --tries 30" << endl;
  outFile << "# Lumi = " << "19785" << " pb-1" << endl;
  outFile << "imax "<<ncat << endl;
  outFile << "jmax 5" << endl; // number of BKG
  outFile << "kmax *" << endl;
  outFile << "---------------" << endl;
  outFile << "shapes data_obs cat0 " << TString(fileBkgName)+".root" << " w_all:data_obs_cat0" << endl;
  outFile << "shapes data_obs cat1 "<< TString(fileBkgName)+".root" << " w_all:data_obs_cat1" << endl;
  outFile << "############## shape with reparametrization" << endl;
  outFile << "shapes mggBkg cat0 " << TString(fileBkgName)+".root" << " w_all:CMS_hgg_bkg_8TeV_cat0" << endl;
  outFile << "shapes mggBkg cat1 "<< TString(fileBkgName)+".root" << " w_all:CMS_hgg_bkg_8TeV_cat1" << endl;
  outFile << "# signal" << endl;
  outFile << "shapes mggSig cat0 " << TString(fileBaseName)+".inputsig.root" << " w_all:CMS_hgg_sig_cat0" << endl;
  outFile << "shapes mggSig cat1 " << TString(fileBaseName)+".inputsig.root" << " w_all:CMS_hgg_sig_cat1" << endl;
  outFile << "# ggh" << endl;
  outFile << "shapes mggHigggh cat0 " << TString(fileHiggsNameggh)+".inputsig.root" << " w_all:CMS_hgg_hig_0_cat0" << endl;
  outFile << "shapes mggHigggh cat1 " << TString(fileHiggsNameggh)+".inputsig.root" << " w_all:CMS_hgg_hig_0_cat1" << endl;
  outFile << "# tth" << endl;
  outFile << "shapes mggHigtth cat0 " << TString(fileHiggsNametth)+".inputsig.root" << " w_all:CMS_hgg_hig_1_cat0" << endl;
  outFile << "shapes mggHigtth cat1 " << TString(fileHiggsNametth)+".inputsig.root" << " w_all:CMS_hgg_hig_1_cat1" << endl;
  outFile << "# vbf" << endl;
  outFile << "shapes mggHigvbf cat0 " << TString(fileHiggsNamevbf)+".inputsig.root" << " w_all:CMS_hgg_hig_2_cat0" << endl;
  outFile << "shapes mggHigvbf cat1 " << TString(fileHiggsNamevbf)+".inputsig.root" << " w_all:CMS_hgg_hig_2_cat1" << endl;
  outFile << "# vh" << endl;
  outFile << "shapes mggHigvh cat0 " << TString(fileHiggsNamevh)+".inputsig.root" << " w_all:CMS_hgg_hig_3_cat0" << endl;
  outFile << "shapes mggHigvh cat1 " << TString(fileHiggsNamevh)+".inputsig.root" << " w_all:CMS_hgg_hig_3_cat1" << endl;
/* CMS_hgg_hig_0_
  outFile << "shapes mggHigggh cat0 " << TString(fileHiggsNameggh)+".inputsig.root" << " w_all:CMS_hgg_sig_ggh_cat0" << endl;
  outFile << "shapes mggHigggh cat1 " << TString(fileHiggsNameggh)+".inputsig.root" << " w_all:CMS_hgg_sig_ggh_cat1" << endl;
  outFile << "# tth" << endl;
  outFile << "shapes mggHigtth cat0 " << TString(fileHiggsNametth)+".inputsig.root" << " w_all:CMS_hgg_sig_tth_cat0" << endl;
  outFile << "shapes mggHigtth cat1 " << TString(fileHiggsNametth)+".inputsig.root" << " w_all:CMS_hgg_sig_tth_cat1" << endl;
  outFile << "# vbf" << endl;
  outFile << "shapes mggHigvbf cat0 " << TString(fileHiggsNamevbf)+".inputsig.root" << " w_all:CMS_hgg_sig_vbf_cat0" << endl;
  outFile << "shapes mggHigvbf cat1 " << TString(fileHiggsNamevbf)+".inputsig.root" << " w_all:CMS_hgg_sig_vbf_cat1" << endl;
  outFile << "# vh" << endl;
  outFile << "shapes mggHigvh cat0 " << TString(fileHiggsNamevh)+".inputsig.root" << " w_all:CMS_hgg_sig_vh_cat0" << endl;
  outFile << "shapes mggHigvh cat1 " << TString(fileHiggsNamevh)+".inputsig.root" << " w_all:CMS_hgg_sig_vh_cat1" << endl;
*/
  outFile << "---------------" << endl;
  /////////////////////////////////////
  if(addHiggs) { //
  outFile << "bin cat0 cat1 " << endl;
  cout<<"here"<<endl;
  outFile << "observation "<< data[0]->sumEntries() <<" " << data[1]->sumEntries() <<" "<< endl;
  outFile << "------------------------------" << endl;
  outFile << "bin cat0 cat0 cat0 cat0 cat0 cat0"
                                <<" cat1 cat1 cat1 cat1 cat1 cat1" << endl;
  outFile << "process mggSig mggBkg mggHigggh mggHigtth mggHigvbf mggHigvh"
                            <<" mggSig mggBkg mggHigggh mggHigtth mggHigvbf mggHigvh" << endl;
  outFile << "process 0 1 2 3 4 5 "
                                 <<" 0 1 2 3 4 5 " << endl;
  outFile << "rate "
<<" "<<sigToFit[0]->sumEntries()<<" "<<1<<" "<<higToFitggh[0]->sumEntries()<<" "<<higToFittth[0]->sumEntries()<<" "<<higToFitvbf[0]->sumEntries()<<" "<<higToFitvh[0]->sumEntries()
<<" "<<sigToFit[1]->sumEntries()<<" "<<1<<" "<<higToFitggh[1]->sumEntries()<<" "<<higToFittth[1]->sumEntries()<<" "<<higToFitvbf[1]->sumEntries()<<" "<<higToFitvh[1]->sumEntries()
<<" "<<endl;
  outFile << "     " << endl;
  outFile << "############## Total normalisation" << endl;
  outFile << "lumi_8TeV lnN "
<< "1.026 - 1.026 1.026 1.026 1.026 "
<< "1.026 - 1.026 1.026 1.026 1.026 " << endl;
  outFile << "     " << endl;
  outFile << "############## Photon selection normalisation uncertainties " << endl;
  outFile << "DiphoTrigger lnN "
<< "1.01 - 1.010 1.010 1.010 1.010 "
<< "1.01 - 1.010 1.010 1.010 1.010 "
<< "# Trigger efficiency" << endl;
  outFile << "CMS_hgg_eff_g lnN "
   << "1.010 - 1.010 1.010 1.010 1.010 "
   << "1.010 - 1.010 1.010 1.010 1.010 "
   << "# photon selection accep." << endl;
  outFile << "     " << endl;
  outFile << "############## Jet selection and phase space cuts normalisation uncertainties " << endl;
  outFile << "Mjj_PTj_cut_acceptance lnN "
<< "1.015 - 1.015 1.015 1.015 1.015 "
<< "1.015 - 1.015 1.015 1.015 1.015 "
<<"# JER and JES " << endl;
  outFile << "btag_eff lnN "
<< "1.046 - 1.046 1.046 1.046 1.046 "
<< "0.988 - 0.988 0.988 0.988 0.988 "
<<"# b tag efficiency uncertainty" << endl;
  outFile << "maajj_cut_acceptance lnN "
   << "1.02 - 1.02 1.02 1.02 1.02 "
   << "1.02 - 1.02 1.02 1.02 1.02 " << endl;
  outFile << "     " << endl;
  outFile << "############## Theory uncertainties on SM Higgs production " << endl;
  outFile << "PDF lnN "
   << " - - 0.931/1.075 0.919/1.081 0.972/1.026 0.976/1.024 "
   << " - - 0.931/1.075 0.919/1.081 0.972/1.026 0.976/1.024 " << endl;
  outFile << "QCD_scale lnN "
   << " - - 0.922/1.072 0.907/1.038 0.998/1.002 0.980/1.020 "
   << " - - 0.922/1.072 0.907/1.038 0.998/1.002 0.980/1.020 " << endl;
  outFile << "gg_migration lnN "
   << " - - 1.25 1.25 1.08 1.08 "
   << " - - 1.25 1.25 1.08 1.08 # UEPS" << endl;
  outFile << "gluonSplitting lnN "
   << " - - 1.40 1.40 1.40 1.40 "
   << " - - 1.40 1.40 1.40 1.40 " << endl;
  outFile << "     " << endl;
  outFile << "############## Signal parametric shape uncertainties " << endl;
  outFile << "CMS_hgg_sig_m0_absShift param 1 0.0057 # displacement of the dipho mean error = sqrt(0.45^ 2 + 0.35^ 2) " << endl;
  outFile << "CMS_hgg_sig_sigmaScale param 1 0.22 # optimistic estimative of resolution uncertainty " << endl;
  //
/*  outFile << "# Parametric shape uncertainties, entered by hand. they act on higgs" << endl;
  outFile << "CMS_hgg_hig_m0_0_absShift param 1 0.0057 # displacement of the dipho mean error = sqrt(0.45^ 2 + 0.35^ 2)" << endl;
  outFile << "CMS_hgg_hig_0_sigmaScale param 1 0.22 # optimistic estimative of resolution uncertainty " << endl;
  //
  outFile << "CMS_hgg_hig_m0_1_absShift param 1 0.0057 # displacement of the dipho mean error = sqrt(0.45^ 2 + 0.35^ 2)" << endl;
  outFile << "CMS_hgg_hig_1_sigmaScale param 1 0.22 # optimistic estimative of resolution uncertainty " << endl;
  //
  outFile << "CMS_hgg_hig_m0_2_absShift param 1 0.0057 # displacement of the dipho mean error = sqrt(0.45^ 2 + 0.35^ 2)" << endl;
  outFile << "CMS_hgg_hig_2_sigmaScale param 1 0.22 # optimistic estimative of resolution uncertainty " << endl;
  //
  outFile << "CMS_hgg_hig_m0_3_absShift param 1 0.0057 # displacement of the dipho mean error = sqrt(0.45^ 2 + 0.35^ 2)" << endl;
  outFile << "CMS_hgg_hig_3_sigmaScale param 1 0.22 # optimistic estimative of resolution uncertainty " << endl;
  //
*/
  outFile << "############## for mgg fit - slopes" << endl;
  outFile << "CMS_hgg_bkg_8TeV_cat0_norm flatParam # Normalization uncertainty on background slope" << endl;
  outFile << "CMS_hgg_bkg_8TeV_cat1_norm flatParam # Normalization uncertainty on background slope" << endl;

  outFile << "CMS_hgg_bkg_8TeV_slope1_cat0 flatParam # Mean and absolute uncertainty on background slope" << endl;
  outFile << "CMS_hgg_bkg_8TeV_slope1_cat1 flatParam # Mean and absolute uncertainty on background slope" << endl;

  } // if ncat ==2
  /////////////////////////////////////
  outFile.close();
  cout << "Write data card in: " << filename << " file" << endl;
  return;
} // close write full datacard
//////////////////////////////////////////////////
///// datacards witout higgs
//////////////////////////////////////////////////
// with reparametrization of BKG
void MakeDataCardREP(RooWorkspace* w, const char* fileBaseName, const char* fileBkgName) {
  TString cardDir = "datacards/";
  const Int_t ncat = NCAT;
  RooDataSet* data[ncat];
  RooDataSet* sigToFit[ncat];
  for (int c = 0; c < ncat; ++c) {
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    sigToFit[c] = (RooDataSet*) w->data(TString::Format("Sig_cat%d",c));
  }
  //RooRealVar* lumi = w->var("lumi");
  cout << "======== Expected Events Number =====================" << endl;
  //cout << ".........Measured Data for L = " << lumi->getVal() << " pb-1 ............................" << endl;
  cout << "#Events data: " << w->data("Data")->sumEntries() << endl;
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events data cat%d: ",c) << data[c]->sumEntries() << endl;
  }
  //cout << ".........Expected Signal for L = " << lumi->getVal() << " pb-1 ............................" << endl;
  cout << "#Events Signal: " << w->data("Data")->sumEntries() << endl;
  Float_t siglikeErr[ncat];
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events Signal cat%d: ",c) << sigToFit[c]->sumEntries() << endl;
    siglikeErr[c]=0.6*sigToFit[c]->sumEntries();
  }
  cout << "====================================================" << endl;
  TString filename(cardDir+TString(fileBaseName)+"rep.txt");
  ofstream outFile(filename);
  //outFile << "#CMS-HGG DataCard for Unbinned Limit Setting, " << lumi->getVal() << " pb-1 " << endl;
  outFile << "#Run with: combine -d hgg.mH260.0.shapes-Unbinned.txt -U -m 130 -H ProfileLikelihood -M MarkovChainMC --rMin=0 --rMax=20.0 -b 3000 -i 50000 --optimizeSim=1 --tries 30" << endl;
  //outFile << "# Lumi = " << lumi->getVal() << " pb-1" << endl;
  outFile << "imax "<<ncat << endl;
  outFile << "jmax 1" << endl;
  outFile << "kmax *" << endl;
  outFile << "---------------" << endl;
cout<<"here"<<endl;
  outFile << "# the name after w_all is the name of the rooextpdf we want to use, we have both saved" << endl;
  outFile << "# BKG" << endl;
  outFile << "shapes data_obs cat0 " << TString(fileBkgName)+".root" << " w_all:data_obs_cat0" << endl;
  outFile << "shapes data_obs cat1 "<< TString(fileBkgName)+".root" << " w_all:data_obs_cat1" << endl;
  outFile << "############## shape with reparametrization" << endl;
  outFile << "shapes mtotBkg cat0 " << TString(fileBkgName)+".root" << " w_all:CMS_bkg_8TeV_cat0" << endl;
  outFile << "shapes mtotBkg cat1 "<< TString(fileBkgName)+".root" << " w_all:CMS_bkg_8TeV_cat1" << endl;
  outFile << "# signal" << endl;
  outFile << "shapes mtotSig cat0 " << TString(fileBaseName)+".inputsig.root" << " w_all:CMS_sig_cat0" << endl;
  outFile << "shapes mtotSig cat1 " << TString(fileBaseName)+".inputsig.root" << " w_all:CMS_sig_cat1" << endl;
  outFile << "---------------" << endl;
  /////////////////////////////////////
  outFile << "bin cat0 cat1 " << endl;
  outFile << "observation "
<< data[0]->sumEntries() << " "
<< data[1]->sumEntries() << " "
<< endl;
  outFile << "------------------------------" << endl;
  outFile << "bin cat0 cat0 cat1 cat1" << endl;
  outFile << "process mtotSig mtotBkg mtotSig mtotBkg" << endl;
  outFile << "process 0 1 0 1" << endl;
  outFile << "rate "
<< " " << sigToFit[0]->sumEntries() << " " << 1
<< " " << sigToFit[1]->sumEntries() << " " << 1
<< " " << endl;
  outFile << "--------------------------------" << endl;
  outFile << "lumi_8TeV lnN "
<< "1.026 - "
<< "1.026 - " << endl;
  outFile << "DiphoTrigger lnN "
<< "1.01 - "
<< "1.01 - "
<< "# Trigger efficiency" << endl;
  outFile << "############## jet" << endl;
  outFile << "Mjj_acceptance lnN "
<< "1.015 - "
<< "1.015 - "
<<"# JER and JES " << endl;
  outFile << "btag_eff lnN "
<< "1.046 - "
<< "0.988 - "
<<"# b tag efficiency uncertainty" << endl;
  outFile << "############## photon " << endl;
  outFile << "Photon_selectios_accep lnN "
<< "1.01 - "
<< "1.01 - "
<<"# photon acceptance" << endl;
  outFile << "CMS_hgg_eff_g lnN "
   << "1.010 - "
   << "1.010 - "
   << "# photon selection accep." << endl;
  outFile << "############## for mtot fit" << endl;
  outFile << "maajj_acceptance lnN "
   << "1.02 - "
   << "1.02 - "
   << "# photon energy resolution" << endl;
  outFile << "############## normalization floating" << endl;
  //outFile << "# Parametric shape uncertainties, entered by hand. they act on both higgs/signal " << endl;
  //
  outFile << "# Parametric shape uncertainties, entered by hand." << endl;
  outFile << "CMS_hgg_sig_m0_absShift param 1 0.0045 # displacement of the dipho mean" << endl;
  outFile << "CMS_hgg_sig_sigmaScale param 1 0.22 # optimistic estimative of resolution uncertainty " << endl;
  //
  //outFile << "CMS_hgg_sig_m0_absShift param 1 0.006 # displacement of the dipho mean" << endl;
  //outFile << "CMS_hgg_sig_sigmaScale param 1 0.30 # optimistic estimative of resolution uncertainty " << endl;
  outFile << "############## for mtot fit - slopes" << endl;
  outFile << "############## with reparametrization" << endl;
  outFile << "CMS_bkg_8TeV_cat0_norm flatParam # Normalization uncertainty on background slope" << endl;
  outFile << "CMS_bkg_8TeV_cat1_norm flatParam # Normalization uncertainty on background slope" << endl;
  outFile << "CMS_hgg_bkg_8TeV_slope1_cat0 flatParam # Mean and absolute uncertainty on background slope" << endl;
  outFile << "CMS_hgg_bkg_8TeV_slope1_cat1 flatParam # Mean and absolute uncertainty on background slope" << endl;
  outFile << "CMS_hbb_bkg_8TeV_slope1_cat0 flatParam # Mean and absolute uncertainty on background slope" << endl;
  outFile << "CMS_hbb_bkg_8TeV_slope1_cat1 flatParam # Mean and absolute uncertainty on background slope" << endl;
  //outFile << "CMS_hbb_bkg_8TeV_slope2_cat0 flatParam # Mean and absolute uncertainty on background slope" << endl;
  //outFile << "CMS_hbb_bkg_8TeV_slope2_cat1 flatParam # Mean and absolute uncertainty on background slope" << endl;
  //outFile << "CMS_hbb_bkg_8TeV_slope3_cat0 flatParam # Mean and absolute uncertainty on background slope" << endl;
  //outFile << "CMS_hbb_bkg_8TeV_slope3_cat1 flatParam # Mean and absolute uncertainty on background slope" << endl;
  //outFile << "CMS_hbb_bkg_8TeV_slope4_cat0 flatParam # Mean and absolute uncertainty on background slope" << endl;
  //outFile << "CMS_hbb_bkg_8TeV_slope4_cat1 flatParam # Mean and absolute uncertainty on background slope" << endl;
  //outFile << "CMS_hbb_bkg_8TeV_slope5_cat0 flatParam # Mean and absolute uncertainty on background slope" << endl;
  //outFile << "CMS_hbb_bkg_8TeV_slope5_cat1 flatParam # Mean and absolute uncertainty on background slope" << endl;
  outFile.close();
  cout << "Write data card in: " << filename << " file" << endl;
  return;
} // close write datacard with rep
//////////////////////////////////////////////////
// withou reparametrization of BKG
void MakeDataCardLNU(RooWorkspace* w, const char* fileBaseName, const char* fileBkgName) {
  TString cardDir = "datacards/";
  const Int_t ncat = NCAT;
  RooDataSet* data[ncat];
  RooDataSet* sigToFit[ncat];
  for (int c = 0; c < ncat; ++c) {
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    sigToFit[c] = (RooDataSet*) w->data(TString::Format("Sig_cat%d",c));
  }
  //RooRealVar* lumi = w->var("lumi");
  cout << "======== Expected Events Number =====================" << endl;
  //cout << ".........Measured Data for L = " << lumi->getVal() << " pb-1 ............................" << endl;
  cout << "#Events data: " << w->data("Data")->sumEntries() << endl;
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events data cat%d: ",c) << data[c]->sumEntries() << endl;
  }
  //cout << ".........Expected Signal for L = " << lumi->getVal() << " pb-1 ............................" << endl;
  cout << "#Events Signal: " << w->data("Data")->sumEntries() << endl;
  Float_t siglikeErr[ncat];
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events Signal cat%d: ",c) << sigToFit[c]->sumEntries() << endl;
    siglikeErr[c]=0.6*sigToFit[c]->sumEntries();
  }
  cout << "====================================================" << endl;
  TString filename(cardDir+TString(fileBaseName)+"lnu.txt");
  ofstream outFile(filename);
  //outFile << "#CMS-HGG DataCard for Unbinned Limit Setting, " << lumi->getVal() << " pb-1 " << endl;
  outFile << "#Run with: combine -d hgg.mH260.0.shapes-Unbinned.txt -U -m 130 -H ProfileLikelihood -M MarkovChainMC --rMin=0 --rMax=20.0 -b 3000 -i 50000 --optimizeSim=1 --tries 30" << endl;
  //outFile << "# Lumi = " << lumi->getVal() << " pb-1" << endl;
  outFile << "imax "<<ncat << endl;
  outFile << "jmax 1" << endl;
  outFile << "kmax *" << endl;
  outFile << "---------------" << endl;
cout<<"here"<<endl;
  outFile << "# the name after w_all is the name of the rooextpdf we want to use, we have both saved" << endl;
  outFile << "# BKG" << endl;
  outFile << "shapes data_obs cat0 " << TString(fileBkgName)+".root" << " w_all:data_obs_cat0" << endl;
  outFile << "shapes data_obs cat1 "<< TString(fileBkgName)+".root" << " w_all:data_obs_cat1" << endl;
  outFile << "############## shape with reparametrization" << endl;
  outFile << "shapes mtotBkg cat0 " << TString(fileBkgName)+".root" << " w_all:mggBkg_cat0" << endl;
  outFile << "shapes mtotBkg cat1 "<< TString(fileBkgName)+".root" << " w_all:mggBkg_cat1" << endl;
  outFile << "# signal" << endl;
  outFile << "shapes mtotSig cat0 " << TString(fileBaseName)+".inputsig.root" << " w_all:CMS_hgg_sig_cat0" << endl;
  outFile << "shapes mtotSig cat1 " << TString(fileBaseName)+".inputsig.root" << " w_all:CMS_hgg_sig_cat1" << endl;
  outFile << "---------------" << endl;
  /////////////////////////////////////
  /////////////////////////////////////
  outFile << "bin cat0 cat1 " << endl;
  outFile << "observation "
<< data[0]->sumEntries() << " "
<< data[1]->sumEntries() << " "
<< endl;
  outFile << "------------------------------" << endl;
  outFile << "bin cat0 cat0 cat1 cat1" << endl;
  outFile << "process mtotSig mtotBkg mtotSig mtotBkg" << endl;
  outFile << "process 0 1 0 1" << endl;
  outFile << "rate "
<< " " << sigToFit[0]->sumEntries() << " " << data[0]->sumEntries()
<< " " << sigToFit[1]->sumEntries() << " " << data[1]->sumEntries()
<< " " << endl;
  outFile << "--------------------------------" << endl;
  outFile << "lumi_8TeV lnN "
<< "1.022 - "
<< "1.022 - " << endl;
  outFile << "############## jet" << endl;
  outFile << "Mjj_acceptance lnN "
<< "1.015 - "
<< "1.015 - "
<<"# JER and JES " << endl;
  outFile << "Photon_selectios_accep lnN "
<< "1.08 - "
<< "1.08 - "
<<"# photon acceptance" << endl;
  outFile << "btag_eff lnN "
<< "1.06 - "
<< "1.03 - "
<<"# b tag efficiency uncertainty" << endl;
  outFile << "############## photon " << endl;
  outFile << "CMS_hgg_eff_g lnN "
   << "1.010 - "
   << "1.010 - "
   << "# photon selection accep." << endl;
  outFile << "DiphoTrigger lnN "
<< "1.01 - "
<< "1.01 - "
<< "# Trigger efficiency" << endl;
  outFile << "############## for mtot fit" << endl;
  outFile << "maa_acceptance lnN "
   << "1.10 - "
   << "1.10 - "
   << "# photon energy resolution" << endl;
  outFile << "mggBkg lnU "
   << " - 2 "
   << " - 2 "
   << "# photon energy resolution" << endl;
  outFile << "############## normalization floating" << endl;
  outFile << "# Parametric shape uncertainties, entered by hand. they act on both higgs/signal " << endl;
  outFile << "CMS_hgg_sig_m0_absShift param 1 0.0045 # displacement of the dipho mean 0.0075 = sqrt(0.006**2 + 0.004**2)" << endl;
  outFile << "CMS_hgg_sig_sigmaScale param 1 0.22 # optimistic estimative of resolution uncertainty " << endl;
  outFile << "############## for mtot fit - slopes" << endl;
  outFile << "############## with reparametrization" << endl;
  outFile << "mgg_bkg_8TeV_slope1_cat0 flatParam # Mean and absolute uncertainty on background slope" << endl;
  outFile << "mgg_bkg_8TeV_slope1_cat1 flatParam # Mean and absolute uncertainty on background slope" << endl;
//  outFile << "#mgg_bkg_8TeV_slope2_cat0 flatParam # Mean and absolute uncertainty on background slope" << endl;
//  outFile << "mgg_bkg_8TeV_slope2_cat1 flatParam # Mean and absolute uncertainty on background slope" << endl;
//  outFile << "#mgg_bkg_8TeV_slope3_cat0 flatParam # Mean and absolute uncertainty on background slope" << endl;
//  outFile << "mgg_bkg_8TeV_slope3_cat1 flatParam # Mean and absolute uncertainty on background slope" << endl;
//  outFile << "#mgg_bkg_8TeV_slope4_cat1 flatParam # Mean and absolute uncertainty on background slope" << endl;
//($\sim$ 125 GeV)  outFile << "#mgg_bkg_8TeV_slope5_cat1 flatParam # Mean and absolute uncertainty on background slope" << endl;
  outFile.close();
  cout << "Write data card in: " << filename << " file" << endl;
  return;
} // close write datacard without rep
////////////////////////////////////////////////////////////////////////////////
void MakeDataCardonecat(RooWorkspace* w, const char* fileBaseName, const char* fileBkgName , 
const char* fileHiggsNameggh, const char* fileHiggsNametth, const char* fileHiggsNamevbf, const char* fileHiggsNamevh) {
  TString cardDir = "datacards/";
  const Int_t ncat = NCAT;
  RooDataSet* data[ncat];
  RooDataSet* sigToFit[ncat];
  //RooDataSet* higToFit[ncat];
  RooDataSet* higToFitggh[ncat];
  RooDataSet* higToFittth[ncat];
  RooDataSet* higToFitvbf[ncat];
  RooDataSet* higToFitvh[ncat];
  for (int c = 0; c < ncat; ++c) {
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    sigToFit[c] = (RooDataSet*) w->data(TString::Format("Sig_cat%d",c));
    //higToFit[c] = (RooDataSet*) w->data(TString::Format("Hig_cat%d",c));
    higToFitggh[c] = (RooDataSet*) w->data(TString::Format("Hig_%d_cat%d",0,c));
    higToFittth[c] = (RooDataSet*) w->data(TString::Format("Hig_%d_cat%d",1,c));
    higToFitvbf[c] = (RooDataSet*) w->data(TString::Format("Hig_%d_cat%d",2,c));
    higToFitvh[c] = (RooDataSet*) w->data(TString::Format("Hig_%d_cat%d",3,c));
  }
  //RooRealVar* lumi = w->var("lumi");
  cout << "======== Expected Events Number =====================" << endl;
  cout << ".........Measured Data for L = " << "19785" << " pb-1 ............................" << endl;
  cout << "#Events data: " << w->data("Data")->sumEntries() << endl;
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events data cat%d: ",c) << data[c]->sumEntries() << endl;
  }
  // cout << ".........Expected Signal for L = " << lumi->getVal() << " pb-1 ............................" << endl;
  cout << ".........Expected Signal for L = " << "19785" << " pb-1 ............................" << endl;
  cout << "#Events Signal: " << w->data("Data")->sumEntries() << endl;
  Float_t siglikeErr[ncat];
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events Signal cat%d: ",c) << sigToFit[c]->sumEntries() << endl;
    siglikeErr[c]=0.6*sigToFit[c]->sumEntries();
  }
  cout << "====================================================" << endl;
  TString filename(cardDir+TString(fileBaseName)+"onecat.txt");
  ofstream outFile(filename);

  //outFile << "#CMS-HGG DataCard for Unbinned Limit Setting, " << lumi->getVal() << " pb-1 " << endl;
  outFile << "#Run with: combine -d hgg.mH260.0.shapes-Unbinned.txt -U -m 130 -H ProfileLikelihood -M MarkovChainMC --rMin=0 --rMax=20.0 -b 3000 -i 50000 --optimizeSim=1 --tries 30" << endl;
  // outFile << "# Lumi = " << lumi->getVal() << " pb-1" << endl;
  outFile << "# Lumi = " << "19785" << " pb-1" << endl;
  outFile << "imax 1" << endl;
  outFile << "jmax 5" << endl; // number of BKG
  outFile << "kmax *" << endl;
  outFile << "---------------" << endl;
  outFile << "shapes data_obs cat0 " << TString(fileBkgName)+".root" << " w_all:data_obs_cat0" << endl;
  outFile << "shapes data_obs cat1 "<< TString(fileBkgName)+".root" << " w_all:data_obs_cat1" << endl;
  outFile << "############## shape with reparametrization" << endl;
  outFile << "shapes mggBkg cat0 " << TString(fileBkgName)+".root" << " w_all:CMS_hgg_bkg_8TeV_cat0" << endl;
  outFile << "shapes mggBkg cat1 "<< TString(fileBkgName)+".root" << " w_all:CMS_hgg_bkg_8TeV_cat1" << endl;
  outFile << "# signal" << endl;
  outFile << "shapes mggSig cat0 " << TString(fileBaseName)+".inputsig.root" << " w_all:CMS_hgg_sig_cat0" << endl;
  outFile << "shapes mggSig cat1 " << TString(fileBaseName)+".inputsig.root" << " w_all:CMS_hgg_sig_cat1" << endl;
  outFile << "# ggh" << endl;
  outFile << "shapes mggHigggh cat0 " << TString(fileHiggsNameggh)+".inputsig.root" << " w_all:CMS_hgg_hig_0_cat0" << endl;
  outFile << "# tth" << endl;
  outFile << "shapes mggHigtth cat0 " << TString(fileHiggsNametth)+".inputsig.root" << " w_all:CMS_hgg_hig_1_cat0" << endl;
  outFile << "# vbf" << endl;
  outFile << "shapes mggHigvbf cat0 " << TString(fileHiggsNamevbf)+".inputsig.root" << " w_all:CMS_hgg_hig_2_cat0" << endl;
  outFile << "# vh" << endl;
  outFile << "shapes mggHigvh cat0 " << TString(fileHiggsNamevh)+".inputsig.root" << " w_all:CMS_hgg_hig_3_cat0" << endl;
/*
  outFile << "# ggh" << endl;
  outFile << "shapes mggHigggh cat0 " << TString(fileHiggsNameggh)+".inputsig.root" << " w_all:CMS_hgg_sig_ggh_cat0" << endl;
  outFile << "shapes mggHigggh cat1 " << TString(fileHiggsNameggh)+".inputsig.root" << " w_all:CMS_hgg_sig_ggh_cat1" << endl;
  outFile << "# tth" << endl;
  outFile << "shapes mggHigtth cat0 " << TString(fileHiggsNametth)+".inputsig.root" << " w_all:CMS_hgg_sig_tth_cat0" << endl;
  outFile << "shapes mggHigtth cat1 " << TString(fileHiggsNametth)+".inputsig.root" << " w_all:CMS_hgg_sig_tth_cat1" << endl;
  outFile << "# vbf" << endl;
  outFile << "shapes mggHigvbf cat0 " << TString(fileHiggsNamevbf)+".inputsig.root" << " w_all:CMS_hgg_sig_vbf_cat0" << endl;
  outFile << "shapes mggHigvbf cat1 " << TString(fileHiggsNamevbf)+".inputsig.root" << " w_all:CMS_hgg_sig_vbf_cat1" << endl;
  outFile << "# vh" << endl;
  outFile << "shapes mggHigvh cat0 " << TString(fileHiggsNamevh)+".inputsig.root" << " w_all:CMS_hgg_sig_vh_cat0" << endl;
  outFile << "shapes mggHigvh cat1 " << TString(fileHiggsNamevh)+".inputsig.root" << " w_all:CMS_hgg_sig_vh_cat1" << endl;
*/
  outFile << "---------------" << endl; //CMS_hgg_hig_0_
  /////////////////////////////////////
  if(addHiggs) { //
  outFile << "bin cat0 " << endl;
  cout<<"here"<<endl;
  outFile << "observation "<< data[0]->sumEntries() <<" "<< endl;
  outFile << "------------------------------" << endl;
  outFile << "bin cat0 cat0 cat0 cat0 cat0 cat0" << endl;
  outFile << "process mggSig mggBkg mggHigggh mggHigtth mggHigvbf mggHigvh"<< endl;
  outFile << "process 0 1 2 3 4 5 " << endl;
  outFile << "rate "
<<" "<<sigToFit[0]->sumEntries()<<" "<<1<<" "<<higToFitggh[0]->sumEntries()<<" "<<higToFittth[0]->sumEntries()<<" "<<higToFitvbf[0]->sumEntries()<<" "<<higToFitvh[0]->sumEntries()<<" "<<endl;
  outFile << "--------------------------------" << endl;
  outFile << "############## Total normalisation" << endl;
  outFile << "lumi_8TeV lnN "
<< "1.026 - 1.026 1.026 1.026 1.026 "<< endl;
  outFile << "     " << endl;
  outFile << "############## Photon selection normalisation uncertainties " << endl;
  outFile << "DiphoTrigger lnN "
<< "1.01 - 1.010 1.010 1.010 1.010 "
<< "# Trigger efficiency" << endl;
  outFile << "CMS_hgg_eff_g lnN "
   << "1.010 - 1.010 1.010 1.010 1.010 "
   << "# photon selection accep." << endl;
  outFile << "     " << endl;
  outFile << "############## Jet selection and phase space cuts normalisation uncertainties " << endl;
  outFile << "Mjj_PTj_cut_acceptance lnN "
<< "1.015 - 1.015 1.015 1.015 1.015 "
<<"# JER and JES " << endl;
  outFile << "btag_eff lnN "
<< "1.046 - 1.046 1.046 1.046 1.046 "
<<"# b tag efficiency uncertainty" << endl;
  outFile << "maajj_cut_acceptance lnN "
   << "1.02 - 1.02 1.02 1.02 1.02 "<< endl;
  outFile << "     " << endl;
  outFile << "############## Theory uncertainties on SM Higgs production " << endl;
  outFile << "PDF lnN "
   << " - - 0.931/1.075 0.919/1.081 0.972/1.026 0.976/1.024 " << endl;
  outFile << "QCD_scale lnN "
   << " - - 0.922/1.072 0.907/1.038 0.998/1.002 0.980/1.020 "<< endl;
  outFile << "gg_migration lnN "
   << " - - 1.25 1.25 1.08 1.08 # UEPS" << endl;
  outFile << "gluonSplitting lnN "
   << " - - 1.40 1.40 1.40 1.40 "<< endl;
  outFile << "     " << endl;
  outFile << "############## Signal parametric shape uncertainties " << endl;
  outFile << "CMS_hgg_sig_m0_absShift param 1 0.0057 # displacement of the dipho mean error = sqrt(0.45^ 2 + 0.35^ 2) " << endl;
  outFile << "CMS_hgg_sig_sigmaScale param 1 0.22 # optimistic estimative of resolution uncertainty " << endl;  
  outFile << "############## for mgg fit - slopes" << endl;
  outFile << "CMS_hgg_bkg_8TeV_cat0_norm flatParam # Normalization uncertainty on background slope" << endl;

  outFile << "CMS_hgg_bkg_8TeV_slope1_cat0 flatParam # Mean and absolute uncertainty on background slope" << endl;
  } // if ncat ==2
  /////////////////////////////////////

  outFile.close();
  cout << "Write data card in: " << filename << " file" << endl;
  return;
} // close write datacard one cat
////////////////////////////////////////////////////////////////////////////////
void MakeDataCardonecatnohiggs(RooWorkspace* w, const char* fileBaseName, const char* fileBkgName , const char* fileHiggsName, const char* fileHiggsNametth, const char* fileHiggsNamevbf, const char* fileHiggsNamevh) {
  TString cardDir = "datacards/";
  const Int_t ncat = NCAT;
  RooDataSet* data[ncat];
  RooDataSet* sigToFit[ncat];
  RooDataSet* higToFit[ncat];
  for (int c = 0; c < ncat; ++c) {
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    sigToFit[c] = (RooDataSet*) w->data(TString::Format("Sig_cat%d",c));
    higToFit[c] = (RooDataSet*) w->data(TString::Format("Hig_cat%d",c));
  }
  //RooRealVar* lumi = w->var("lumi");
  cout << "======== Expected Events Number =====================" << endl;
  cout << ".........Measured Data for L = " << "19785" << " pb-1 ............................" << endl;
  cout << "#Events data: " << w->data("Data")->sumEntries() << endl;
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events data cat%d: ",c) << data[c]->sumEntries() << endl;
  }
// cout << ".........Expected Signal for L = " << lumi->getVal() << " pb-1 ............................" << endl;
  cout << ".........Expected Signal for L = " << "19785" << " pb-1 ............................" << endl;
  cout << "#Events Signal: " << w->data("Data")->sumEntries() << endl;
  Float_t siglikeErr[ncat];
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events Signal cat%d: ",c) << sigToFit[c]->sumEntries() << endl;
    siglikeErr[c]=0.6*sigToFit[c]->sumEntries();
  }
  cout << "====================================================" << endl;
  TString filename(cardDir+TString(fileBaseName)+"onecatnohiggs.txt");
  ofstream outFile(filename);

  //outFile << "#CMS-HGG DataCard for Unbinned Limit Setting, " << lumi->getVal() << " pb-1 " << endl;
  outFile << "#Run with: combine -d hgg.mH260.0.shapes-Unbinned.txt -U -m 130 -H ProfileLikelihood -M MarkovChainMC --rMin=0 --rMax=20.0 -b 3000 -i 50000 --optimizeSim=1 --tries 30" << endl;
// outFile << "# Lumi = " << lumi->getVal() << " pb-1" << endl;
  outFile << "# Lumi = " << "19785" << " pb-1" << endl;
  outFile << "imax 1" << endl;
  outFile << "jmax 1" << endl; // number of BKG
  outFile << "kmax *" << endl;
  outFile << "---------------" << endl;

cout<<"here"<<endl;
  outFile << "shapes data_obs cat0 " << TString(fileBkgName)+".root" << " w_all:data_obs_cat0" << endl;
  outFile << "############## shape with reparametrization" << endl;
  outFile << "shapes mggBkg cat0 " << TString(fileBkgName)+".root" << " w_all:CMS_hgg_bkg_8TeV_cat0" << endl;
  outFile << "# signal" << endl;
  outFile << "shapes mggSig cat0 " << TString(fileBaseName)+".inputsig.root" << " w_all:CMS_hgg_sig_cat0" << endl;
  outFile << "# signal" << endl;
  outFile << "shapes mggHig cat0 " << TString(fileHiggsName)+".inputsig.root" << " w_all:CMS_hgg_sig_cat0" << endl;

  outFile << "---------------" << endl;
  /////////////////////////////////////
  if(addHiggs) { //
  outFile << "bin cat0 " << endl;
  outFile << "observation "
<< data[0]->sumEntries() << " " << endl;
  outFile << "------------------------------" << endl;
  outFile << "bin cat0 cat0 " << endl;
  outFile << "process mggSig mggBkg " << endl;
  outFile << "process 0 1 " << endl;
  outFile << "rate "
<< " " << sigToFit[0]->sumEntries() << " " << 1
<< " " << endl;
  outFile << "--------------------------------" << endl;
  outFile << "lumi_8TeV lnN "
<< "1.022 - " << endl;
  outFile << "############## jet" << endl;
  outFile << "Mjj_acceptance lnN "
<< "1.015 - "
<<"# JER and JES " << endl;
  outFile << "btag_eff lnN "
<< "1.046 - "
<<"# b tag efficiency uncertainty" << endl;
  outFile << "############## photon " << endl;
  outFile << "CMS_hgg_eff_g lnN "
   << "1.010 - "
   << "# photon selection accep." << endl;
  outFile << "DiphoTrigger lnN "
<< "1.01 - "
<< "# Trigger efficiency" << endl;
  outFile << "############## for mtot fit" << endl;
  outFile << "maa_acceptance lnN "
   << "1.10 - "
   << "# photon energy resolution" << endl;
  outFile << "# Parametric shape uncertainties, entered by hand. they act on signal " << endl;
  outFile << "CMS_hgg_sig_m0_absShift param 1 0.0045 # displacement of the dipho mean" << endl;
  outFile << "CMS_hgg_sig_sigmaScale param 1 0.22 # optimistic estimative of resolution uncertainty " << endl;
  outFile << "# Parametric shape uncertainties, entered by hand. they act on higgs " << endl;
  outFile << "CMS_hgg_hig_m0_absShift param 1 0.045 # displacement of the dipho mean" << endl;
  outFile << "CMS_hgg_hig_sigmaScale param 1 0.22 # optimistic estimative of resolution uncertainty " << endl;
  outFile << "############## for mgg fit - slopes" << endl;
  outFile << "CMS_hgg_bkg_8TeV_cat0_norm flatParam # Normalization uncertainty on background slope" << endl;
  outFile << "CMS_hgg_bkg_8TeV_cat1_norm flatParam # Normalization uncertainty on background slope" << endl;

  outFile << "CMS_hgg_bkg_8TeV_slope1_cat0 flatParam # Mean and absolute uncertainty on background slope" << endl;


  } // if ncat ==2
  /////////////////////////////////////

  outFile.close();
  cout << "Write data card in: " << filename << " file" << endl;
  return;
} // close write datacard one cat



void style(){
  TStyle *defaultStyle = new TStyle("defaultStyle","Default Style");
//  defaultStyle->SetOptStat(0000);
//  defaultStyle->SetOptFit(000); 
//  defaultStyle->SetPalette(1);
  /////// pad ////////////
  defaultStyle->SetPadBorderMode(1);
  defaultStyle->SetPadBorderSize(1);
  defaultStyle->SetPadColor(0);
  defaultStyle->SetPadTopMargin(0.08);
  defaultStyle->SetPadBottomMargin(0.15);
  defaultStyle->SetPadLeftMargin(0.20);
  defaultStyle->SetPadRightMargin(0.04);
  /////// canvas /////////
  defaultStyle->SetCanvasBorderMode(0);
  defaultStyle->SetCanvasColor(0);
//  defaultStyle->SetCanvasDefH(600);
//  defaultStyle->SetCanvasDefW(600);
  /////// frame //////////
  defaultStyle->SetFrameBorderMode(0);
  defaultStyle->SetFrameBorderSize(1);
  defaultStyle->SetFrameFillColor(0); 
  defaultStyle->SetFrameLineColor(1);
  /////// label //////////
//  defaultStyle->SetLabelOffset(0.005,"XY");
  defaultStyle->SetLabelSize(0.05,"XY");
  defaultStyle->SetLabelFont(42,"XY");
  /////// title //////////
//  defaultStyle->SetTitleOffset(1.1,"X");
//  defaultStyle->SetTitleSize(0.01,"X");
  defaultStyle->SetTitleOffset(1.25,"Y");
//  defaultStyle->SetTitleSize(0.05,"Y");
  defaultStyle->SetTitleFont(42, "XYZ");
  /////// various ////////
  defaultStyle->SetNdivisions(303,"Y");
  //defaultStyle->SetTitleFillStyle(10, "Z");

//  defaultStyle->SetLegendBorderSize(0);  // For the axis titles:

//    defaultStyle->SetTitleColor(1, "XYZ");
//    defaultStyle->SetTitleFont(42, "XYZ");
    defaultStyle->SetTitleSize(0.06, "XYZ");
 
    // defaultStyle->SetTitleYSize(Float_t size = 0.02);
    //defaultStyle->SetTitleXOffset(0.9);
    //defaultStyle->SetTitleYOffset(1.05);
    // defaultStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

    // For the axis labels:
    defaultStyle->SetLabelColor(1, "XYZ");
    defaultStyle->SetLabelFont(42, "XYZ");
   // defaultStyle->SetLabelOffset(0.007, "XYZ");
    defaultStyle->SetLabelSize(0.045, "XYZ");

    // For the axis:
//    defaultStyle->SetAxisColor(1, "XYZ");
    defaultStyle->SetStripDecimals(kTRUE);
    defaultStyle->SetTickLength(0.03, "XYZ");
    defaultStyle->SetNdivisions(510, "XYZ");
    defaultStyle->SetPadTickX(1);   // To get tick marks on the opposite side of the frame
    defaultStyle->SetPadTickY(1);
    defaultStyle->cd();
  return;
}

void SetParamNames(RooWorkspace* w) { // not used it if Workspaces are created => float fit
  const Int_t ncat = NCAT;
  //****************************//
  // mgg signal all categories
  //****************************//
  RooRealVar* mgg_sig_m0 = w->var("mgg_sig_m0");
  RooRealVar* mgg_sig_sigma = w->var("mgg_sig_sigma");
  RooRealVar* mgg_sig_alpha = w->var("mgg_sig_alpha");
  RooRealVar* mgg_sig_n = w->var("mgg_sig_n");
  RooRealVar* mgg_sig_gsigma = w->var("mgg_sig_gsigma");
  RooRealVar* mgg_sig_frac = w->var("mgg_sig_frac");
  mgg_sig_m0 ->SetName("m_{0}");
  mgg_sig_sigma ->SetName("#sigma_{CB}");
  mgg_sig_alpha ->SetName("#alpha");
  mgg_sig_n ->SetName("n");
  mgg_sig_gsigma->SetName("#sigma_G");
  mgg_sig_frac ->SetName("f_G");
  mgg_sig_m0 ->setUnit("GeV");
  mgg_sig_sigma ->setUnit("GeV");
  mgg_sig_gsigma->setUnit("GeV");
  //****************************//
  // mgg background
  //****************************//
  RooRealVar* mgg_bkg_8TeV_slope1 = w->var("mgg_bkg_8TeV_slope1");
  mgg_bkg_8TeV_slope1 ->SetName("a_{B}");
  mgg_bkg_8TeV_slope1 ->setUnit("1/GeV");
  RooRealVar* mgg_bkg_8TeV_slope2 = w->var("mgg_bkg_8TeV_slope2");
  mgg_bkg_8TeV_slope2 ->SetName("a_{B}");
  mgg_bkg_8TeV_slope2 ->setUnit("1/GeV");
  //****************************//
  // mgg per category
  //****************************//
  for (int c = 0; c < ncat; ++c) {
    mgg_sig_m0 = (RooRealVar*) w->var(TString::Format("mgg_sig_m0_cat%d",c));
    mgg_sig_sigma = (RooRealVar*) w->var(TString::Format("mgg_sig_sigma_cat%d",c));
    mgg_sig_alpha = (RooRealVar*) w->var(TString::Format("mgg_sig_alpha_cat%d",c));
    mgg_sig_n = (RooRealVar*) w->var(TString::Format("mgg_sig_n_cat%d",c));
    mgg_sig_gsigma = (RooRealVar*) w->var(TString::Format("mgg_sig_gsigma_cat%d",c));
    mgg_sig_frac = (RooRealVar*) w->var(TString::Format("mgg_sig_frac_cat%d",c));
    mgg_sig_m0 ->SetName("m_{0}");
    mgg_sig_sigma ->SetName("#sigma_{CB}");
    mgg_sig_alpha ->SetName("#alpha");
    mgg_sig_n ->SetName("n");
    mgg_sig_gsigma ->SetName("#sigma_{G}");
    mgg_sig_frac ->SetName("f_{G}");
    mgg_sig_m0 ->setUnit("GeV");
    mgg_sig_sigma ->setUnit("GeV");
    mgg_sig_gsigma ->setUnit("GeV");
    mgg_bkg_8TeV_slope1 = w->var(TString::Format("mgg_bkg_8TeV_slope1_cat%d",c));
    mgg_bkg_8TeV_slope1 ->SetName("p_{B}^{1}");
    mgg_bkg_8TeV_slope1 ->setUnit("1/GeV");
    mgg_bkg_8TeV_slope2 = w->var(TString::Format("mgg_bkg_8TeV_slope2_cat%d",c));
    mgg_bkg_8TeV_slope2 ->SetName("p_{B}^{2}");
    mgg_bkg_8TeV_slope2 ->setUnit("1/GeV^{2}");
// RooRealVar* mgg_bkg_8TeV_frac = w->var(TString::Format("mgg_bkg_8TeV_frac_cat%d",c));
// mgg_bkg_8TeV_frac ->SetName("f");
  }
} // close setparameters

