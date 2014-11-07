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
bool addHiggs=true;
void AddSigData(RooWorkspace*, Float_t);
void AddHigData(RooWorkspace*, Float_t,int);
void AddBkgData(RooWorkspace*);
void SigModelFit(RooWorkspace*, Float_t);
void HigModelFit(RooWorkspace*, Float_t, int);
void MakePlots(RooWorkspace*, Float_t);
void MakePlotsHiggs(RooWorkspace* w, Float_t Mass);
void MakeSigWS(RooWorkspace* w, const char* filename);
void MakeHigWS(RooWorkspace* w, const char* filename,int);
void MakeBkgWS(RooWorkspace* w, const char* filename);//,
// const char* filenameh0, const char* filenameh1, const char* filenameh2, const char* filenameh4);
void MakeDataCard(RooWorkspace* w, const char* filename, const char* filename1,
                  const char*, const char*, const char*, const char*, const char*);
void MakeDataCardonecatnohiggs(RooWorkspace* w, const char* filename, const char* filename1, const char* filename2);
void SetParamNames(RooWorkspace*);
void SetConstantParams(const RooArgSet* params);

RooFitResult* fitresult[NCAT]; // container for the fit results
RooFitResult* BkgModelFitBernstein(RooWorkspace*, Bool_t);

Bool_t doblinding = false; //True if you want to blind

RooArgSet* defineVariables()
{
  RooRealVar* mgg = new RooRealVar("mgg","M(#gamma#gamma)",100,180,"GeV");
  //RooRealVar* mtot = new RooRealVar("mtot","M(#gamma#gammajj)",200,1600,"GeV");
  //RooRealVar* mjj = new RooRealVar("mjj","M(jj)",100,1600,"GeV");
  RooRealVar* evWeight = new RooRealVar("evWeight","HqT x PUwei",0,100,"");
  RooCategory* cut_based_ct = new RooCategory("cut_based_ct","event category 4") ;
  //
  cut_based_ct->defineType("cat4_0",0);
  cut_based_ct->defineType("cat4_1",1);
  //
  RooArgSet* ntplVars = new RooArgSet(*mgg, * cut_based_ct, *evWeight);
  ntplVars->add(*mgg);
  //ntplVars->add(*mtot);
  //ntplVars->add(*mjj);
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
  TString fileHiggsNamebbh(TString::Format("hgg.hig.mH%.1f_8TeV.bbh", mass));
  TString fileBkgName(TString::Format("hgg.inputbkg_8TeV", mass));
  TString card_name("models_test.rs"); // put the model parameters here!
  HLFactory hlf("HLFactory", card_name, false);
  RooWorkspace* w = hlf.GetWs();
  RooFitResult* fitresults;
  bool cutbased=true;
  // the minitree to be addeed
  //
  TString hhiggsggh = "/afs/cern.ch/user/h/hebda/public/forRadion/v35_fitToMgg_resSearch_withKinFit/ggh_m125_powheg_8TeV_m350.root";
  TString hhiggstth = "/afs/cern.ch/user/h/hebda/public/forRadion/v35_fitToMgg_resSearch_withKinFit/tth_m125_8TeV_m350.root";
  TString hhiggsvbf = "/afs/cern.ch/user/h/hebda/public/forRadion/v35_fitToMgg_resSearch_withKinFit/vbf_m125_8TeV_m350.root";
  TString hhiggsvh = "/afs/cern.ch/user/h/hebda/public/forRadion/v35_fitToMgg_resSearch_withKinFit/wzh_m125_8TeV_zh_m350.root";
  TString hhiggsbbh = "/afs/cern.ch/user/h/hebda/public/forRadion/v35_fitToMgg_resSearch_withKinFit/bbh_m125_8TeV_m350.root";
  //
  TString ssignal = "/afs/cern.ch/user/h/hebda/public/forRadion/v35_fitToMgg_resSearch_withKinFit/Radion_m350_8TeV_m350.root ";
  TString ddata = "/afs/cern.ch/user/h/hebda/public/forRadion/v35_fitToMgg_resSearch_withKinFit/Data_m350.root";
  //
  // TString hhiggs = "MiniTrees/OlivierOc13/v16_base_mgg_0_massCutVersion0/02013-11-05-Radion_m350_8TeV_nm_m350.root";
  // TString ssignal = "MiniTrees/OlivierOc13/v16_base_mgg_0_massCutVersion0/02013-11-05-Radion_m350_8TeV_nm_m350.root";
  // TString ddata = "MiniTrees/OlivierOc13/v16_base_mgg_0_massCutVersion0/02013-11-05-Data_m350.root";
  //
  // TString hhiggs = "MiniTrees/OlivierOc13/v15_regkin_mgg_0_massCutVersion0/02013-10-30-Radion_m350_8TeV_nm_m350.root";
  // TString ssignal = "MiniTrees/OlivierOc13/v15_regkin_mgg_0_massCutVersion0/02013-10-30-Radion_m350_8TeV_nm_m350.root";
  // TString ddata = "MiniTrees/OlivierOc13/v15_regkin_mgg_0_massCutVersion0/02013-10-30-Data_m350.root";
  //
  cout<<"Signal: "<<ssignal<<endl;
  cout<<"Data: "<<ddata<<endl;
  //

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
  cout<<"=============================== Higgs: ============================= "<<hhiggsvh<<endl;
  AddHigData(w, mass,hhiggsvh,3);
  HigModelFit(w, mass,3); // constructing higgs pdf
  MakeHigWS(w, fileHiggsNamevh,3);
  cout<<"HIGGS ADDED"<<endl;
  //
  cout<<"=============================== Higgs: ============================= "<<hhiggsbbh<<endl;
  AddHigData(w, mass,hhiggsbbh,4);
  HigModelFit(w, mass,4); // constructing higgs pdf
  MakeHigWS(w, fileHiggsNamebbh,4);
  cout<<"HIGGS ADDED"<<endl;
  MakePlotsHiggs(w, mass);
  //


  AddBkgData(w,ddata);
  w->Print("v");
  cout<<"BKG ADDED"<<endl;
  bool dobands=true;
  fitresults = BkgModelFitBernstein(w, dobands); // this is berestein 3
  MakeBkgWS(w, fileBkgName);
  // construct the models to fit
  //
  // MakeDataCardonecat(w, fileBaseName, fileBkgName, fileHiggsNameggh, fileHiggsNametth, fileHiggsNamevbf, fileHiggsNamevh);
  MakeDataCardonecatnohiggs(w, fileBaseName, fileBkgName);
  MakeDataCard(w, fileBaseName, fileBkgName, fileHiggsNameggh, fileHiggsNametth, fileHiggsNamevbf, fileHiggsNamevh, fileHiggsNamebbh);
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
  Float_t Lum = 19712.0; // pb-1
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
					       *w->var("mgg"),
					       mainCut+TString::Format(" && cut_based_ct==%d ",0)+cut0+cutj0);
  w->import(*sigToFit[0],Rename(TString::Format("Sig_cat%d",0)));
  //
  sigToFit[1] = (RooDataSet*) sigScaled.reduce(
					       *w->var("mgg"),
					       mainCut+TString::Format(" && cut_based_ct==%d ",1)+cut1+cutj1);
  w->import(*sigToFit[1],Rename(TString::Format("Sig_cat%d",1)));
  // Create full signal data set without categorization
  RooDataSet* sigToFitAll = (RooDataSet*) sigScaled.reduce(*w->var("mgg"),mainCut);
  cout << "======================================================================" << endl;
  w->import(*sigToFitAll,Rename("Sig"));
  // here we print the number of entries on the different categories
  cout << "========= the number of entries on the different categories ==========" << endl;
  cout << "---- one channel: " << sigScaled.sumEntries() << endl;
  for (int c = 0; c < ncat; ++c) {
    Float_t nExpEvt = sigToFitAll[c]->sumEntries();
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
  RooDataSet Data("Data","dataset",dataTree,*ntplVars,"","weightVar");
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
					   *w->var("mgg"),
					   mainCut+TString::Format(" && cut_based_ct==%d",0)+cut0+cutj0);
  if(doblinding){ dataToPlot[0] = (RooDataSet*) Data.reduce(
					    *w->var("mgg"),
					    mainCut+TString::Format(" && cut_based_ct==%d",0)
					    +TString::Format(" && (mgg > 130 || mgg < 120)")// blind
					    +cut0+cutj0);
  }else{

                  dataToPlot[0] = (RooDataSet*) Data.reduce(
					    *w->var("mgg"),
					    mainCut+TString::Format(" && cut_based_ct==%d",0)
					    +cut0+cutj0);

  }
   
  dataToFit[1] = (RooDataSet*) Data.reduce(
					   *w->var("mgg"),
					   mainCut+TString::Format(" && cut_based_ct==%d",1)+cut1);
  if(doblinding){ dataToPlot[1] = (RooDataSet*) Data.reduce(
					    *w->var("mgg"),
					    mainCut+TString::Format(" && cut_based_ct==%d",1)
					    +TString::Format(" && (mgg > 130 || mgg < 120)") // blind
					    +cut1);
  }else{
                  dataToPlot[1] = (RooDataSet*) Data.reduce(
					    *w->var("mgg"),
					    mainCut+TString::Format(" && cut_based_ct==%d",1)
					    +cut1);  
  }

  for (int c = 0; c < ncat; ++c) {
    w->import(*dataToFit[c],Rename(TString::Format("Data_cat%d",c)));
    w->import(*dataToPlot[c],Rename(TString::Format("Dataplot_cat%d",c)));
  }
  // Create full data set without categorization
  RooDataSet* data = (RooDataSet*) Data.reduce(*w->var("mgg"),mainCut);
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
  // fit range
  /*
    const int minsigfit =mass - 120, maxsigfit=mass +120;
    RooDataSet* sigToFit[ncat];
    RooAbsPdf* mtotSig[ncat];
    // fit range
    Float_t minMassFit2(minfit2),minMassFit1(minsigfit),maxMassFit(maxsigfit);
  */
  Float_t minSigFit(115),maxSigFit(135);
  for (int c = 0; c < ncat; ++c) {
    // import sig and data from workspace
    sigToFit[c] = (RooDataSet*) w->data(TString::Format("Sig_cat%d",c));
    mggSig[c] = (RooAbsPdf*) w->pdf(TString::Format("mggSig_cat%d",c));
    ((RooRealVar*) w->var(TString::Format("mgg_sig_m0_cat%d",c)))->setVal(MASS);
    //RooRealVar* peak = w->var(TString::Format("mgg_sig_m0_cat%d",c));
    //peak->setVal(MASS);
    cout << "OK up to now..." <<MASS<< endl;
    // Fit model as M(x|y) to D(x,y)

    mggSig[c]->fitTo(*sigToFit[c],Range(minSigFit,maxSigFit),SumW2Error(kTRUE));
    cout << "old = " << ((RooRealVar*) w->var(TString::Format("mgg_sig_m0_cat%d",c)))->getVal() << endl;

    double mPeak = ((RooRealVar*) w->var(TString::Format("mgg_sig_m0_cat%d",c)))->getVal()+0.6; // shift the peak
    ((RooRealVar*) w->var(TString::Format("mgg_sig_m0_cat%d",c)))->setVal(mPeak); // shift the peak

    cout << "mPeak = " << mPeak << endl;
    cout << "new mPeak position = " << ((RooRealVar*) w->var(TString::Format("mgg_sig_m0_cat%d",c)))->getVal() << endl;

    // IMPORTANT: fix all pdf parameters to constant, why?
    w->defineSet(TString::Format("SigPdfParam_cat%d",c),
		 RooArgSet(
			   *w->var(TString::Format("mgg_sig_m0_cat%d",c)),
			   *w->var(TString::Format("mgg_sig_sigma_cat%d",c)),
			   *w->var(TString::Format("mgg_sig_alpha_cat%d",c)),
			   *w->var(TString::Format("mgg_sig_n_cat%d",c)),
			   *w->var(TString::Format("mgg_sig_gsigma_cat%d",c)),
			   *w->var(TString::Format("mgg_sig_frac_cat%d",c))) );
    SetConstantParams(w->set(TString::Format("SigPdfParam_cat%d",c)));
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
  // fit range
  Float_t minSigFit(115),maxSigFit(135);
  for (int c = 0; c < ncat; ++c) {
    // import sig and data from workspace
    higToFit[c] = (RooDataSet*) w->data(TString::Format("Hig_%d_cat%d",higgschannel,c));
    mggHig[c] = (RooAbsPdf*) w->pdf(TString::Format("mggHig_%d_cat%d",higgschannel,c));
    //RooRealVar* peak = w->var(TString::Format("mgg_hig_m0_cat%d",c));
    //peak->setVal(MASS);
    ((RooRealVar*) w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)))->setVal(MASS);
    cout << "OK up to now..." <<MASS<< endl;
    cout << "old = " << ((RooRealVar*) w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)))->getVal() << endl;

    double mPeak = ((RooRealVar*) w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)))->getVal()+0.6; // shift the peak
    ((RooRealVar*) w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)))->setVal(mPeak); // shift the peak

    cout << "mPeak = " << mPeak << endl;
    cout << "new mPeak position = " << ((RooRealVar*) w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)))->getVal() << endl;

    // Fit model as M(x|y) to D(x,y)
    mggHig[c]->fitTo(*higToFit[c],Range(minSigFit,maxSigFit),SumW2Error(kTRUE));
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
} // close higgs model fit
////////////////////////////////////////////////////////////
// BKG model berestein 3
RooFitResult* BkgModelFitBernstein(RooWorkspace* w, Bool_t dobands) {
  const Int_t ncat = NCAT;
  std::vector<TString> catdesc;
  catdesc.push_back("2 btag");
  catdesc.push_back("1 btag");
  catdesc.push_back("cat 2");
  catdesc.push_back("cat 3");
  //******************************************//
  // Fit background with model pdfs
  //******************************************//
  // retrieve pdfs and datasets from workspace to fit with pdf models
  RooDataSet* data[ncat];
  RooDataSet* dataplot[ncat]; // the data
  RooBernstein* mggBkg[ncat];// the polinomial of 4* order
  RooGaussian* Higgs[ncat]; // the higgs to sum
  RooPlot* plotmggBkg[ncat];
  RooDataSet* sigToFit0[ncat];
  RooDataSet* sigToFit1[ncat];
  RooDataSet* sigToFit2[ncat];
  RooDataSet* sigToFit3[ncat];
  RooAbsPdf* mggSig[ncat];
  Float_t minMassFit(100),maxMassFit(180);
  // Fit data with background pdf for data limit
  RooRealVar* mgg = w->var("mgg");
  mgg->setUnit("GeV");
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
    RooFormulaVar *p1mod = new RooFormulaVar(
					     TString::Format("p1mod_cat%d",c),
					     "","@0*@0",
					     *w->var(TString::Format("mgg_bkg_8TeV_slope1_cat%d",c)));
    // if(c==1){RooFormulaVar *p2mod = new RooFormulaVar(
    //TString::Format("p2mod_cat%d",c)
    //,"","@0*@0",
    //*w->var(TString::Format("mgg_bkg_8TeV_slope2_cat%d",c)));
    // RooFormulaVar *p3mod = new RooFormulaVar(
    //TString::Format("p3mod_cat%d",c)
    //,"","@0*@0",
    //*w->var(TString::Format("mgg_bkg_8TeV_slope3_cat%d",c)));
    /*
      RooFormulaVar *p4mod = new RooFormulaVar(
      TString::Format("p4mod_cat%d",c)
      ,"","@0*@0",
      *w->var(TString::Format("mgg_bkg_8TeV_slope4_cat%d",c)));
      RooFormulaVar *p5mod = new RooFormulaVar(
      TString::Format("p5mod_cat%d",c)
      ,"","@0*@0",
      *w->var(TString::Format("mgg_bkg_8TeV_slope5_cat%d",c)));
      */
    // }
    ////////////////////////////////////////////////////////////////////
    /* RooAbsPdf* mggBkgTmp0 = 0; // declare a empty pdf
    // adding pdf's, using the variables
    if(c==0){
    mggBkgTmp0 = new RooBernstein( // fill the pdf with the floating parameters
    TString::Format("mggBkgTmp0_cat%d",c),
    "", *mgg,
    RooArgList(RooConst(1.0),*p1mod));
    }
    if(c==1){
    mggBkgTmp0 = new RooBernstein( // fill the pdf with the floating parameters
    TString::Format("mggBkgTmp0_cat%d",c),
    "", *mgg,
    RooArgList(RooConst(1.0),*p1mod, *p2mod, *p3mod));//, *p4mod, *p5mod));
    }
    */
    RooAbsPdf* mggBkgTmp0 = new RooGenericPdf( // if exp function
					      TString::Format("DijetBackground_%d",c),
					      "1./pow(@0,@1)",
					      RooArgList(*mgg, *p1mod));
    // we first wrap the normalization of mggBkgTmp0
    w->factory(TString::Format("mgg_bkg_8TeV_norm_cat%d[1.0,0.0,100000]",c));
    RooExtendPdf mggBkgTmp( // we copy the pdf? normalized
			   TString::Format("mggBkg_cat%d",c),
			   "",*mggBkgTmp0,
			   *w->var(TString::Format("mgg_bkg_8TeV_norm_cat%d",c)) // normalization only on full bkg
			   );
    fitresult[c] = mggBkgTmp.fitTo( // fit with normalized pdf,and return values
				   *data[c], // bkg
				   Strategy(1), // MINUIT strategy
				   Minos(kFALSE), // interpretation on the errors, nonlinearities
				   Range(minMassFit,maxMassFit),
				   SumW2Error(kTRUE),
				   Save(kTRUE));
    w->import(mggBkgTmp);
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
    mggBkgTmp.plotOn(
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
    plotmggBkg[c]->SetTitle("CMS preliminary 19.7/fb");
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
      mgg->setRange("errRange",minMassFit,maxMassFit);
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
    ////////////////////////////////////////////////////////// plot higgs
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

    sigToFit4[c] = (RooDataSet*) w->data(TString::Format("Hig_4_cat%d",c));
    double norm4 = 1.0*sigToFit1[c]->sumEntries(); //
    mggSig4[c] = (RooAbsPdf*) w->pdf(TString::Format("mggHig_4_cat%d",c));
    // we are not constructing signal pdf, this is constructed on sig to fit function...
    mggSig4[c] ->plotOn(
			plotmggBkg[c],
			Normalization(norm4,RooAbsPdf::NumEvent),
			DrawOption("F"),
			LineColor(kBlue),FillStyle(1001),FillColor(19));
    mggSig4[c]->plotOn(
		       plotmggBkg[c],
		       Normalization(norm4,RooAbsPdf::NumEvent),LineColor(kBlue),LineStyle(1));

    //////////////////////////////////////////////////////////
    plotmggBkg[c]->Draw("SAME");
    if(c==0)plotmggBkg[c]->SetMinimum(0.005); // no error bar in bins with zero events
    if(c==1)plotmggBkg[c]->SetMinimum(0.01); // no error bar in bins with zero events
    if(c==0)plotmggBkg[c]->SetMaximum(5.3); // no error bar in bins with zero events
    if(c==1)plotmggBkg[c]->SetMaximum(20); // no error bar in bins with zero events
    // plotmggBkg[c]->SetMinimum(0.005); // no error bar in bins with zero events
    //plotmggBkg[c]->SetLogy(0);
    cout << "!!!!!!!!!!!!!!!!!" << endl;
    TLegend *legmc = new TLegend(0.40,0.72,0.62,0.9);
    TLegend *legmcH = new TLegend(0.66,0.72,0.94,0.9);
    legmc->AddEntry(plotmggBkg[c]->getObject(2),"Data ","LPE"); // not...
    legmc->AddEntry(plotmggBkg[c]->getObject(1),"Fit","L");
    if(dobands)legmc->AddEntry(twosigma,"two sigma ","F"); // not...
    if(dobands)legmc->AddEntry(onesigma,"one sigma","F");
    legmcH->AddEntry(plotmggBkg[c]->getObject(3),"ggH ","LPE"); // not...
    legmcH->AddEntry(plotmggBkg[c]->getObject(5),"ttH ","LPE"); // not...
    legmcH->AddEntry(plotmggBkg[c]->getObject(7),"VBF ","LPE"); // not...
    legmcH->AddEntry(plotmggBkg[c]->getObject(9),"VH ","LPE"); // not...
    legmcH->AddEntry(plotmggBkg[c]->getObject(11),"bbH ","LPE"); // not...
    legmc->SetHeader(" 350 GeV");
    legmcH->SetHeader(" Higgs");
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();
    legmcH->Draw();
    TLatex *lat2 = new TLatex(minMassFit+1.5,0.75*plotmggBkg[c]->GetMaximum(),catdesc.at(c));
    lat2->Draw();
    //
    ctmp->SaveAs(TString::Format("databkgoversig_cat%d.pdf",c));
    ctmp->SaveAs(TString::Format("databkgoversig_cat%d.png",c));

    if(c==0)plotmggBkg[c]->SetMaximum(100); // no error bar in bins with zero events
    if(c==1)plotmggBkg[c]->SetMaximum(1000); // no error bar in bins with zero events
    ctmp->SetLogy(1);
    ctmp->SaveAs(TString::Format("databkgoversig_cat%d_log.pdf",c));
    ctmp->SaveAs(TString::Format("databkgoversig_cat%d_log.png",c));
    // ctmp->SaveAs(TString::Format("databkgoversig_cat%d.C",c));
  } // close to each category
  RooBernstein mggBkgAll("mggBkgAll", "", *mgg,
			 RooArgList(RooConst(1.0),
				    *w->var("mgg_bkg_8TeV_slope1"),
				    //*w->var("mgg_bkg_8TeV_slope2"),
				    *w->var("mgg_bkg_8TeV_slope3")));
  w->import(mggBkgAll);
  RooFitResult* fitresults;
  fitresults = w->pdf("mggBkgAll")->fitTo( // save results to workspace
					  *w->data("Data"),
					  Range(minMassFit,maxMassFit),
					  SumW2Error(kTRUE), Save(kTRUE));
  fitresults->Print();
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
  RooAbsPdf* mggSigPdf[ncat];
  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");
  for (int c = 0; c < ncat; ++c) {
    mggSigPdf[c] = (RooAbsPdf*) w->pdf(TString::Format("mggSig_cat%d",c));
    wAll->import(*w->pdf(TString::Format("mggSig_cat%d",c)));
  }
  // (2) Systematics on energy scale and resolution
  // 1,1,1 statistical to be treated on the datacard
  wAll->factory("CMS_hgg_sig_m0_absShift[1,1,1]");
  wAll->factory("prod::CMS_hgg_sig_m0_cat0(mgg_sig_m0_cat0, CMS_hgg_sig_m0_absShift)");
  wAll->factory("prod::CMS_hgg_sig_m0_cat1(mgg_sig_m0_cat1, CMS_hgg_sig_m0_absShift)");
  // (3) Systematics on resolution
  wAll->factory("CMS_hgg_sig_sigmaScale[1,1,1]");
  wAll->factory("prod::CMS_hgg_sig_sigma_cat0(mgg_sig_sigma_cat0, CMS_hgg_sig_sigmaScale)");
  //wAll->factory("prod::CMS_hgg_sig_sigma_cat0(mgg_sig_sigma_cat0, CMS_hgg_sig_sigmaScale)");

  wAll->factory("prod::CMS_hgg_sig_sigma_cat1(mgg_sig_sigma_cat1, CMS_hgg_sig_sigmaScale)");
  wAll->factory("prod::CMS_hgg_sig_gsigma_cat0(mgg_sig_gsigma_cat0, CMS_hgg_sig_sigmaScale)");
  wAll->factory("prod::CMS_hgg_sig_gsigma_cat1(mgg_sig_gsigma_cat1, CMS_hgg_sig_sigmaScale)");
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
					       TString::Format("EDIT::CMS_hgg_sig_cat%d(mggSig_cat%d,",c,c) +
					       TString::Format(" mgg_sig_m0_cat%d=CMS_hgg_sig_m0_cat%d, ", c,c) +
					       TString::Format(" mgg_sig_sigma_cat%d=CMS_hgg_sig_sigma_cat%d, ", c,c) +
					       // TString::Format(" mgg_sig_alpha_cat%d=CMS_hgg_sig_alpha_cat%d, ", c,c) +
					       // TString::Format(" mgg_sig_n_cat%d=CMS_hgg_sig_n_cat%d, ", c,c) +
					       // TString::Format(" mgg_sig_frac_cat%d=CMS_hgg_sig_frac_cat%d, ", c,c) +
					       TString::Format(" mgg_sig_gsigma_cat%d=CMS_hgg_sig_gsigma_cat%d)", c,c)
					       );
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
  RooAbsPdf* mggBkgPdf[ncat];
  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");
  for (int c = 0; c < ncat; ++c) {
    cout<<"here"<<endl;
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));

    //RooDataHist* dataBinned = data[c]->binnedClone(); // Uncomment if you want to use wights in the limits

    mggBkgPdf[c] = (RooAbsPdf*) w->pdf(TString::Format("mggBkg_cat%d",c));
    wAll->import(*data[c], Rename(TString::Format("data_obs_cat%d",c)));// Comment if you want to use wights in the limits

    //wAll->import(*dataBinned, Rename(TString::Format("data_obs_cat%d",c))); // Uncomment if you want to use wights in the limits

    cout<<"here"<<endl;
    wAll->import(*w->pdf(TString::Format("mggBkg_cat%d",c)));
    cout<<"here"<<endl;
    wAll->factory(
		  TString::Format("CMS_hgg_bkg_8TeV_cat%d_norm[%g,0.0,100000.0]",
				  c, w->var(TString::Format("mgg_bkg_8TeV_norm_cat%d",c))->getVal()));
    cout<<"here"<<endl;
    wAll->factory(
		  TString::Format("CMS_hgg_bkg_8TeV_slope1_cat%d[%g,-10,10]",
				  c, w->var(TString::Format("mgg_bkg_8TeV_slope1_cat%d",c))->getVal()));
    cout<<"here"<<endl;
  } // close ncat
  // (2) do reparametrization of background
  for (int c = 0; c < ncat; ++c){
    if(c==0) wAll->factory(
			   TString::Format("EDIT::CMS_hgg_bkg_8TeV_cat%d(mggBkg_cat%d,",c,c) +
			   TString::Format(" mgg_bkg_8TeV_norm_cat%d=CMS_hgg_bkg_8TeV_cat%d_norm,", c,c)+
			   TString::Format(" mgg_bkg_8TeV_slope1_cat%d=CMS_hgg_bkg_8TeV_slope1_cat%d)", c,c)
			   //TString::Format(" mgg_bkg_8TeV_slope2_cat%d=CMS_hgg_bkg_8TeV_slope2_cat%d,", c,c)+
			   //TString::Format(" mgg_bkg_8TeV_slope3_cat%d=CMS_hgg_bkg_8TeV_slope3_cat%d)", c,c)
			   );
    if(c==1) wAll->factory(
			   TString::Format("EDIT::CMS_hgg_bkg_8TeV_cat%d(mggBkg_cat%d,",c,c) +
			   TString::Format(" mgg_bkg_8TeV_norm_cat%d=CMS_hgg_bkg_8TeV_cat%d_norm,", c,c)+
			   TString::Format(" mgg_bkg_8TeV_slope1_cat%d=CMS_hgg_bkg_8TeV_slope1_cat%d)", c,c)//+
			   //TString::Format(" mgg_bkg_8TeV_slope2_cat%d=CMS_hgg_bkg_8TeV_slope2_cat%d,", c,c)+
			   //TString::Format(" mgg_bkg_8TeV_slope3_cat%d=CMS_hgg_bkg_8TeV_slope3_cat%d)", c,c)
			   // TString::Format(" mgg_bkg_8TeV_slope4_cat%d=CMS_hgg_bkg_8TeV_slope4_cat%d,", c,c)+
			   // TString::Format(" mgg_bkg_8TeV_slope5_cat%d=CMS_hgg_bkg_8TeV_slope5_cat%d)", c,c)
			   );
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
  catdesc.push_back(" 2 btag");
  catdesc.push_back(" 1 btag");
  catdesc.push_back("cat 2");
  catdesc.push_back("cat 3");
  // retrieve data sets from the workspace
  // RooDataSet* dataAll = (RooDataSet*) w->data("Data");
  RooDataSet* signalAll = (RooDataSet*) w->data("Sig");
  //RooDataSet* higgsAll = (RooDataSet*) w->data("Hig");
  // blinded dataset
  // RooDataSet* data[ncat];
  RooDataSet* sigToFit[ncat];
  RooAbsPdf* mggGaussSig[ncat];
  RooAbsPdf* mggCBSig[ncat];
  RooAbsPdf* mggSig[ncat];
  //
  RooAbsPdf* mggBkg[ncat];
  for (int c = 0; c < ncat; ++c) {
    // data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    sigToFit[c] = (RooDataSet*) w->data(TString::Format("Sig_cat%d",c));
    mggGaussSig[c] = (RooAbsPdf*) w->pdf(TString::Format("mggGaussSig_cat%d",c));
    mggCBSig[c] = (RooAbsPdf*) w->pdf(TString::Format("mggCBSig_cat%d",c));
    mggSig[c] = (RooAbsPdf*) w->pdf(TString::Format("mggSig_cat%d",c));
    mggBkg[c] = (RooAbsPdf*) w->pdf(TString::Format("mggBkg_cat%d",c));
  } // close categories
  RooRealVar* mgg = w->var("mgg");
  mgg->setUnit("GeV");
  RooAbsPdf* mggGaussSigAll = w->pdf("mggGaussSig");
  RooAbsPdf* mggCBSigAll = w->pdf("mggCBSig");
  RooAbsPdf* mggSigAll = w->pdf("mggSig");
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
  RooPlot* plotmggAll = mgg->frame(Range(minSigFit,maxSigFit),Bins(nBinsMass));
  signalAll->plotOn(plotmggAll);
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
    sigToFit[c]->plotOn(plotmgg[c],LineColor(kWhite),MarkerColor(kWhite));
    mggSig[c] ->plotOn(plotmgg[c]);
    double chi2n = plotmgg[c]->chiSquare(0) ;
    cout << "------------------------- Experimentakl chi2 = " << chi2n << endl;
    mggSig[c] ->plotOn(
		       plotmgg[c],
		       Components(TString::Format("GaussSig_cat%d",c)),
		       LineStyle(kDashed),LineColor(kGreen));
    mggSig[c] ->plotOn(
		       plotmgg[c],
		       Components(TString::Format("CBSig_cat%d",c)),
		       LineStyle(kDashed),LineColor(kRed));
    mggSig[c] ->paramOn(plotmgg[c]);
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
    legmc->AddEntry(plotmgg[c]->getObject(5),"Simulation","LPE");
    legmc->AddEntry(plotmgg[c]->getObject(1),"Parametric Model","L");
    legmc->AddEntry(plotmgg[c]->getObject(2),"Gaussian Outliers","L");
    legmc->AddEntry(plotmgg[c]->getObject(3),"Crystal Ball component","L");
    legmc->SetHeader(" ");
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();
    // float effS = effSigma(hist);
    TLatex *lat = new TLatex(
			     minSigFit+0.5,0.85*plotmgg[c]->GetMaximum(),
			     " Resonance - 350 GeV");
    lat->Draw();
    TLatex *lat2 = new TLatex(
			      minSigFit+1.5,0.75*plotmgg[c]->GetMaximum(),catdesc.at(c));
    lat2->Draw();
    ///////
    char myChi2buffer[50];
    sprintf(myChi2buffer,"#chi^{2}/ndof = %f",chi2n);
    TLatex* latex = new TLatex(0.52, 0.7, myChi2buffer);
    latex -> SetNDC();
    latex -> SetTextFont(42);
    latex -> SetTextSize(0.04);
    //latex -> Draw("same");
    ctmp->SaveAs(TString::Format("sigmodel_cat%d.pdf",c));
    ctmp->SaveAs(TString::Format("sigmodel_cat%d.png",c));
    //ctmp->SaveAs(TString::Format("sigmodel_cat%d.C",c));
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
  TString component[5] = {"ggH","ttH","VBF","VH", "bbH"};

  //    higToFit[c] = (RooDataSet*) w->data(TString::Format("Hig_%d_cat%d",higgschannel,c));

  for (int d = 0; d < 5; ++d){

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
      sigToFit[c]->plotOn(plotmgg[c],LineColor(kWhite),MarkerColor(kWhite));
      mggSig[c] ->plotOn(plotmgg[c]);
      double chi2n = plotmgg[c]->chiSquare(0) ;
      cout << "------------------------- Experimentakl chi2 = " << chi2n << endl;
      mggSig[c] ->plotOn(
			 plotmgg[c],
			 Components(TString::Format("GaussHig_%d_cat%d",d,c)),
			 LineStyle(kDashed),LineColor(kGreen));
      mggSig[c] ->plotOn(
			 plotmgg[c],
			 Components(TString::Format("CBHig_%d_cat%d",d,c)),
			 LineStyle(kDashed),LineColor(kRed));
      mggSig[c] ->paramOn(plotmgg[c]);
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
      legmc->AddEntry(plotmgg[c]->getObject(2),"Gaussian Outliers","L");
      legmc->AddEntry(plotmgg[c]->getObject(3),"Crystal Ball component","L");
      legmc->SetHeader(" ");
      legmc->SetBorderSize(0);
      legmc->SetFillStyle(0);
      legmc->Draw();
      // float effS = effSigma(hist);
      TLatex *lat = new TLatex(
			       minSigFit+0.5,0.85*plotmgg[c]->GetMaximum(),
			       " Resonance - 350 GeV");
      lat->Draw();
      TLatex *lat2 = new TLatex(
				minSigFit+1.5,0.75*plotmgg[c]->GetMaximum(),catdesc.at(c));
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
  //higTree->AddVar(scaleWeightVar);
  // common preselection cut
  TString mainCut("1");
  // one channel with right weights
  //RooDataSet higScaled1(
  RooDataSet higScaled(
		       "higScaled1",
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
					       *w->var("mgg"),
					       mainCut+TString::Format(" && cut_based_ct==%d ",0)+cut0+cutj0);
  w->import(*higToFit[0],Rename(TString::Format("Hig_%d_cat%d",higgschannel,0)));
  //
  higToFit[1] = (RooDataSet*) higScaled.reduce(
					       *w->var("mgg"),
					       mainCut+TString::Format(" && cut_based_ct==%d ",1)+cut1+cutj1);
  w->import(*higToFit[1],Rename(TString::Format("Hig_%d_cat%d",higgschannel,1))); // Create full signal data set without categorization
  RooDataSet* higToFitAll = (RooDataSet*) higScaled->reduce(*w->var("mgg"),mainCut);
  w->import(*higToFitAll,Rename("Hig"));
  // here we print the number of entries on the different categories
  cout << "========= the number of entries on the different categories ==========" << endl;
  cout << "---- one channel: " << higScaled->sumEntries() << endl;
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
  wAll->factory(TString::Format("CMS_hgg_hig_%d_m0_absShift[1,1,1]",higgschannel));
  wAll->factory(TString::Format("prod::CMS_hgg_hig_m0_%d_cat0(mgg_hig_m0_%d_cat0, CMS_hgg_hig_%d_m0_absShift)",higgschannel,higgschannel,higgschannel));
  wAll->factory(TString::Format("prod::CMS_hgg_hig_m0_%d_cat1(mgg_hig_m0_%d_cat1, CMS_hgg_hig_%d_m0_absShift)",higgschannel,higgschannel,higgschannel));
  // (3) Systematics on resolution
  wAll->factory(TString::Format("CMS_hgg_hig_%d_sigmaScale[1,1,1]",higgschannel));
  wAll->factory(TString::Format("prod::CMS_hgg_hig_sigma_%d_cat0(mgg_hig_sigma_%d_cat0, CMS_hgg_hig_%d_sigmaScale)",higgschannel,higgschannel,higgschannel));

  wAll->factory(TString::Format("prod::CMS_hgg_hig_sigma_%d_cat1(mgg_hig_sigma_%d_cat1, CMS_hgg_hig_%d_sigmaScale)",higgschannel,higgschannel,higgschannel));
  wAll->factory(TString::Format("prod::CMS_hgg_hig_gsigma_%d_cat0(mgg_hig_gsigma_%d_cat0, CMS_hgg_hig_%d_sigmaScale)",higgschannel,higgschannel,higgschannel));
  wAll->factory(TString::Format("prod::CMS_hgg_hig_gsigma_%d_cat1(mgg_hig_gsigma_%d_cat1, CMS_hgg_hig_%d_sigmaScale)",higgschannel,higgschannel,higgschannel));
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
  TString filename(wsDir+TString(fileHiggsName)+".inputhig.root");
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
void MakeDataCard(RooWorkspace* w, const char* fileBaseName, const char* fileBkgName , const char* fileHiggsNameggh, const char* fileHiggsNametth, const char* fileHiggsNamevbf, const char* fileHiggsNamevh, const char* fileHiggsNamebbh) {
  TString cardDir = "datacards/";
  const Int_t ncat = NCAT;
  RooDataSet* data[ncat];
  RooDataSet* sigToFit[ncat];
  RooDataSet* higToFitggh[ncat];
  RooDataSet* higToFittth[ncat];
  RooDataSet* higToFitvbf[ncat];
  RooDataSet* higToFitvh[ncat];
  RooDataSet* higToFitbbh[ncat];
  for (int c = 0; c < ncat; ++c) {
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    sigToFit[c] = (RooDataSet*) w->data(TString::Format("Sig_cat%d",c));
    //
    higToFitggh[c] = (RooDataSet*) w->data(TString::Format("Hig_%d_cat%d",0,c));
    higToFittth[c] = (RooDataSet*) w->data(TString::Format("Hig_%d_cat%d",1,c));
    higToFitvbf[c] = (RooDataSet*) w->data(TString::Format("Hig_%d_cat%d",2,c));
    higToFitvh[c] = (RooDataSet*) w->data(TString::Format("Hig_%d_cat%d",3,c));
    higToFitbbh[c] = (RooDataSet*) w->data(TString::Format("Hig_%d_cat%d",4,c));
  } // close cat
  ////////////////////////////////////////////////////////////////////////////////////
  //RooRealVar* lumi = w->var("lumi");
  cout << "======== Expected Events Number =====================" << endl;
  cout << ".........Measured Data for L = " << "19712" << " pb-1 ............................" << endl;
  cout << "#Events data: " << w->data("Data")->sumEntries() << endl;
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events data cat%d: ",c) << data[c]->sumEntries() << endl;
  }
  cout << ".........Expected Signal for L = " << "19712" << " pb-1 ............................" << endl;
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
  outFile << "#Run with: combine -d hgg.mH350.0.shapes-Unbinned.txt -U -m 130 -H ProfileLikelihood -M MarkovChainMC --rMin=0 --rMax=20.0 -b 3500 -i 50000 --optimizeSim=1 --tries 30" << endl;
  outFile << "# Lumi = " << "19712" << " pb-1" << endl;
  outFile << "imax "<<ncat << endl;
  outFile << "jmax 6" << endl; // number of BKG
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
  outFile << "shapes mggHigggh cat0 " << TString(fileHiggsNameggh)+".inputhig.root" << " w_all:CMS_hgg_hig_0_cat0" << endl;
  outFile << "shapes mggHigggh cat1 " << TString(fileHiggsNameggh)+".inputhig.root" << " w_all:CMS_hgg_hig_0_cat1" << endl;
  outFile << "# tth" << endl;
  outFile << "shapes mggHigtth cat0 " << TString(fileHiggsNametth)+".inputhig.root" << " w_all:CMS_hgg_hig_1_cat0" << endl;
  outFile << "shapes mggHigtth cat1 " << TString(fileHiggsNametth)+".inputhig.root" << " w_all:CMS_hgg_hig_1_cat1" << endl;
  outFile << "# vbf" << endl;
  outFile << "shapes mggHigvbf cat0 " << TString(fileHiggsNamevbf)+".inputhig.root" << " w_all:CMS_hgg_hig_2_cat0" << endl;
  outFile << "shapes mggHigvbf cat1 " << TString(fileHiggsNamevbf)+".inputhig.root" << " w_all:CMS_hgg_hig_2_cat1" << endl;
  outFile << "# vh" << endl;
  outFile << "shapes mggHigvh cat0 " << TString(fileHiggsNamevh)+".inputhig.root" << " w_all:CMS_hgg_hig_3_cat0" << endl;
  outFile << "shapes mggHigvh cat1 " << TString(fileHiggsNamevh)+".inputhig.root" << " w_all:CMS_hgg_hig_3_cat1" << endl;
  outFile << "# bbh" << endl;
  outFile << "shapes mggHigbbh cat0 " << TString(fileHiggsNamebbh)+".inputhig.root" << " w_all:CMS_hgg_hig_4_cat0" << endl;
  outFile << "shapes mggHigbbh cat1 " << TString(fileHiggsNamebbh)+".inputhig.root" << " w_all:CMS_hgg_hig_4_cat1" << endl;
  outFile << "---------------" << endl;
  /////////////////////////////////////
  if(addHiggs) { //
    outFile << "bin cat0 cat1 " << endl;
    cout<<"here"<<endl;
    outFile << "observation "<< data[0]->sumEntries() <<" " << data[1]->sumEntries() <<" "<< endl;
    outFile << "------------------------------" << endl;
    outFile << "bin cat0 cat0 cat0 cat0 cat0 cat0 cat0"
	    <<" cat1 cat1 cat1 cat1 cat1 cat1 cat1" << endl;
    outFile << "process mggSig mggBkg mggHigggh mggHigtth mggHigvbf mggHigvh mggHigbbh"
	    <<" mggSig mggBkg mggHigggh mggHigtth mggHigvbf mggHigvh mggHigbbh" << endl;
    outFile << "process 0 1 2 3 4 5 6"
	    <<" 0 1 2 3 4 5 6 " << endl;
    outFile << "rate "
	    <<" "<<sigToFit[0]->sumEntries()<<" "<<1<<" "<<higToFitggh[0]->sumEntries()<<" "<<higToFittth[0]->sumEntries()<<" "<<higToFitvbf[0]->sumEntries()<<" "<<higToFitvh[0]->sumEntries()<< " "<<higToFitbbh[0]->sumEntries()
	    <<" "<<sigToFit[1]->sumEntries()<<" "<<1<<" "<<higToFitggh[1]->sumEntries()<<" "<<higToFittth[1]->sumEntries()<<" "<<higToFitvbf[1]->sumEntries()<<" "<<higToFitvh[1]->sumEntries()<<" "<<higToFitbbh[1]->sumEntries()
	    <<" "<<endl;
    outFile << " " << endl;
    outFile << "############## Total normalisation" << endl;
    outFile << "lumi_8TeV lnN "
	    << "1.026 - 1.026 1.026 1.026 1.026 1.026 "
	    << "1.026 - 1.026 1.026 1.026 1.026 1.026 " << endl;
    outFile << " " << endl;
    outFile << "############## Photon selection normalisation uncertainties " << endl;
    outFile << "DiphoTrigger lnN "
	    << "1.01 - 1.010 1.010 1.010 1.010 1.010 "
	    << "1.01 - 1.010 1.010 1.010 1.010 1.010 "
	    << "# Trigger efficiency" << endl;
    outFile << "CMS_hgg_eff_g lnN "
	    << "1.010 - 1.010 1.010 1.010 1.010 1.010 "
	    << "1.010 - 1.010 1.010 1.010 1.010 1.010 "
	    << "# photon selection accep." << endl;
    outFile << " " << endl;
    outFile << "############## Jet selection and phase space cuts normalisation uncertainties " << endl;
    outFile << "Mjj_PTj_cut_acceptance lnN "
	    << "1.015 - 1.015 1.015 1.015 1.015 1.015 "
	    << "1.015 - 1.015 1.015 1.015 1.015 1.015 "
	    <<"# JER and JES " << endl;
    outFile << "btag_eff lnN "
	    << "1.046 - 1.046 1.046 1.046 1.046 1.046 "
	    << "0.988 - 0.988 0.988 0.988 0.988 0.988 "
	    <<"# b tag efficiency uncertainty" << endl;
    outFile << "maajj_cut_acceptance lnN "
	    << "1.02 - 1.02 1.02 1.02 1.02 1.02 "
	    << "1.02 - 1.02 1.02 1.02 1.02 1.02 " << endl;
    outFile << " " << endl;
    outFile << "############## Theory uncertainties on SM Higgs production " << endl;
    outFile << "PDF lnN "
	    << " - - 0.931/1.075 0.919/1.081 0.972/1.026 0.976/1.024 0.976/1.024 "
	    << " - - 0.931/1.075 0.919/1.081 0.972/1.026 0.976/1.024 0.976/1.024 " << endl;
    outFile << "QCD_scale lnN "
	    << " - - 0.922/1.072 0.907/1.038 0.998/1.002 0.980/1.020 0.980/1.020 "
	    << " - - 0.922/1.072 0.907/1.038 0.998/1.002 0.980/1.020 0.980/1.020 " << endl;
    outFile << "gg_migration lnN "
	    << " - - 1.25 1.25 1.08 1.08 1.08 "
	    << " - - 1.25 1.25 1.08 1.08 1.08 # UEPS" << endl;
    outFile << "gluonSplitting lnN "
	    << " - - 1.40 1.40 1.40 1.40 1.40 "
	    << " - - 1.40 1.40 1.40 1.40 1.40 " << endl;
    outFile << " " << endl;
    outFile << "############## Signal parametric shape uncertainties " << endl;
    outFile << "CMS_hgg_sig_m0_absShift param 1 0.0057 # displacement of the dipho mean error = sqrt(0.45^ 2 + 0.35^ 2) " << endl;
    outFile << "CMS_hgg_sig_sigmaScale param 1 0.22 # optimistic estimative of resolution uncertainty " << endl;
    //
    outFile << "# Parametric shape uncertainties, entered by hand. they act on higgs" << endl;
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
    outFile << "CMS_hgg_hig_m0_4_absShift param 1 0.0057 # displacement of the dipho mean error = sqrt(0.45^ 2 + 0.35^ 2)" << endl;
    outFile << "CMS_hgg_hig_4_sigmaScale param 1 0.22 # optimistic estimative of resolution uncertainty " << endl;
    //
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




void MakeDataCardonecatnohiggs(RooWorkspace* w, const char* fileBaseName, const char* fileBkgName) {
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
  cout << ".........Measured Data for L = " << "19712" << " pb-1 ............................" << endl;
  cout << "#Events data: " << w->data("Data")->sumEntries() << endl;
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events data cat%d: ",c) << data[c]->sumEntries() << endl;
  }
  // cout << ".........Expected Signal for L = " << lumi->getVal() << " pb-1 ............................" << endl;
  cout << ".........Expected Signal for L = " << "19712" << " pb-1 ............................" << endl;
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
  outFile << "#Run with: combine -d hgg.mH350.0.shapes-Unbinned.txt -U -m 130 -H ProfileLikelihood -M MarkovChainMC --rMin=0 --rMax=20.0 -b 3500 -i 50000 --optimizeSim=1 --tries 30" << endl;
  // outFile << "# Lumi = " << lumi->getVal() << " pb-1" << endl;
  outFile << "# Lumi = " << "19712" << " pb-1" << endl;
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
  defaultStyle->SetOptStat(0000);
  defaultStyle->SetOptFit(000);
  defaultStyle->SetPalette(1);
  /////// pad ////////////
  defaultStyle->SetPadBorderMode(1);
  defaultStyle->SetPadBorderSize(1);
  defaultStyle->SetPadColor(0);
  defaultStyle->SetPadTopMargin(0.05);
  defaultStyle->SetPadBottomMargin(0.13);
  defaultStyle->SetPadLeftMargin(0.13);
  defaultStyle->SetPadRightMargin(0.02);
  /////// canvas /////////
  defaultStyle->SetCanvasBorderMode(0);
  defaultStyle->SetCanvasColor(0);
  defaultStyle->SetCanvasDefH(600);
  defaultStyle->SetCanvasDefW(600);
  /////// frame //////////
  defaultStyle->SetFrameBorderMode(0);
  defaultStyle->SetFrameBorderSize(1);
  defaultStyle->SetFrameFillColor(0);
  defaultStyle->SetFrameLineColor(1);
  /////// label //////////
  defaultStyle->SetLabelOffset(0.005,"XY");
  defaultStyle->SetLabelSize(0.05,"XY");
  defaultStyle->SetLabelFont(42,"XY");
  /////// title //////////
  defaultStyle->SetTitleOffset(1.1,"X");
  defaultStyle->SetTitleSize(0.01,"X");
  defaultStyle->SetTitleOffset(1.25,"Y");
  defaultStyle->SetTitleSize(0.05,"Y");
  defaultStyle->SetTitleFont(42, "XYZ");
  /////// various ////////
  defaultStyle->SetNdivisions(505,"Y");
  defaultStyle->SetLegendBorderSize(0); // For the axis titles:

  defaultStyle->SetTitleColor(1, "XYZ");
  defaultStyle->SetTitleFont(42, "XYZ");
  defaultStyle->SetTitleSize(0.06, "XYZ");
 
  // defaultStyle->SetTitleYSize(Float_t size = 0.02);
  defaultStyle->SetTitleXOffset(0.9);
  defaultStyle->SetTitleYOffset(1.05);
  // defaultStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

  // For the axis labels:
  defaultStyle->SetLabelColor(1, "XYZ");
  defaultStyle->SetLabelFont(42, "XYZ");
  defaultStyle->SetLabelOffset(0.007, "XYZ");
  defaultStyle->SetLabelSize(0.04, "XYZ");

  // For the axis:
  defaultStyle->SetAxisColor(1, "XYZ");
  defaultStyle->SetStripDecimals(kTRUE);
  defaultStyle->SetTickLength(0.03, "XYZ");
  defaultStyle->SetNdivisions(510, "XYZ");
  defaultStyle->SetPadTickX(1); // To get tick marks on the opposite side of the frame
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

