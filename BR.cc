
//#include <iostream>
//#include <fstream>
//#include <string>
//using namespace std;
//int main
{
   int ncat=3;
  // comparing WP's
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
  defaultStyle->SetPadLeftMargin(0.15);
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
  defaultStyle->SetLegendBorderSize(0);  // For the axis titles:

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
    defaultStyle->SetPadTickX(1);   // To get tick marks on the opposite side of the frame
    defaultStyle->SetPadTickY(1);
    defaultStyle->cd();

 bool HH=false;
 bool base = true;
 bool low=false;
 bool obs=false;
 bool oldExp= true;

//  TLegend *leg = new TLegend(0.65,0.5,0.99,0.94);
  TLegend *leg = new TLegend(0.55,0.5,0.99,0.94);
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  double radCX[16]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
  double radMASS[16]={260,270,300,350,400,450,500,550,600,650,700,800,900,1000,1100,1200};
  if(!HH) double br=1;//1./(0.577*0.00228); 
  if(HH) double br=1./(2*0.577*0.00228); // HH
//////////////////////////////////////////
// draw the radion line // MR      rad_CX(fb)       grav_CX(fb)
   //ntuple = new TNtuple("ntuple","NTUPLE","x:y:z");
   FILE *fp = fopen("CX_runnig_aS.data","r");
   ntuple = new TNtuple("ntuple","NTUPLE","x:y"); 
   float x0,y0;
   char line[127];
   while (fgets(&line,127,fp)) {
      sscanf(&line[0],"%f %f",&x0,&y0);
      printf("x0=%f, y0=%f\n",x0,y0);
      ntuple->Fill(x0,y0); 
      //cout<<x0<<endl;
      //ntuple->Fill(x0,y0); 
   }
   cout<<ntuple.GetNvar()<<" "<<ntuple.GetVar1()<<ntuple.GetVar2()<<endl;
  ////////////////////////////////////////////////////////////////////////
  // graviton
  cout<<"Bulk graviton"<<endl;
  FILE *fpg = fopen("CX_grav_asfeaseability.data","r");
  ntupleg = new TNtuple("ntupleg","NTUPLE","x:y"); 
  float x0,y0;
  char line[127];
  while (fgets(&line,127,fpg)) {
      sscanf(&line[0],"%f %f",&x0,&y0);
      printf("x0=%f, y0=%f\n",x0,y0);
      ntupleg->Fill(x0,y0); 
      //cout<<x0<<endl;
      //ntuple->Fill(x0,y0); 
  }
  cout<<ntupleg.GetNvar()<<" "<<ntupleg.GetVar1()<<ntupleg.GetVar2()<<endl;
///////////////////////////////////////////////
// souvick NMSSM
// mh2[GeV] $t_\beta=1.5$ $t_\beta$=2 $t_\beta$=2.5& $t_\beta$=3 $t_\beta$=3.5 $t_\beta$ =4 \\ \hline


//return 0;
//   TGraph *radion = new TGraph(ntuple->GetSelectedRows(),
//                                   ntuple->GetV2(), ntuple->GetV1());
///////////////////////////////////////
// draw the BR flag -- Main results
TFile* f1 = new TFile("/afs/cern.ch/work/b/bmarzocc/LimitTrees/CMSSW_6_1_1/src/Limit_codes/radlim_CSV_WP4/radlim260_CSV/higgsCombineTest.Asymptotic.mH125.mR260_higgs.root");
TFile* f2 = new TFile("/afs/cern.ch/work/b/bmarzocc/LimitTrees/CMSSW_6_1_1/src/Limit_codes/radlim_CSV_WP4/radlim270_CSV/higgsCombineTest.Asymptotic.mH125.mR270_higgs.root");
TFile* f3 = new TFile("/afs/cern.ch/work/b/bmarzocc/LimitTrees/CMSSW_6_1_1/src/Limit_codes/radlim_CSV_WP4/radlim300_CSV/higgsCombineTest.Asymptotic.mH125.mR300_higgs.root");
TFile* f4 = new TFile("/afs/cern.ch/work/b/bmarzocc/LimitTrees/CMSSW_6_1_1/src/Limit_codes/radlim_CSV_WP4/radlim350_CSV/higgsCombineTest.Asymptotic.mH125.mR350_higgs.root");
TFile* f5 = new TFile("/afs/cern.ch/work/b/bmarzocc/LimitTrees/CMSSW_6_1_1/src/Limit_codes/radlim_CSV_WP4/radlim400_CSV/higgsCombineTest.Asymptotic.mH125.mR400.root");
TFile* f6 = new TFile("/afs/cern.ch/work/b/bmarzocc/LimitTrees/CMSSW_6_1_1/src/Limit_codes/radlim_CSV_WP4/radlim450_CSV/higgsCombineTest.Asymptotic.mH125.mR450.root");
TFile* f7 = new TFile("/afs/cern.ch/work/b/bmarzocc/LimitTrees/CMSSW_6_1_1/src/Limit_codes/radlim_CSV_WP4/radlim500_CSV/higgsCombineTest.Asymptotic.mH125.mR500.root");
 TFile* f8 = new TFile("/afs/cern.ch/work/b/bmarzocc/LimitTrees/CMSSW_6_1_1/src/Limit_codes/radlim_CSV_WP4/radlim550_CSV/higgsCombineTest.Asymptotic.mH125.mR550.root");
 TFile* f9 = new TFile("/afs/cern.ch/work/b/bmarzocc/LimitTrees/CMSSW_6_1_1/src/Limit_codes/radlim_CSV_WP4/radlim600_CSV/higgsCombineTest.Asymptotic.mH125.mR600.root");
 TFile* f10 = new TFile("/afs/cern.ch/work/b/bmarzocc/LimitTrees/CMSSW_6_1_1/src/Limit_codes/radlim_CSV_WP4/radlim650_CSV/higgsCombineTest.Asymptotic.mH125.mR650.root");
 TFile* f11 = new TFile("/afs/cern.ch/work/b/bmarzocc/LimitTrees/CMSSW_6_1_1/src/Limit_codes/radlim_CSV_WP4/radlim700_CSV/higgsCombineTest.Asymptotic.mH125.mR700.root");
 TFile* f12 = new TFile("/afs/cern.ch/work/b/bmarzocc/LimitTrees/CMSSW_6_1_1/src/Limit_codes/radlim_CSV_WP4/radlim800_CSV/higgsCombineTest.Asymptotic.mH125.mR800.root");
 TFile* f13 = new TFile("/afs/cern.ch/work/b/bmarzocc/LimitTrees/CMSSW_6_1_1/src/Limit_codes/radlim_CSV_WP4/radlim900_CSV/higgsCombineTest.Asymptotic.mH125.mR900.root");
 TFile* f14 = new TFile("/afs/cern.ch/work/b/bmarzocc/LimitTrees/CMSSW_6_1_1/src/Limit_codes/radlim_CSV_WP4/radlim1000_CSV/higgsCombineTest.Asymptotic.mH125.mR1000.root");
 TFile* f15 = new TFile("/afs/cern.ch/work/b/bmarzocc/LimitTrees/CMSSW_6_1_1/src/Limit_codes/radlim_CSV_WP4/radlim1100_CSV/higgsCombineTest.Asymptotic.mH125.mR1100.root");
 TFile* f16 = new TFile("/afs/cern.ch/work/b/bmarzocc/LimitTrees/CMSSW_6_1_1/src/Limit_codes/radlim_CSV_WP4/radlim1100_CSV/higgsCombineTest.Asymptotic.mH125.mR1100.root");
// TFile* f5 = new TFile("radlim1500_CSV/higgsCombineTest.Asymptotic.mR1500_CSV.root");

 TTree          *fChain1 = (TTree*) f1->Get("limit;1");
 TTree          *fChain2 = (TTree*) f2->Get("limit;1");
 TTree          *fChain3 = (TTree*) f3->Get("limit;1");
 TTree          *fChain4 = (TTree*) f4->Get("limit;1");
 TTree          *fChain5 = (TTree*) f5->Get("limit;1");
 TTree          *fChain6 = (TTree*) f6->Get("limit;1");
 TTree          *fChain7 = (TTree*) f7->Get("limit;1");
 TTree          *fChain8 = (TTree*) f8->Get("limit;1");
 TTree          *fChain9 = (TTree*) f9->Get("limit;1");
 TTree          *fChain10 = (TTree*) f10->Get("limit;1");
 TTree          *fChain11 = (TTree*) f11->Get("limit;1");
 TTree          *fChain12 = (TTree*) f12->Get("limit;1");
 TTree          *fChain13 = (TTree*) f13->Get("limit;1");
 TTree          *fChain14 = (TTree*) f14->Get("limit;1");
 TTree          *fChain15 = (TTree*) f15->Get("limit;1");
 TTree          *fChain16 = (TTree*) f16->Get("limit;1");
// TTree          *fChain5 = (TTree*) f5->Get("limit;1");

  fChain1->SetMakeClass(1);
  fChain2->SetMakeClass(1);
  fChain3->SetMakeClass(1);
  fChain4->SetMakeClass(1);
  fChain5->SetMakeClass(1);
  fChain6->SetMakeClass(1);
  fChain7->SetMakeClass(1);
  fChain8->SetMakeClass(1);
  fChain9->SetMakeClass(1);
  fChain10->SetMakeClass(1);
  fChain11->SetMakeClass(1);
  fChain12->SetMakeClass(1);
  fChain13->SetMakeClass(1);
  fChain14->SetMakeClass(1);
  fChain15->SetMakeClass(1);
  fChain16->SetMakeClass(1);
//  fChain5->SetMakeClass(1);

  TBranch        *b_limit; 
  Double_t         limit;


 fChain1->SetBranchAddress("limit", &limit, &b_limit);
 fChain2->SetBranchAddress("limit", &limit, &b_limit); 
 fChain3->SetBranchAddress("limit", &limit, &b_limit);
 fChain4->SetBranchAddress("limit", &limit, &b_limit); 
 fChain5->SetBranchAddress("limit", &limit, &b_limit);
 fChain6->SetBranchAddress("limit", &limit, &b_limit); 
 fChain7->SetBranchAddress("limit", &limit, &b_limit);
 fChain8->SetBranchAddress("limit", &limit, &b_limit);
 fChain9->SetBranchAddress("limit", &limit, &b_limit);
 fChain10->SetBranchAddress("limit", &limit, &b_limit); 
 fChain11->SetBranchAddress("limit", &limit, &b_limit);
 fChain12->SetBranchAddress("limit", &limit, &b_limit);
 fChain13->SetBranchAddress("limit", &limit, &b_limit);
 fChain14->SetBranchAddress("limit", &limit, &b_limit); 
 fChain15->SetBranchAddress("limit", &limit, &b_limit);
 fChain16->SetBranchAddress("limit", &limit, &b_limit);
// fChain5->SetBranchAddress("limit", &limit, &b_limit);
 int nmass=15;
 float rad[6][nmass];
 float radOld[6][20];

// cout << "Radion270 " <<endl;
  int j=0; 
  for (int i = 0; i<6; i++){
    fChain1->GetTree()->GetEntry(i);
    rad[i][j]=limit*radCX[j]*br;
    cout <<"MX = "<< radMASS[j]<<" centrality "<< i<<" limit = " << limit << endl; 
    }
  int j=1; 
  for (int i = 0; i<6; i++){
    fChain2->GetTree()->GetEntry(i);
    rad[i][j]=limit*radCX[j]*br;
    cout <<"MX = "<< radMASS[j]<<" centrality "<< i<<" limit = " << limit << endl; 
    }
  int j=2; 
  for (int i = 0; i<6; i++){
    fChain3->GetTree()->GetEntry(i);
    rad[i][j]=limit*radCX[j]*br;
    cout <<"MX = "<< radMASS[j]<<" centrality "<< i<<" limit = " << limit << endl; 
    }
  int j=3; 
  for (int i = 0; i<6; i++){
    fChain4->GetTree()->GetEntry(i);
    rad[i][j]=limit*radCX[j]*br;
    cout <<"MX = "<< radMASS[j]<<" centrality "<< i<<" limit = " << limit << endl; 
    }
  int j=4; 
  for (int i = 0; i<6; i++){
    fChain5->GetTree()->GetEntry(i);
    rad[i][j]=limit*radCX[j]*br;
    cout <<"MX = "<< radMASS[j]<<" centrality "<< i<<" limit = " << limit << endl; 
    }
  int j=5; 
  for (int i = 0; i<6; i++){
    fChain6->GetTree()->GetEntry(i);
    rad[i][j]=limit*radCX[j]*br;
    cout <<"MX = "<< radMASS[j]<<" centrality "<< i<<" limit = " << limit << endl; 
    }
  int j=6; 
  for (int i = 0; i<6; i++){
    fChain7->GetTree()->GetEntry(i);
    rad[i][j]=limit*radCX[j]*br;
    cout <<"MX = "<< radMASS[j]<<" centrality "<< i<<" limit = " << limit << endl; 
    }
  int j=7; 
  for (int i = 0; i<6; i++){
    fChain8->GetTree()->GetEntry(i);
    rad[i][j]=limit*radCX[j]*br;
    cout <<"MX = "<< radMASS[j]<<" centrality "<< i<<" limit = " << limit << endl; 
    }
  int j=8; 
  for (int i = 0; i<6; i++){
    fChain9->GetTree()->GetEntry(i);
    rad[i][j]=limit*radCX[j]*br;
    cout <<"MX = "<< radMASS[j]<<" centrality "<< i<<" limit = " << limit << endl; 
    }
  int j=9; 
  for (int i = 0; i<6; i++){
    fChain10->GetTree()->GetEntry(i);
    rad[i][j]=limit*radCX[j]*br;
    cout <<"MX = "<< radMASS[j]<<" centrality "<< i<<" limit = " << limit << endl; 
    }
  int j=10; 
  for (int i = 0; i<6; i++){
    fChain11->GetTree()->GetEntry(i);
    rad[i][j]=limit*radCX[j]*br;
    cout <<"MX = "<< radMASS[j]<<" centrality "<< i<<" limit = " << limit << endl; 
    }
  int j=11; 
  for (int i = 0; i<6; i++){
    fChain12->GetTree()->GetEntry(i);
    rad[i][j]=limit*radCX[j]*br;
    cout <<"MX = "<< radMASS[j]<<" centrality "<< i<<" limit = " << limit << endl; 
    }
  int j=12; 
  for (int i = 0; i<6; i++){
    fChain13->GetTree()->GetEntry(i);
    rad[i][j]=limit*radCX[j]*br;
    cout <<"MX = "<< radMASS[j]<<" centrality "<< i<<" limit = " << limit << endl; 
    }
  int j=13; 
  for (int i = 0; i<6; i++){
    fChain14->GetTree()->GetEntry(i);
    rad[i][j]=limit*radCX[j]*br;
    cout <<"MX = "<< radMASS[j]<<" centrality "<< i<<" limit = " << limit << endl; 
    }
  int j=14; 
  for (int i = 0; i<6; i++){
    fChain15->GetTree()->GetEntry(i);
    rad[i][j]=limit*radCX[j]*br;
    cout <<"MX = "<< radMASS[j]<<" centrality "<< i<<" limit = " << limit << endl; 
    }
  
  //PAS expected results
  radOld[2][0] = 2.12;
  radOld[2][1] = 2.40;
  radOld[2][2] = 2.73;
  radOld[2][3] = 2.23;
  radOld[2][4] = 1.87;
  radOld[2][5] = 1.42;
  radOld[2][6] = 0.97;
  radOld[2][7] = 0.80;
  radOld[2][8] = 0.69;
  radOld[2][9] = 0.60;
  radOld[2][10] = 0.54;
  radOld[2][11] = 0.46;
  radOld[2][12] = 0.43;
  radOld[2][13] = 0.43;
  radOld[2][14] = 0.48;

  
  
  // we do a plot r*MR
  TMultiGraph *mg = new TMultiGraph();
  if(!HH) mg->SetMinimum(0.01);
  if(HH) mg->SetMinimum(10); // HH
  if(HH && low) mg->SetMaximum(3500); // HH 1000000
  if(!HH && low) mg->SetMaximum(11); // 10000
  if(HH && !low) mg->SetMaximum(100000); // HH 
  if(!HH && !low) mg->SetMaximum(1000); // 
  if(!HH) mg->SetTitle("#splitline{#scale[1.0]{CX*BR(HH)*2*BR(#gamma #gamma)*BR(bb) (fb)}}{#scale[0.8]{CMS preliminary 19.7/fb  }}");
  if(HH) mg->SetTitle("#splitline{#scale[1.0]{CX*BR(HH)(fb) - SM Higgs BR}}{#scale[0.8]{CMS preliminary 19.7/fb  }}");
  //if(HH) mg->SetTitle("Radion > HH > #gamma #gamma bb~");
//  float radMASS[nmass]={300,500,700,1000};

   TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);
   Double_t x[nmass], yobs[nmass], y2up[nmass], y1up[nmass], y1down[nmass], y2down[nmass], ymean[nmass],ymeanOld[nmass];
   //
   int nmax;
   if(low) nmax= 4;//nmass;//nmass
   if(!low) nmax= nmass;//nmass
   c1.cd();
   if(low) c1->SetLogy(0);
   if(!low) c1->SetLogy(1);
   c1->SetGrid();
  if(HH) ntuple->Draw("y*(9/1)*0.25 : x");
  if(!HH)  ntuple->Draw("y*0.25*0.577*0.00228*2*9 : x");
   TGraphErrors *radion10 = 
	new TGraphErrors(ntuple->GetSelectedRows(), ntuple->GetV2(), ntuple->GetV1());
   radion10->SetLineColor(kRed); 
   radion10->SetLineWidth(3);
  if(HH)  ntuple->Draw("y*0.25 : x");
  if(!HH)  ntuple->Draw("y*0.25*0.577*0.00228*2 : x");
   TGraphErrors *radion = 
	new TGraphErrors(ntuple->GetSelectedRows(), ntuple->GetV2(), ntuple->GetV1());//point,radionx,radiony);
   radion->SetLineWidth(3);
   if(base){
   // RS graviton lambda =1
//   if(HH)  ntupleg->Draw("y : x");
//   if(!HH)  ntupleg->Draw("y*0.577*0.00228*2 : x");
   TGraphErrors *graviton3 = new TGraphErrors(1); 
// 1.6580e+01
   if(HH && base)
   for(int j=0;j<100;j++) graviton3->SetPoint(j, 250*(1+j), (3000*(0.2*3.83/(250*(1+j))))*(3000*(0.2*3.83/(250*(1+j))))*16.58);
   if(!HH && base)
   for(int j=0;j<100;j++) graviton3->SetPoint(j, 250*(1+j), (3000*(0.2*3.83/(250*(1+j))))*(3000*(0.2*3.83/(250*(1+j))))*16.58*0.577*0.00228*2 );
//	new TGraphErrors(ntupleg->GetSelectedRows(), ntupleg->GetV2(), ntupleg->GetV1());
   graviton3->SetLineColor(kBlue); 
   graviton3->SetLineWidth(3);
   graviton3->SetLineStyle(2);
   // bulk graviton
   if(HH)  ntupleg->Draw("y : x");
   if(!HH)  ntupleg->Draw("y*0.577*0.00228*2 : x");
   TGraphErrors *bulk3 = 
	new TGraphErrors(ntupleg->GetSelectedRows(), ntupleg->GetV2(), ntupleg->GetV1());//point,radionx,radiony);
   bulk3->SetLineColor(kBlack); 
   bulk3->SetLineWidth(3);
   bulk3->SetLineStyle(2);
   } //if base plt graviton
   c1->Clear();
//   mg->Draw("AP");
//return 0;
   // draw radion CX
//  int point = ntuple->GetEntri();
  //float radionx[point]=ntuple->GetVar1();
  //float radiony[point]=ntuple->GetVar2();
  //cout<<radionx[1]<<endl;
//   ntuple->SetLineColor(kRed);
//   ntuple->SetLineWidth(3);
//   
   for (Int_t i=0;i<nmax;i++) {
//     x[i] = radMASS[i];
     yobs[i] = rad[5][i];
     y2up[i] = rad[0][i];             
     y1up[i] = rad[1][i];  
     ymean[i] = rad[2][i]; 
     ymeanOld[i] = radOld[2][i]; 
     y1down[i] = rad[3][i];  
     y2down[i] =  rad[4][i];  
   }

  TGraphErrors *grobs = new TGraphErrors(1);
  grobs->SetMarkerStyle(kFullDotLarge); 
  grobs->SetMarkerColor(kBlue); 
  grobs->SetLineColor(kBlue);
  grobs->SetLineWidth(3);
//  grobs->SetLineStyle(3);
  TGraphErrors *gr2up = new TGraphErrors(1);
//  gr2up->SetLineColor(1);
//  gr2up->SetLineWidth(1);
  TGraphErrors *gr1up = new TGraphErrors(1);
//  gr1up->SetLineColor(kRed);
//  gr1up->SetLineWidth(1);
  TGraphErrors *grmean = new TGraphErrors(1);
  grmean->SetLineColor(1);
  grmean->SetLineWidth(3);
  grmean->SetLineStyle(3);

  TGraphErrors *grmeanOld = new TGraphErrors(1);
  grmeanOld->SetLineColor(2);
  grmeanOld->SetLineWidth(3);
  grmeanOld->SetLineStyle(3);

  TGraphErrors *gr1down = new TGraphErrors(1);
//  gr1down->SetLineColor(kYellow);
//  gr1down->SetLineWidth(1);
  TGraphErrors *gr2down = new TGraphErrors(1);
//  gr2down->SetLineColor(kGreen);
//  gr2down->SetLineWidth(1);
 
  
  for(int j=0;j<nmax;j++){
    grobs->SetPoint(j, radMASS[j], yobs[j]);
    gr2up->SetPoint(j, radMASS[j], y2up[j]);
    gr1up->SetPoint(j, radMASS[j], y1up[j]);
    grmean->SetPoint(j, radMASS[j], ymean[j]);
    grmeanOld->SetPoint(j, radMASS[j], ymeanOld[j]);
    gr1down->SetPoint(j, radMASS[j], y1down[j]);    
    gr2down->SetPoint(j, radMASS[j], y2down[j]);
  }


   mg->Add(gr2up);//->Draw("same");
   mg->Add(gr1up);//->Draw("same");
   mg->Add(grmean,"L*");//->Draw("same,AC*");
   mg->Add(grmeanOld,"L*");//->Draw("same,AC*");
   mg->Add(gr1down);//->Draw("same,AC*");
   mg->Add(gr2down);//->Draw("same,AC*");
   if(obs) mg->Add(grobs,"L,P");//->Draw("AC*");



  mg->Draw("AP");
  mg->GetXaxis()->SetRangeUser(100,1400);
//   mg->GetYaxis()->SetTitle("CX 95% excluded");
  mg->GetXaxis()->SetTitle("M(GeV)");
  mg->GetYaxis()->SetTitleSize(0.05);
  mg->GetXaxis()->SetTitleSize(0.042);
  mg->GetYaxis()->SetLabelSize(0.05);
  mg->GetYaxis()->SetTitleOffset(1.05);
  mg->GetYaxis()->CenterTitle(true);
  mg->GetXaxis()->SetTitleOffset(1.01);
  mg->GetXaxis()->CenterTitle(true);
//  mg->Draw("AP,same");
//  mg->GetXaxis()->SetMoreLogLabels(true);
 //TLatex *lat  = new TLatex(90,90,"#splitline{#scale[1.0]{CMS Preliminary}}{#scale[0.8]{#sqrt{s} = 8 TeV L = 5.3 fb^{-1}}}");
    //lat->Draw("same");
  // histo to shade
   int n=nmax;

   TGraph *grgreen = new TGraph(2*n);
   TGraph *gryellow = new TGraph(2*n);
   for (i=0;i<n;i++) {
      grgreen->SetPoint(i,radMASS[i],y2up[i]);
      grgreen->SetPoint(n+i,radMASS[n-i-1],y2down[n-i-1]);
      //
      gryellow->SetPoint(i,radMASS[i],y1up[i]);
      gryellow->SetPoint(n+i,radMASS[n-i-1],y1down[n-i-1]);

      cout<<" observed "<<radMASS[i]<<" "<<y2down[i]<<endl; 
   }
//   gryellow->SetFillStyle(3013);
   grgreen->SetFillColor(kYellow);
   grgreen->Draw("f"); 
//   gryellow->SetFillStyle(3013);
   gryellow->SetFillColor(kGreen);
   gryellow->Draw("f"); 
  grmean->Draw("L,same");
  if(oldExp) grmeanOld->Draw("L,same");
  if(obs) grobs->Draw("L,P,E,same");
  
  /////////////////////////////////////////////////////////////
  // cross checks
  if(!HH){
  int check1mass[3] = {400,450,500};
  Double_t check1limit[3]={2.0078,1.6641,1.3008};
  TGraphErrors *gcheck1 = new TGraphErrors(1);
  gcheck1->SetMarkerStyle(kFullDotLarge); 
  gcheck1->SetLineColor(kRed+2);
  gcheck1->SetMarkerColor(kRed+2);
  gcheck1->SetLineWidth(3);
//  gcheck1->SetLineStyle(3);
  for(int j=0;j<3;j++) gcheck1->SetPoint(j, check1mass[j], check1limit[j]*br);
//  if(!low) if(!base)gcheck1->Draw("L,P");
  // sidebands
  int check2mass[5] = {270,300,350,400,450};
  Double_t check2limit[5]={2.5391,2.6172,1.7891,1.3750,1.0898};
  TGraphErrors *gcheck2 = new TGraphErrors(1);
  gcheck2->SetMarkerStyle(kFullDotLarge); 
  gcheck2->SetLineColor(kRed-7);
  gcheck2->SetMarkerColor(kRed-7);
  gcheck2->SetLineWidth(3);
  for(int j=0;j<5;j++) gcheck2->SetPoint(j, check2mass[j], check2limit[j]*br);
//  if(!base)gcheck2->Draw("L,P");
  //
  int check4mass[10] = {400,450,500,550,600,650,700,900,1000,1100};
  Double_t check4limit[10]={1.73 , 1.96 , 1.48 , 1.24 , 1.00 , 0.96 , 0.87 , 0.77 , 0.82 , 1.01};
  TGraphErrors *gcheck4 = new TGraphErrors(1);
  gcheck4->SetMarkerStyle(kFullDotLarge); 
  gcheck4->SetLineColor(kViolet);
  gcheck4->SetMarkerColor(kViolet);
  gcheck4->SetLineWidth(3);
//  gcheck1->SetLineStyle(3);
  for(int j=0;j<10;j++) gcheck4->SetPoint(j, check4mass[j], check4limit[j]*br);
  if(!low)if(!base) gcheck4->Draw("L,P");
  //
  //
  int check7mass[10] = {400,450,500,550,600,650,700,900,1000,1100};
  Double_t check7limit[10]={2.95 , 2.09 , 1.51 , 1.21 , 0.91 , 0.85 , 0.73 , 0.50 , 0.48 , 0.49 };
  TGraphErrors *gcheck7 = new TGraphErrors(1);
  gcheck7->SetMarkerStyle(kFullDotLarge); 
  gcheck7->SetLineColor(kRed+3);
  gcheck7->SetMarkerColor(kRed+3);
  gcheck7->SetLineWidth(3);
//  gcheck1->SetLineStyle(3);
  for(int j=0;j<10;j++) gcheck7->SetPoint(j, check7mass[j], check7limit[j]*br);
  if(!low)if(!base) gcheck7->Draw("L,P");
  } // close cross checks
  //
  //////////////////////////////////////////////////////////////////////
  // 2btag w higgs
  int check3mass[6] = {270,300,350,400,450,500};
  Double_t check3limit[6]={3.20 , 3.67 , 2.78 , 2.18 , 1.79 , 1.41};
  TGraphErrors *gcheck3 = new TGraphErrors(1);
  gcheck3->SetMarkerStyle(kFullDotLarge); 
  gcheck3->SetLineColor(kCyan);
  gcheck3->SetMarkerColor(kCyan);
  gcheck3->SetLineWidth(3);
//  gcheck1->SetLineStyle(3);
  for(int j=0;j<3;j++) gcheck3->SetPoint(j, check3mass[j], check3limit[j]*br);
//  if(low) if(!base) gcheck3->Draw("L,P");
  // 2btag wo higgs
  int check5mass[6] = {270,300,350,400,450,500};
  Double_t check5limit[6]={3.1406,3.6094,2.7266,2.3359,1.7734,1.3945};
  TGraphErrors *gcheck5 = new TGraphErrors(1);
  gcheck5->SetMarkerStyle(kFullDotLarge); 
  gcheck5->SetLineColor(kMagenta);
  gcheck5->SetMarkerColor(kMagenta);
  gcheck5->SetLineWidth(3);
//  gcheck1->SetLineStyle(3);
  for(int j=0;j<6;j++) gcheck5->SetPoint(j, check5mass[j], check5limit[j]*br);
  if(!low) if(!base) gcheck5->Draw("L,P");
  // wo higgs
  int check6mass[6] = {270,300,350,400,450,500};
  Double_t check6limit[6]={ 2.85 , 3.27 , 2.37 , 1.87 , 1.56 , 1.20};
  TGraphErrors *gcheck6 = new TGraphErrors(1);
  gcheck6->SetMarkerStyle(kFullDotLarge); 
  gcheck6->SetLineColor(kBlue);
  gcheck6->SetMarkerColor(kBlue);
  gcheck6->SetLineWidth(3);
//  gcheck1->SetLineStyle(3);
  for(int j=0;j<3;j++) gcheck6->SetPoint(j, check6mass[j], check6limit[j]*br);
//  if(low) if(!base) gcheck6->Draw("L,P");
  ///////////////////////////////////////////
if(!oldExp) radion->Draw("L,same");
if(!low && !oldExp)  radion10->Draw("L,same");
if(base && !oldExp)  graviton3->Draw("L,same");
if(base && !oldExp)  bulk3->Draw("L,same");
//   mg->Add(radion,"L,P");//->Draw("AC*");
//   mg->Add(radion10,"L,P");//->Draw("AC*");
if(base)  leg->SetHeader("WED: kl = 35, k/Mpl = 0.2, elementary top");
  
  if(!base && low) {leg->AddEntry(grmean, "Expected (m#gamma #gamma fit wo/ Higgs)", "L"); }
  else leg->AddEntry(grmean, "Expected v35", "L");
  
  if(oldExp) leg->AddEntry(grmeanOld, "Expected PAS", "L");

  if(obs) leg->AddEntry(grobs, "Observed", "L,P");
  leg->AddEntry(grgreen, "2 sigma", "f");
  leg->AddEntry(gryellow, "1 sigma", "f");
//  leg->AddEntry((TObject*)0, "WED: kl = 35, k/Mpl = 0.2, elementary top", "");
if(!oldExp) leg->AddEntry(radion, "RS radion (#LambdaR = 3 TeV)", "L");
if(!low && !oldExp)  leg->AddEntry(radion10, "RS radion (#LambdaR = 1 TeV)", "L");
if(base && !oldExp)  leg->AddEntry(graviton3, "\"RS1\" Graviton (DY + gluon fusion) ", "L");
if(base && !oldExp)  leg->AddEntry(bulk3, " Bulk Graviton *", "L");
//  leg->AddEntry((TObject*)0, "* k/Mpl = 0.2, elementary top", "");

  if(!base){
//  if(!low) leg->AddEntry(gcheck1, "m#gamma #gamma  fit", "LP");
  if(!low)leg->AddEntry(gcheck4, " 4 body fit - 2 btag only (w/ kinfit)", "LP");
  if(!low)leg->AddEntry(gcheck7, " 4 body fit - wo/ kinfit ", "LP");
  if(low)leg->AddEntry(gcheck5, "m#gamma #gamma  fit - 2 btag only (wo/ Higgs)", "LP");
//  if(low)leg->AddEntry(gcheck3, "2btag only (m#gamma #gamma fit w/ Higgs)", "LP");
  if(!low)leg->AddEntry(gcheck5, "m#gamma #gamma fit - 2btag only ", "LP");
//  if(low)leg->AddEntry(gcheck6, "m#gamma #gamma fit w/ Higgs", "LP");
//  leg->AddEntry(gcheck2, "#gamma #gamma Sidebands (stat only)", "LP");
  }
//  leg->SetHeader("categories combined");
 // leg->SetHeader("#splitline{#scale[1.0]{CMS Preliminary}}{categories combined | M jet | 19.6/pb}");
  leg->Draw();
 
  if(!low && !base){
  TLine l(400,0.01,400,350);
  l->Draw();
  }
//  mg->GetYaxis()->SetRangeUser(0.005,2);
//  mg->GetXaxis()->SetLimits(1,20);
  if(HH && low ) c1->SaveAs("WP4_cutbased_HH_low.png"); // HH
  if(HH && low ) c1->SaveAs("WP4_cutbased_HH_low.pdf"); // HH
  if(!HH && low && !base) c1->SaveAs("WP4_cutbased_low.png");
  if(!HH && low && !base) c1->SaveAs("WP4_cutbased_low.pdf");
  if(HH && !low ) c1->SaveAs("WP4_cutbased_HH.png"); // HH
  if(HH && !low ) c1->SaveAs("WP4_cutbased_HH.pdf"); // HH
  if(!HH && !low && !base) c1->SaveAs("WP4_cutbased.png");
  if(!HH && !low && !base) c1->SaveAs("WP4_cutbased.pdf");
  //
  if(!HH && low && base) c1->SaveAs("WP4_cutbased_low_base.png"); // HH
  if(!HH && low && base) c1->SaveAs("WP4_cutbased_low_base.pdf"); // HH
  if(!HH && !low && base && !oldExp) c1->SaveAs("WP4_cutbased_base.png");
  if(!HH && !low && base && !oldExp) c1->SaveAs("WP4_cutbased_base.pdf");
  if(!HH && !low && base && oldExp) c1->SaveAs("WP4_cutbased_base_oldExp.png");
  if(!HH && !low && base && oldExp) c1->SaveAs("WP4_cutbased_base_oldExp.pdf");
   //c1->SaveAs("radionlim/Mcut_cutbased.png");
   return c1;



}
