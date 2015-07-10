mtot[320,1200];

mtot_sig_m0_cat0[810, 800, 830];
mtot_sig_sigma_cat0[5., 3., 15.0];
mtot_sig_alpha_cat0[-1.5, -4.0, -0.5];
mtot_sig_n_cat0[4.0, 0.0, 20];
mtot_sig_gsigma_cat0[30., 1., 40.0];
mtot_sig_frac_cat0[0.4, 0.0, 0.5];

mtot_sig_m0_cat1[820, 810, 830];
mtot_sig_sigma_cat1[5.0, 0.5, 15.0];
mtot_sig_alpha_cat1[-1.5, -3.0, -0.5];
mtot_sig_n_cat1[10.0, 0.4, 20];
mtot_sig_gsigma_cat1[25.0, 1.0, 40.0];
mtot_sig_frac_cat1[0.4, 0.0, 0.6];

mtotGaussSig_cat0 = Gaussian(mtot, mtot_sig_m0_cat0, mtot_sig_gsigma_cat0);
mtotCBSig_cat0 = CBShape(mtot, mtot_sig_m0_cat0, mtot_sig_sigma_cat0, mtot_sig_alpha_cat0, mtot_sig_n_cat0);
mtotSig_cat0 = AddPdf(mtotGaussSig_cat0, mtotCBSig_cat0, mtot_sig_frac_cat0);

mtotGaussSig_cat1 = Gaussian(mtot, mtot_sig_m0_cat1, mtot_sig_gsigma_cat1);
mtotCBSig_cat1 = CBShape(mtot, mtot_sig_m0_cat1, mtot_sig_sigma_cat1, mtot_sig_alpha_cat1, mtot_sig_n_cat1);
mtotSig_cat1 = AddPdf(mtotGaussSig_cat1, mtotCBSig_cat1, mtot_sig_frac_cat1);

mtot_bkg_8TeV_slope1_cat0[1.9, -10., 10.];

mtot_bkg_8TeV_slope1_cat1[1.9, -10., 10.];
mtot_bkg_8TeV_slope2_cat1[1000,800, 1200];

wei[1,0,10];
sqrtS[8000., 8000., 8000., 8000.]
