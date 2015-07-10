mtot[320,1200];

mtot_sig_m0_cat0[605, 580, 650];
mtot_sig_sigma_cat0[8., 5., 15.0];
mtot_sig_alpha_cat0[-0.1, -4.0, 4.0];
mtot_sig_n_cat0[4.0, 0.0, 20];
mtot_sig_gsigma_cat0[20., 10., 35.0];
mtot_sig_frac_cat0[0.5, 0.3, 1.0];

mtot_sig_m0_cat1[630, 550, 650];
mtot_sig_sigma_cat1[15.0, 0.5, 30.0];
mtot_sig_alpha_cat1[-1.0, -1.0, 3.0];
mtot_sig_n_cat1[10.0, 0.4, 20];
mtot_sig_gsigma_cat1[15.0, 1.0, 20.0];
mtot_sig_frac_cat1[0.8, 0.5, 0.9];

mtotGaussSig_cat0 = Gaussian(mtot, mtot_sig_m0_cat0, mtot_sig_gsigma_cat0);
mtotCBSig_cat0 = CBShape(mtot, mtot_sig_m0_cat0, mtot_sig_sigma_cat0, mtot_sig_alpha_cat0, mtot_sig_n_cat0);
mtotSig_cat0 = AddPdf(mtotGaussSig_cat0, mtotCBSig_cat0, mtot_sig_frac_cat0);

mtotGaussSig_cat1 = Gaussian(mtot, mtot_sig_m0_cat1, mtot_sig_gsigma_cat1);
mtotCBSig_cat1 = CBShape(mtot, mtot_sig_m0_cat1, mtot_sig_sigma_cat1, mtot_sig_alpha_cat1, mtot_sig_n_cat1);
mtotSig_cat1 = AddPdf(mtotGaussSig_cat1, mtotCBSig_cat1, mtot_sig_frac_cat1);

mtot_bkg_8TeV_slope1_cat0[1.9, -100., 100.];

mtot_bkg_8TeV_slope1_cat1[1.9, -100., 100.];
mtot_bkg_8TeV_slope2_cat1[1000,800, 1200];

wei[1,0,10];
sqrtS[8000., 8000., 8000., 8000.]
