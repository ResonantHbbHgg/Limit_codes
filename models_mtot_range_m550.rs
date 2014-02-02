mtot[300,1200];

mtot_sig_m0_cat0[550, 500, 600];
mtot_sig_sigma_cat0[6., 3., 10.0];
mtot_sig_alpha_cat0[-0.5, -4.0, 4.0];
mtot_sig_n_cat0[4.0, 0.1, 20];
mtot_sig_gsigma_cat0[15., 1., 25.0];
mtot_sig_frac_cat0[0.1, 0.0, 0.4];

mtot_sig_m0_cat1[550, 500, 600];
mtot_sig_sigma_cat1[20.0, 0.5, 30.0];
mtot_sig_alpha_cat1[-0.5, -1.0, 3.0];
mtot_sig_n_cat1[10.0, 0.4, 20];
mtot_sig_gsigma_cat1[8.0, 1.0, 20.0];
mtot_sig_frac_cat1[0.5, 0.4, 0.8];

mtotGaussSig_cat0 = Gaussian(mtot, mtot_sig_m0_cat0, mtot_sig_gsigma_cat0);
mtotCBSig_cat0 = CBShape(mtot, mtot_sig_m0_cat0, mtot_sig_sigma_cat0, mtot_sig_alpha_cat0, mtot_sig_n_cat0);
mtotSig_cat0 = AddPdf(mtotGaussSig_cat0, mtotCBSig_cat0, mtot_sig_frac_cat0);

mtotGaussSig_cat1 = Gaussian(mtot, mtot_sig_m0_cat1, mtot_sig_gsigma_cat1);
mtotCBSig_cat1 = CBShape(mtot, mtot_sig_m0_cat1, mtot_sig_sigma_cat1, mtot_sig_alpha_cat1, mtot_sig_n_cat1);
mtotSig_cat1 = AddPdf(mtotGaussSig_cat1, mtotCBSig_cat1, mtot_sig_frac_cat1);

mtot_bkg_8TeV_slope1_cat0[1.9, -5., 5.];
mtot_bkg_8TeV_slope2_cat0[0.,0., 0.];
mtot_bkg_8TeV_slope3_cat0[0.,-10.0, 10.0];

mtot_bkg_8TeV_slope1_cat1[1.9, -5., 5.];
mtot_bkg_8TeV_slope2_cat1[1000000,800000, 1200000];
mtot_bkg_8TeV_slope3_cat1[0.,-10.0, 10.0];

mtot_bkg_8TeV_slope_cat2[10.0,-500.0, 500.0];
mtot_bkg_8TeV_slope1_cat2[1., -500.0, 500.0];
mtot_bkg_8TeV_slope2_cat2[.0,-100.0, 100.0];
mtot_bkg_8TeV_slope3_cat2[0.,-10.0, 10.0];

wei[1,0,10];
sqrtS[8000., 8000., 8000., 8000.]
