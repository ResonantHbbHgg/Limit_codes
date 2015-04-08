mgg[100,180];


mgg_bkg_8TeV_slope1[9,-10.0, 20.0];
mgg_bkg_8TeV_slope1_cat0[9,-10.0, 20.0];
mgg_bkg_8TeV_slope1_cat1[9,-10.0, 20.0];
mgg_bkg_8TeV_slope1_cat2[9,-10.0, 20.0];
mgg_bkg_8TeV_slope1_cat3[9,-10.0, 20.0];

mgg_sig_m0_cat0[124.2, 123, 125];
mgg_sig_sigma_cat0[2.0, 1.0, 3.0];
mgg_sig_alpha_cat0[1.0, 0.0, 2.5];
mgg_sig_n_cat0[2.0, 1.0, 5.0];
mgg_sig_gsigma_cat0[1.2, 0.8, 1.6];
mgg_sig_frac_cat0[0.6, 0.4, 1.0];

mgg_sig_m0_cat1[124.2, 123, 125];
mgg_sig_sigma_cat1[2.0, 1.0, 3.0];
mgg_sig_alpha_cat1[1.0, 0.0, 2.5];
mgg_sig_n_cat1[2.0, 1.5, 10];
mgg_sig_gsigma_cat1[1.2, 0.8, 1.6];
mgg_sig_frac_cat1[0.6, 0.4, 1.0];

mgg_sig_m0_cat2[124.2, 123, 125];
mgg_sig_sigma_cat2[2.0, 1.0, 3.0];
mgg_sig_alpha_cat2[1.0, 0.0, 2.5];
mgg_sig_n_cat2[2.0, 1.0, 5.0];
mgg_sig_gsigma_cat2[1.2, 0.8, 1.6];
mgg_sig_frac_cat2[0.6, 0.4, 1.0];

mgg_sig_m0_cat3[124.2, 123, 125];
mgg_sig_sigma_cat3[2.0, 1.0, 3.0];
mgg_sig_alpha_cat3[1.0, 0.0, 2.5];
mgg_sig_n_cat3[2.0, 1.5, 10];
mgg_sig_gsigma_cat3[1.2, 0.8, 1.6];
mgg_sig_frac_cat3[0.6, 0.4, 1.0];

mggGaussSig_cat0 = Gaussian(mgg, mgg_sig_m0_cat0, mgg_sig_gsigma_cat0);
mggCBSig_cat0 = CBShape(mgg, mgg_sig_m0_cat0, mgg_sig_sigma_cat0, mgg_sig_alpha_cat0, mgg_sig_n_cat0);
mggSig_cat0 = AddPdf(mggGaussSig_cat0, mggCBSig_cat0, mgg_sig_frac_cat0);

mggGaussSig_cat1 = Gaussian(mgg, mgg_sig_m0_cat1, mgg_sig_gsigma_cat1);
mggCBSig_cat1 = CBShape(mgg, mgg_sig_m0_cat1, mgg_sig_sigma_cat1, mgg_sig_alpha_cat1, mgg_sig_n_cat1);
mggSig_cat1 = AddPdf(mggGaussSig_cat1, mggCBSig_cat1, mgg_sig_frac_cat1);

mggGaussSig_cat2 = Gaussian(mgg, mgg_sig_m0_cat2, mgg_sig_gsigma_cat2);
mggCBSig_cat2 = CBShape(mgg, mgg_sig_m0_cat2, mgg_sig_sigma_cat2, mgg_sig_alpha_cat2, mgg_sig_n_cat2);
mggSig_cat2 = AddPdf(mggGaussSig_cat2, mggCBSig_cat2, mgg_sig_frac_cat2);

mggGaussSig_cat3 = Gaussian(mgg, mgg_sig_m0_cat3, mgg_sig_gsigma_cat3);
mggCBSig_cat3 = CBShape(mgg, mgg_sig_m0_cat3, mgg_sig_sigma_cat3, mgg_sig_alpha_cat3, mgg_sig_n_cat3);
mggSig_cat3 = AddPdf(mggGaussSig_cat3, mggCBSig_cat3, mgg_sig_frac_cat3);

mgg_hig_m0_0_cat0[124.2, 123, 125];
mgg_hig_sigma_0_cat0[2.0, 1.0, 3.0];
mgg_hig_alpha_0_cat0[1.0, 0.0, 2.5];
mgg_hig_n_0_cat0[2.0, 1.0, 5.0];
mgg_hig_gsigma_0_cat0[1.2, 0.8, 1.6];
mgg_hig_frac_0_cat0[0.6, 0.4, 1.0];

mgg_hig_m0_0_cat1[124.2, 123, 125];
mgg_hig_sigma_0_cat1[2.0, 1.0, 3.0];
mgg_hig_alpha_0_cat1[1.0, 0.0, 2.5];
mgg_hig_n_0_cat1[2.0, 1.5, 10];
mgg_hig_gsigma_0_cat1[1.2, 0.8, 1.6];
mgg_hig_frac_0_cat1[0.6, 0.4, 1.0];

mgg_hig_m0_0_cat2[124.2, 123, 125];
mgg_hig_sigma_0_cat2[2.0, 1.0, 3.0];
mgg_hig_alpha_0_cat2[1.0, 0.0, 2.5];
mgg_hig_n_0_cat2[2.0, 1.0, 5.0];
mgg_hig_gsigma_0_cat2[1.2, 0.8, 1.6];
mgg_hig_frac_0_cat2[0.6, 0.4, 1.0];

mgg_hig_m0_0_cat3[124.2, 123, 125];
mgg_hig_sigma_0_cat3[2.0, 1.0, 3.0];
mgg_hig_alpha_0_cat3[1.0, 0.0, 2.5];
mgg_hig_n_0_cat3[2.0, 1.5, 10];
mgg_hig_gsigma_0_cat3[1.2, 0.8, 1.6];
mgg_hig_frac_0_cat3[0.6, 0.4, 1.0];

mggGaussHig_0_cat0 = Gaussian(mgg, mgg_hig_m0_0_cat0, mgg_hig_gsigma_0_cat0);
mggCBHig_0_cat0 = CBShape(mgg, mgg_hig_m0_0_cat0, mgg_hig_sigma_0_cat0, mgg_hig_alpha_0_cat0, mgg_hig_n_0_cat0);
mggHig_0_cat0 = AddPdf(mggGaussHig_0_cat0, mggCBHig_0_cat0, mgg_hig_frac_0_cat0);

mggGaussHig_0_cat1 = Gaussian(mgg, mgg_hig_m0_0_cat1, mgg_hig_gsigma_0_cat1);
mggCBHig_0_cat1 = CBShape(mgg, mgg_hig_m0_0_cat1, mgg_hig_sigma_0_cat1, mgg_hig_alpha_0_cat1, mgg_hig_n_0_cat1);
mggHig_0_cat1 = AddPdf(mggGaussHig_0_cat1, mggCBHig_0_cat1, mgg_hig_frac_0_cat1);

mggGaussHig_0_cat2 = Gaussian(mgg, mgg_hig_m0_0_cat2, mgg_hig_gsigma_0_cat2);
mggCBHig_0_cat2 = CBShape(mgg, mgg_hig_m0_0_cat2, mgg_hig_sigma_0_cat2, mgg_hig_alpha_0_cat2, mgg_hig_n_0_cat2);
mggHig_0_cat2 = AddPdf(mggGaussHig_0_cat2, mggCBHig_0_cat2, mgg_hig_frac_0_cat2);

mggGaussHig_0_cat3 = Gaussian(mgg, mgg_hig_m0_0_cat3, mgg_hig_gsigma_0_cat3);
mggCBHig_0_cat3 = CBShape(mgg, mgg_hig_m0_0_cat3, mgg_hig_sigma_0_cat3, mgg_hig_alpha_0_cat3, mgg_hig_n_0_cat3);
mggHig_0_cat3 = AddPdf(mggGaussHig_0_cat3, mggCBHig_0_cat3, mgg_hig_frac_0_cat3);

mgg_hig_m0_1_cat0[124.2, 123, 125];
mgg_hig_sigma_1_cat0[2.0, 1.0, 3.0];
mgg_hig_alpha_1_cat0[1.0, 0.0, 2.5];
mgg_hig_n_1_cat0[2.0, 1.0, 5.0];
mgg_hig_gsigma_1_cat0[1.2, 0.8, 1.6];
mgg_hig_frac_1_cat0[0.6, 0.4, 1.0];

mgg_hig_m0_1_cat1[124.2, 123, 125];
mgg_hig_sigma_1_cat1[2.0, 1.0, 3.0];
mgg_hig_alpha_1_cat1[1.0, 0.0, 2.5];
mgg_hig_n_1_cat1[2.0, 1.5, 10];
mgg_hig_gsigma_1_cat1[1.2, 0.8, 1.6];
mgg_hig_frac_1_cat1[0.6, 0.4, 1.0];

mgg_hig_m0_1_cat2[124.2, 123, 125];
mgg_hig_sigma_1_cat2[2.0, 1.0, 3.0];
mgg_hig_alpha_1_cat2[1.0, 0.0, 2.5];
mgg_hig_n_1_cat2[2.0, 1.0, 5.0];
mgg_hig_gsigma_1_cat2[1.2, 0.8, 1.6];
mgg_hig_frac_1_cat2[0.6, 0.4, 1.0];

mgg_hig_m0_1_cat3[124.2, 123, 125];
mgg_hig_sigma_1_cat3[2.0, 1.0, 3.0];
mgg_hig_alpha_1_cat3[1.0, 0.0, 2.5];
mgg_hig_n_1_cat3[2.0, 1.5, 10];
mgg_hig_gsigma_1_cat3[1.2, 0.8, 1.6];
mgg_hig_frac_1_cat3[0.6, 0.4, 1.0];

mggGaussHig_1_cat0 = Gaussian(mgg, mgg_hig_m0_1_cat0, mgg_hig_gsigma_1_cat0);
mggCBHig_1_cat0 = CBShape(mgg, mgg_hig_m0_1_cat0, mgg_hig_sigma_1_cat0, mgg_hig_alpha_1_cat0, mgg_hig_n_1_cat0);
mggHig_1_cat0 = AddPdf(mggGaussHig_1_cat0, mggCBHig_1_cat0, mgg_hig_frac_1_cat0);

mggGaussHig_1_cat1 = Gaussian(mgg, mgg_hig_m0_1_cat1, mgg_hig_gsigma_1_cat1);
mggCBHig_1_cat1 = CBShape(mgg, mgg_hig_m0_1_cat1, mgg_hig_sigma_1_cat1, mgg_hig_alpha_1_cat1, mgg_hig_n_1_cat1);
mggHig_1_cat1 = AddPdf(mggGaussHig_1_cat1, mggCBHig_1_cat1, mgg_hig_frac_1_cat1);

mggGaussHig_1_cat2 = Gaussian(mgg, mgg_hig_m0_1_cat2, mgg_hig_gsigma_1_cat2);
mggCBHig_1_cat2 = CBShape(mgg, mgg_hig_m0_1_cat2, mgg_hig_sigma_1_cat2, mgg_hig_alpha_1_cat2, mgg_hig_n_1_cat2);
mggHig_1_cat2 = AddPdf(mggGaussHig_1_cat2, mggCBHig_1_cat2, mgg_hig_frac_1_cat2);

mggGaussHig_1_cat3 = Gaussian(mgg, mgg_hig_m0_1_cat3, mgg_hig_gsigma_1_cat3);
mggCBHig_1_cat3 = CBShape(mgg, mgg_hig_m0_1_cat3, mgg_hig_sigma_1_cat3, mgg_hig_alpha_1_cat3, mgg_hig_n_1_cat3);
mggHig_1_cat3 = AddPdf(mggGaussHig_1_cat3, mggCBHig_1_cat3, mgg_hig_frac_1_cat3);

mgg_hig_m0_2_cat0[124.2, 123, 125];
mgg_hig_sigma_2_cat0[2.0, 1.0, 3.0];
mgg_hig_alpha_2_cat0[1.0, 0.0, 2.5];
mgg_hig_n_2_cat0[2.0, 1.0, 5.0];
mgg_hig_gsigma_2_cat0[1.2, 0.8, 1.6];
mgg_hig_frac_2_cat0[0.6, 0.4, 1.0];

mgg_hig_m0_2_cat1[124.2, 123, 125];
mgg_hig_sigma_2_cat1[2.0, 1.0, 3.0];
mgg_hig_alpha_2_cat1[1.0, 0.0, 2.5];
mgg_hig_n_2_cat1[2.0, 1.5, 10];
mgg_hig_gsigma_2_cat1[1.2, 0.8, 1.6];
mgg_hig_frac_2_cat1[0.6, 0.4, 1.0];

mgg_hig_m0_2_cat2[124.2, 123, 125];
mgg_hig_sigma_2_cat2[2.0, 1.0, 3.0];
mgg_hig_alpha_2_cat2[1.0, 0.0, 2.5];
mgg_hig_n_2_cat2[2.0, 1.0, 5.0];
mgg_hig_gsigma_2_cat2[1.2, 0.8, 1.6];
mgg_hig_frac_2_cat2[0.6, 0.4, 1.0];

mgg_hig_m0_2_cat3[124.2, 123, 125];
mgg_hig_sigma_2_cat3[2.0, 1.0, 3.0];
mgg_hig_alpha_2_cat3[1.0, 0.0, 2.5];
mgg_hig_n_2_cat3[2.0, 1.5, 10];
mgg_hig_gsigma_2_cat3[1.2, 0.8, 1.6];
mgg_hig_frac_2_cat3[0.6, 0.4, 1.0];

mggGaussHig_2_cat0 = Gaussian(mgg, mgg_hig_m0_2_cat0, mgg_hig_gsigma_2_cat0);
mggCBHig_2_cat0 = CBShape(mgg, mgg_hig_m0_2_cat0, mgg_hig_sigma_2_cat0, mgg_hig_alpha_2_cat0, mgg_hig_n_2_cat0);
mggHig_2_cat0 = AddPdf(mggGaussHig_2_cat0, mggCBHig_2_cat0, mgg_hig_frac_2_cat0);

mggGaussHig_2_cat1 = Gaussian(mgg, mgg_hig_m0_2_cat1, mgg_hig_gsigma_2_cat1);
mggCBHig_2_cat1 = CBShape(mgg, mgg_hig_m0_2_cat1, mgg_hig_sigma_2_cat1, mgg_hig_alpha_2_cat1, mgg_hig_n_2_cat1);
mggHig_2_cat1 = AddPdf(mggGaussHig_2_cat1, mggCBHig_2_cat1, mgg_hig_frac_2_cat1);

mggGaussHig_2_cat2 = Gaussian(mgg, mgg_hig_m0_2_cat2, mgg_hig_gsigma_2_cat2);
mggCBHig_2_cat2 = CBShape(mgg, mgg_hig_m0_2_cat2, mgg_hig_sigma_2_cat2, mgg_hig_alpha_2_cat2, mgg_hig_n_2_cat2);
mggHig_2_cat2 = AddPdf(mggGaussHig_2_cat2, mggCBHig_2_cat2, mgg_hig_frac_2_cat2);

mggGaussHig_2_cat3 = Gaussian(mgg, mgg_hig_m0_2_cat3, mgg_hig_gsigma_2_cat3);
mggCBHig_2_cat3 = CBShape(mgg, mgg_hig_m0_2_cat3, mgg_hig_sigma_2_cat3, mgg_hig_alpha_2_cat3, mgg_hig_n_2_cat3);
mggHig_2_cat3 = AddPdf(mggGaussHig_2_cat3, mggCBHig_2_cat3, mgg_hig_frac_2_cat3);

mgg_hig_m0_3_cat0[124.2, 123, 125];
mgg_hig_sigma_3_cat0[2.0, 1.0, 3.0];
mgg_hig_alpha_3_cat0[1.0, 0.0, 2.5];
mgg_hig_n_3_cat0[2.0, 1.0, 5.0];
mgg_hig_gsigma_3_cat0[1.2, 0.8, 1.6];
mgg_hig_frac_3_cat0[0.6, 0.4, 1.0];

mgg_hig_m0_3_cat1[124.2, 123, 125];
mgg_hig_sigma_3_cat1[2.0, 1.0, 3.0];
mgg_hig_alpha_3_cat1[1.0, 0.0, 2.5];
mgg_hig_n_3_cat1[2.0, 1.5, 10];
mgg_hig_gsigma_3_cat1[1.2, 0.8, 1.6];
mgg_hig_frac_3_cat1[0.6, 0.4, 1.0];

mgg_hig_m0_3_cat2[124.2, 123, 125];
mgg_hig_sigma_3_cat2[2.0, 1.0, 3.0];
mgg_hig_alpha_3_cat2[1.0, 0.0, 2.5];
mgg_hig_n_3_cat2[2.0, 1.0, 5.0];
mgg_hig_gsigma_3_cat2[1.2, 0.8, 1.6];
mgg_hig_frac_3_cat2[0.6, 0.4, 1.0];

mgg_hig_m0_3_cat3[124.2, 123, 125];
mgg_hig_sigma_3_cat3[2.0, 1.0, 3.0];
mgg_hig_alpha_3_cat3[1.0, 0.0, 2.5];
mgg_hig_n_3_cat3[2.0, 1.5, 10];
mgg_hig_gsigma_3_cat3[1.2, 0.8, 1.6];
mgg_hig_frac_3_cat3[0.6, 0.4, 1.0];

mggGaussHig_3_cat0 = Gaussian(mgg, mgg_hig_m0_3_cat0, mgg_hig_gsigma_3_cat0);
mggCBHig_3_cat0 = CBShape(mgg, mgg_hig_m0_3_cat0, mgg_hig_sigma_3_cat0, mgg_hig_alpha_3_cat0, mgg_hig_n_3_cat0);
mggHig_3_cat0 = AddPdf(mggGaussHig_3_cat0, mggCBHig_3_cat0, mgg_hig_frac_3_cat0);

mggGaussHig_3_cat1 = Gaussian(mgg, mgg_hig_m0_3_cat1, mgg_hig_gsigma_3_cat1);
mggCBHig_3_cat1 = CBShape(mgg, mgg_hig_m0_3_cat1, mgg_hig_sigma_3_cat1, mgg_hig_alpha_3_cat1, mgg_hig_n_3_cat1);
mggHig_3_cat1 = AddPdf(mggGaussHig_3_cat1, mggCBHig_3_cat1, mgg_hig_frac_3_cat1);

mggGaussHig_3_cat2 = Gaussian(mgg, mgg_hig_m0_3_cat2, mgg_hig_gsigma_3_cat2);
mggCBHig_3_cat2 = CBShape(mgg, mgg_hig_m0_3_cat2, mgg_hig_sigma_3_cat2, mgg_hig_alpha_3_cat2, mgg_hig_n_3_cat2);
mggHig_3_cat2 = AddPdf(mggGaussHig_3_cat2, mggCBHig_3_cat2, mgg_hig_frac_3_cat2);

mggGaussHig_3_cat3 = Gaussian(mgg, mgg_hig_m0_3_cat3, mgg_hig_gsigma_3_cat3);
mggCBHig_3_cat3 = CBShape(mgg, mgg_hig_m0_3_cat3, mgg_hig_sigma_3_cat3, mgg_hig_alpha_3_cat3, mgg_hig_n_3_cat3);
mggHig_3_cat3 = AddPdf(mggGaussHig_3_cat3, mggCBHig_3_cat3, mgg_hig_frac_3_cat3);

mgg_hig_m0_4_cat0[124.2, 123, 125];
mgg_hig_sigma_4_cat0[2.0, 1.0, 3.0];
mgg_hig_alpha_4_cat0[1.0, 0.0, 2.5];
mgg_hig_n_4_cat0[2.0, 1.0, 5.0];
mgg_hig_gsigma_4_cat0[1.2, 0.8, 1.6];
mgg_hig_frac_4_cat0[0.6, 0.4, 1.0];

mgg_hig_m0_4_cat1[124.2, 123, 125];
mgg_hig_sigma_4_cat1[2.0, 1.0, 3.0];
mgg_hig_alpha_4_cat1[1.0, 0.0, 2.5];
mgg_hig_n_4_cat1[2.0, 1.5, 10];
mgg_hig_gsigma_4_cat1[1.2, 0.8, 1.6];
mgg_hig_frac_4_cat1[0.6, 0.4, 1.0];

mgg_hig_m0_4_cat2[124.2, 123, 125];
mgg_hig_sigma_4_cat2[2.0, 1.0, 3.0];
mgg_hig_alpha_4_cat2[1.0, 0.0, 2.5];
mgg_hig_n_4_cat2[2.0, 1.0, 5.0];
mgg_hig_gsigma_4_cat2[1.2, 0.8, 1.6];
mgg_hig_frac_4_cat2[0.6, 0.4, 1.0];

mgg_hig_m0_4_cat3[124.2, 123, 125];
mgg_hig_sigma_4_cat3[2.0, 1.0, 3.0];
mgg_hig_alpha_4_cat3[1.0, 0.0, 2.5];
mgg_hig_n_4_cat3[2.0, 1.5, 10];
mgg_hig_gsigma_4_cat3[1.2, 0.8, 1.6];
mgg_hig_frac_4_cat3[0.6, 0.4, 1.0];

mggGaussHig_4_cat0 = Gaussian(mgg, mgg_hig_m0_4_cat0, mgg_hig_gsigma_4_cat0);
mggCBHig_4_cat0 = CBShape(mgg, mgg_hig_m0_4_cat0, mgg_hig_sigma_4_cat0, mgg_hig_alpha_4_cat0, mgg_hig_n_4_cat0);
mggHig_4_cat0 = AddPdf(mggGaussHig_4_cat0, mggCBHig_4_cat0, mgg_hig_frac_4_cat0);

mggGaussHig_4_cat1 = Gaussian(mgg, mgg_hig_m0_4_cat1, mgg_hig_gsigma_4_cat1);
mggCBHig_4_cat1 = CBShape(mgg, mgg_hig_m0_4_cat1, mgg_hig_sigma_4_cat1, mgg_hig_alpha_4_cat1, mgg_hig_n_4_cat1);
mggHig_4_cat1 = AddPdf(mggGaussHig_4_cat1, mggCBHig_4_cat1, mgg_hig_frac_4_cat1);

mggGaussHig_4_cat2 = Gaussian(mgg, mgg_hig_m0_4_cat2, mgg_hig_gsigma_4_cat2);
mggCBHig_4_cat2 = CBShape(mgg, mgg_hig_m0_4_cat2, mgg_hig_sigma_4_cat2, mgg_hig_alpha_4_cat2, mgg_hig_n_4_cat2);
mggHig_4_cat2 = AddPdf(mggGaussHig_4_cat2, mggCBHig_4_cat2, mgg_hig_frac_4_cat2);

mggGaussHig_4_cat3 = Gaussian(mgg, mgg_hig_m0_4_cat3, mgg_hig_gsigma_4_cat3);
mggCBHig_4_cat3 = CBShape(mgg, mgg_hig_m0_4_cat3, mgg_hig_sigma_4_cat3, mgg_hig_alpha_4_cat3, mgg_hig_n_4_cat3);
mggHig_4_cat3 = AddPdf(mggGaussHig_4_cat3, mggCBHig_4_cat3, mgg_hig_frac_4_cat3);
