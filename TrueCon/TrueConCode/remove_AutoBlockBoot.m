
%Copyright © 2024 Koten and Schüppen All rights reserved
%Important Notice: This code is not intended for medical applications 
%and does not have legal approval for such use. We strongly recommend 
%using FDA-approved software for any medical purposes. 


function [dat] = remove_AutoBlockBoot(fMRT1,fMRT2,numboot)

s = size(fMRT1);

% Preallocate result matrices
conTestRawTrans = nan(s(2), s(2));
conTestRawTransLow = nan(s(2), s(2));
conTestRawTransUp = nan(s(2), s(2));

conRetestRawTrans = nan(s(2), s(2));
conRetestRawTransLow = nan(s(2), s(2));
conRetestRawTransUp = nan(s(2), s(2));

reldifTestRawTrans = nan(s(2), s(2));
reldifRetestRawTrans = nan(s(2), s(2));
relmaxRawTrans = nan(s(2), s(2));
relmaxRawTransLow = nan(s(2), s(2));
relmaxRawTransUp = nan(s(2), s(2));

contrueTestRawTrans = nan(s(2), s(2));
contrueTestRawTransLow = nan(s(2), s(2));
contrueTestRawTransUp = nan(s(2), s(2));

contrueRetestRawTrans = nan(s(2), s(2));
contrueRetestRawTransLow = nan(s(2), s(2));
contrueRetestRawTransUp = nan(s(2), s(2));

conTestRemAutoTrans = nan(s(2), s(2));
conTestRemAutoTransLow = nan(s(2), s(2));
conTestRemAutoTransUp = nan(s(2), s(2));

conRetestRemAutoTrans = nan(s(2), s(2));
conRetestRemAutoTransLow = nan(s(2), s(2));
conRetestRemAutoTransUp = nan(s(2), s(2));

reldifTestRemAutoTrans = nan(s(2), s(2));
reldifRetestRemAutoTrans = nan(s(2), s(2));
relmaxRemAutoTrans = nan(s(2), s(2));
relmaxRemAutoTransLow = nan(s(2), s(2));
relmaxRemAutoTransUp = nan(s(2), s(2));

contrueTestRemAutoTrans = nan(s(2), s(2));
contrueTestRemAutoTransLow = nan(s(2), s(2));
contrueTestRemAutoTransUp = nan(s(2), s(2));

contrueRetestRemAutoTrans = nan(s(2), s(2));
contrueRetestRemAutoTransLow = nan(s(2), s(2));
contrueRetestRemAutoTransUp = nan(s(2), s(2));

% Run through the k*f size connecteome
for k = 1:s(2)

  ROItargetT = fMRT1(:, k);
  ROItargetR = fMRT2(:, k);

    parfor f = (k + 1):s(2)

    % Get timecourse of seed ROI for current f index
    ROIseedT = fMRT1(:, f);
    ROIseedR = fMRT2(:, f);

    % Estimate test-retest reliability for seed ROI
    relROIseedRaw_trans = atanh(corr(ROIseedT, ROIseedR));
    % Estimate test-retest reliability for target ROI
    relROItargetRaw_trans = atanh(corr(ROItargetT, ROItargetR));

    % Estimate the confidence intervals of con, true con, relmax
    conf = bootstrapBlock(numboot, ROIseedT, ROIseedR, ROItargetT, ROItargetR);

    % Estimate connectivity of raw data and apply first step of Fisher z transform test for test data
    conTestRaw_trans = atanh(corr(ROIseedT, ROItargetT));
    conTestRawTrans(f, k) = conTestRaw_trans;
    conTestRawTransLow(f, k) = conf.ConTestBoot_trans_LB;
    conTestRawTransUp(f, k) = conf.ConTestBoot_trans_UB;

    % Estimate connectivity of raw data and apply first step of Fisher z transform test for retest data
    conRetestRaw_trans = atanh(corr(ROIseedR, ROItargetR));
    conRetestRawTrans(f, k) = conRetestRaw_trans;
    conRetestRawTransLow(f, k) = conf.ConRetestBoot_trans_LB;
    conRetestRawTransUp(f, k) = conf.ConRetestBoot_trans_UB;

    % Estimate true connectivity of raw data
    [reldifTest, contrueTest, reldifRetest, contrueRetest, relmax] = TrueConSimple(conTestRaw_trans, conRetestRaw_trans, relROIseedRaw_trans, relROItargetRaw_trans);

    % Collect data
    reldifTestRawTrans(f, k) = reldifTest;
    contrueTestRawTrans(f, k) = contrueTest;
    contrueTestRawTransLow(f, k) = conf.contrueTestBoot_Trans_LB;
    contrueTestRawTransUp(f, k) = conf.contrueTestBoot_Trans_UB;
    reldifRetestRawTrans(f, k) = reldifRetest;
    contrueRetestRawTrans(f, k) = contrueRetest;
    contrueRetestRawTransLow(f, k) = conf.contrueRetestBoot_Trans_LB;
    contrueRetestRawTransUp(f, k) = conf.contrueRetestBoot_Trans_UB;

    relmaxRawTrans(f, k) = relmax;
    relmaxRawTransLow(f, k) = conf.relmaxBoot_Trans_LB;
    relmaxRawTransUp(f, k) = conf.relmaxBoot_Trans_UB;

    % Remove residual autocorrelation from data
    [remAutoSeedT, remAutoTargetT, ZetaT, ~] = autocorfilter(ROIseedT, ROItargetT);
    [remAutoSeedR, remAutoTargetR, ZetaR, ~] = autocorfilter(ROIseedR, ROItargetR);
    [remAutoRelSeedT, remAutoRelSeedR, ~] = autocorfilter(ROIseedT, ROIseedR);
    [remAutoRelTargetT, remAutoRelTargetR, ~] = autocorfilter(ROItargetT, ROItargetR);

    relROIseed_trans = atanh(corr(remAutoRelSeedT, remAutoRelSeedR));
    relROItarget_trans = atanh(corr(remAutoRelTargetT, remAutoRelTargetR));

    % Estimate the confidence intervals of con, true con, relmax
    conf = bootstrapBlock(numboot, remAutoSeedT, remAutoSeedR, remAutoTargetT, remAutoTargetR, remAutoRelSeedT, remAutoRelSeedR, remAutoRelTargetT, remAutoRelTargetR);

    % Estimate connectivity of data with removed residual autocorrelations and apply first step of Fisher z transform test
    conTest_trans = atanh(corr(remAutoSeedT, remAutoTargetT));
    conTestRemAutoTrans(f, k) = conTest_trans;
    conTestRemAutoTransLow(f, k) = conf.ConTestBoot_trans_LB;
    conTestRemAutoTransUp(f, k) = conf.ConTestBoot_trans_UB;

    conRetest_trans = atanh(corr(remAutoSeedR, remAutoTargetR));
    conRetestRemAutoTrans(f, k) = conRetest_trans;
    conRetestRemAutoTransLow(f, k) = conf.ConRetestBoot_trans_LB;
    conRetestRemAutoTransUp(f, k) = conf.ConRetestBoot_trans_UB;

    % Estimate true connectivity of data that are controlled for residual autocorrelations
    [reldifTest, contrueTest, reldifRetest, contrueRetest, relmax] = TrueConSimple(conTest_trans, conRetest_trans, relROIseed_trans, relROItarget_trans);

    reldifTestRemAutoTrans(f, k) = reldifTest;
    contrueTestRemAutoTrans(f, k) = contrueTest;
    contrueTestRemAutoTransUp(f, k) = conf.contrueTestBoot_Trans_UB;
    contrueTestRemAutoTransLow(f, k) = conf.contrueTestBoot_Trans_LB;
    reldifRetestRemAutoTrans(f, k) = reldifRetest;
    contrueRetestRemAutoTrans(f, k) = contrueRetest;
    contrueRetestRemAutoTransUp(f, k) = conf.contrueRetestBoot_Trans_UB;
    contrueRetestRemAutoTransLow(f, k) = conf.contrueRetestBoot_Trans_LB;

    relmaxRemAutoTrans(f, k) = relmax;
    relmaxRemAutoTransUp(f, k) = conf.relmaxBoot_Trans_UB;
    relmaxRemAutoTransLow(f, k) = conf.relmaxBoot_Trans_LB;
    end
end

% These are stats for the raw data
dat.connection_test_raw_trans=conTestRawTrans;
dat.connection_test_raw_Up_trans=conTestRawTransUp;
dat.connection_test_raw_Low_trans=conTestRawTransLow;
dat.connection_retest_raw_trans=conRetestRawTrans;
dat.connection_retest_raw_Up_trans=conRetestRawTransUp;
dat.connection_retest_raw_Low_trans=conRetestRawTransLow;

dat.reldif_test_raw_trans=reldifTestRawTrans;
dat.reldif_retest_raw_trans=reldifRetestRawTrans;

dat.relmax_raw_trans=relmaxRawTrans;
dat.relmax_raw_Up_trans=relmaxRawTransUp;
dat.relmax_raw_Low_trans=relmaxRawTransLow;

dat.contrue_test_raw_trans=contrueTestRawTrans ;
dat.contrue_test_raw_Up_trans=contrueTestRawTransUp ;
dat.contrue_test_raw_Low_trans=contrueTestRawTransLow ;
dat.contrue_retest_raw_trans=contrueRetestRawTrans ;
dat.contrue_retest_raw_Up_trans=contrueRetestRawTransUp;
dat.contrue_retest_raw_Low_trans=contrueRetestRawTransLow;

% These are stats for the AR(1) removed data
dat.connection_test_remAutoRoi_trans=conTestRemAutoTrans;
dat.connection_test_remAutoRoi_Up_trans=conTestRemAutoTransUp;
dat.connection_test_remAutoRoi_Low_trans=conTestRemAutoTransLow;
dat.connection_retest_remAutoRoi_trans=conRetestRemAutoTrans;
dat.connection_retest_remAutoRoi_Up_trans=conRetestRemAutoTransUp;
dat.connection_retest_remAutoRoi_Low_trans=conRetestRemAutoTransLow;

dat.reldif_test_remAutoRoi_trans=reldifTestRemAutoTrans;
dat.reldif_retest_remAutoRoi_trans=reldifRetestRemAutoTrans;

dat.relmax_remAutoRoi_trans=relmaxRemAutoTrans;
dat.relmax_remAutoRoi_Up_trans=relmaxRemAutoTransUp;
dat.relmax_remAutoRoi_Low_trans=relmaxRemAutoTransLow;

dat.contrue_test_remAutoRoi_trans=contrueTestRemAutoTrans;
dat.contrue_test_remAutoRoi_Up_trans= contrueTestRemAutoTransUp;
dat.contrue_test_remAutoRoi_Low_trans=contrueTestRemAutoTransLow;
dat.contrue_retest_remAutoRoi_trans=contrueRetestRemAutoTrans;
dat.contrue_retest_remAutoRoi_Up_trans=contrueRetestRemAutoTransUp;
dat.contrue_retest_remAutoRoi_Low_trans=contrueRetestRemAutoTransLow;



