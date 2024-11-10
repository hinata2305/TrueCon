%Copyright © 2024 Koten and Schüppen All rights reserved
%Important Notice: This code is not intended for medical applications 
%and does not have legal approval for such use. We strongly recommend 
%using FDA-approved software for any medical purposes. 

% This function estimates the confidence interval called conf of a
% reliability measure (rel) before and after the measure was corrected for serial correlations

function [dat] = remove_AutoBlockBoot_short(fMRT1, fMRT2, numboot)

s = size(fMRT1);

% initialize array to store raw reliability measure
relROIRawTrans = zeros(s(2),1);
% initialize array to store lower bound of confidence interval for raw reliability measure
relROIRawTransLow = zeros(s(2),1);
% initialize array to store upper bound of confidence interval for raw reliability measure
relROIRawTransUp = zeros(s(2),1);
% initialize array to store reliability measure after removing serial autocorrelation
relROIRemAutoTrans = zeros(s(2),1);
% initialize array to store lower bound of confidence interval for reliability measure after removing serial autocorrelation
relROIRemAutoTransLow = zeros(s(2),1);
% initialize array to store upper bound of confidence interval for reliability measure after removing serial autocorrelation
relROIRemAutoTransUp = zeros(s(2),1);


% region of interest loop

parfor k = 1:s(2)


    ROItargetT = fMRT1(:,k); % get test data for region of interest
    ROItargetR = fMRT2(:,k); % get retest data for region of interest

    % estimate confidence interval for raw reliability measure
    conf = bootstrapBlock(numboot, ROItargetT, ROItargetR);
    relROIRawTrans(k) = atanh(corr(ROItargetT, ROItargetR)); % calculate raw reliability measure
    relROIRawTransLow(k) = conf.Rel_Trans_LB; % store lower bound of confidence interval
    relROIRawTransUp(k) = conf.Rel_Trans_UB; % store upper bound of confidence interval

    % remove serial autocorrelation
    [remAutoRelTargetT, remAutoRelTargetR, ~] = autocorfilter(ROItargetT, ROItargetR);

    % estimate confidence interval for reliability measure after removing serial autocorrelation
    conf = bootstrapBlock(numboot, remAutoRelTargetT, remAutoRelTargetR);
    relROIRemAutoTrans(k) = atanh(corr(remAutoRelTargetT, remAutoRelTargetR)); % calculate reliability measure after removing serial autocorrelation
    relROIRemAutoTransLow(k) = conf.Rel_Trans_LB; % store lower bound of confidence interval
    relROIRemAutoTransUp(k) = conf.Rel_Trans_UB; % store upper bound of confidence interval

end

dat.relRoi_raw_trans=relROIRawTrans;
dat.relRoi_raw_Up_trans=relROIRawTransUp;
dat.relRoi_raw_Low_trans=relROIRawTransLow;

dat.relRoi_remAuto_trans=relROIRemAutoTrans;
dat.relRoi_remAuto_Up_trans=relROIRemAutoTransUp;
dat.relRoi_remAuto_Low_trans=relROIRemAutoTransLow;


end