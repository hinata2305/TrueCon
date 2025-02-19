
%Copyright © 2024 Koten and Schüppen All rights reserved
%Important Notice: This code is not intended for medical applications 
%and does not have legal approval for such use. We strongly recommend 
%using FDA-approved software for any medical purposes. 

% This function estimates reliability, connectivity, and detectable connectivity. It performs AR(1) corrections
% and conducts block bootstrapping at the time course level and
% conventional bootstrapping at the population level.

% Find the path where the package is located.  
folder = pwd;
cd(folder);
path = fullfile(folder, 'Stat_Result');

% Load fMRI test and retest data
load fMRI_clean.mat
% Load the subject-wise test-retest reliability of the response time experiment
load relResponseTime.mat
% Load the code names of the brain regions
load fwf_region_name.mat
% Load the code names of the subjects that are relevant
load namelist_n50smaller3mm_tasksmaller3wrong.mat

% Size of the connectivity matrix
numroi = size(Clean_SPM_p_1, 2);
s = [numroi numroi];
% Matrix used to extract the lower triangle of connectivity data
extractor = tril(true(s(2)), -1);

% Threshold span for Dice overlap
l = 250;

% Number of block bootstraps
numboot = 250;
% Number of bootstraps at sample level
NumSimulation = 250;
% Size of sample
nSamp = 50;

% The subject loop
parfor v = 1:nSamp
    
    % Name of subject
    subject = namelist{v};
    % Read in time courses from a test and retest run per subject
    timesec1 = Clean_SPM_p_1(:, :, v);
    timesec2 = Clean_SPM_p_2(:, :, v);
    
    % Estimate test-retest reliability with confidence intervals from raw time courses
    rmi = remove_AutoBlockBoot_short(timesec1, timesec2, numboot);
    
    RelTrans(:, v) = rmi.relRoi_raw_trans;
    RelTransUp(:, v) = rmi.relRoi_raw_Up_trans;
    RelTransLow(:, v) = rmi.relRoi_raw_Low_trans;
    % Estimate test-retest reliability with confidence intervals for time courses where auto
    % correlations are removed
    RelRemAutoTrans(:, v) = rmi.relRoi_remAuto_trans;
    RelRemAutoTransUp(:, v) = rmi.relRoi_remAuto_Up_trans;
    RelRemAutoTransLow(:, v) = rmi.relRoi_remAuto_Low_trans;
    
    % Estimate (true = detectable) connectivity with confidence intervals from raw time courses
    % and time courses that were corrected for serial correlations.
    mri = remove_AutoBlockBoot(timesec1, timesec2, numboot);
    
    % Connectivity test  
    % For convenience, some matrices are transformed into 1D and gathered into 2D
    connection_test_raw_temp = mri.connection_test_raw_trans;
    connectionTestRawTrans(v, :, :) = connection_test_raw_temp;
    connectionTestRaw_Trans2D(v, :) = connection_test_raw_temp(extractor);
    % Upper and lower bound confidence intervals of connectivity
    connectionTestRaw_Up_Trans2D(v, :) = mri.connection_test_raw_Up_trans(extractor);
    connectionTestRaw_Low_Trans2D(v, :) = mri.connection_test_raw_Low_trans(extractor);
    %retest
    connection_retest_raw_temp = mri.connection_retest_raw_trans;
    connectionRetestRawTrans(v, :, :) = connection_retest_raw_temp;
    connectionRetestRaw_Trans2D(v, :) = connection_retest_raw_temp(extractor);
    connectionRetestRaw_Up_Trans2D(v, :) = mri.connection_retest_raw_Up_trans(extractor);
    connectionRetestRaw_Low_Trans2D(v, :) = mri.connection_retest_raw_Low_trans(extractor);
    
    % Detectable connectivity 
    contrue_test_raw_temp = mri.contrue_test_raw_trans;
    contrueTestRawTrans(v, :, :) = contrue_test_raw_temp;
    contrueTestRaw_Trans2D(v, :) = contrue_test_raw_temp(extractor);
    contrueTestRaw_Up_Trans2D(v, :) = mri.contrue_test_raw_Up_trans(extractor);
    contrueTestRaw_Low_Trans2D(v, :) = mri.contrue_test_raw_Low_trans(extractor);
    
    contrue_retest_raw_temp = mri.contrue_retest_raw_trans;
    contrueRetestRawTrans(v, :, :) = contrue_retest_raw_temp;
    contrueRetestRaw_Trans2D(v, :) = contrue_retest_raw_temp(extractor);
    contrueRetestRaw_Up_Trans2D(v, :) = mri.contrue_retest_raw_Up_trans(extractor);
    contrueRetestRaw_Low_Trans2D(v, :) = mri.contrue_retest_raw_Low_trans(extractor);
    
    % Connectivity of data that were corrected for serial autocorrelations
    connection_test_remAutoRoi_temp = mri.connection_test_remAutoRoi_trans;
    connectionTestRemAutoTrans(v, :, :) = connection_test_remAutoRoi_temp;
    connectionTestRemAuto_Trans2D(v, :) = connection_test_remAutoRoi_temp(extractor);
    connectionTestRemAuto_Up_Trans2D(v, :) = mri.connection_test_remAutoRoi_Up_trans(extractor);
    connectionTestRemAuto_Low_Trans2D(v, :) = mri.connection_test_remAutoRoi_Low_trans(extractor);
    
    connection_retest_remAutoRoi_temp = mri.connection_retest_remAutoRoi_trans;
    connectionRetestRemAutoTrans(v, :, :) = connection_retest_remAutoRoi_temp;
    connectionRetestRemAuto_Trans2D(v, :) = connection_retest_remAutoRoi_temp(extractor);
    connectionRetestRemAuto_Up_Trans2D(v, :) = mri.connection_retest_remAutoRoi_Up_trans(extractor);
    connectionRetestRemAuto_Low_Trans2D(v, :) = mri.connection_retest_remAutoRoi_Low_trans(extractor);
    
    % Detectable connectivity of data that were corrected for serial autocorrelations
    contrue_test_remAutoRoi_temp = mri.contrue_test_remAutoRoi_trans;
    contrueTestRemAutoTrans(v, :, :) = contrue_test_remAutoRoi_temp;
    contrueTestRemAuto_Trans2D(v, :) = contrue_test_remAutoRoi_temp(extractor);
    contrueTestRemAuto_Up_Trans2D(v, :) = mri.contrue_test_remAutoRoi_Up_trans(extractor);
    contrueTestRemAuto_Low_Trans2D(v, :) = mri.contrue_test_remAutoRoi_Low_trans(extractor);
    
    contrue_retest_remAutoRoi_temp = mri.contrue_retest_remAutoRoi_trans;
    contrueRetestRemAutoTrans(v, :, :) = contrue_retest_remAutoRoi_temp;
    contrueRetestRemAuto_Trans2D(v, :) = contrue_retest_remAutoRoi_temp(extractor);
    contrueRetestRemAuto_Up_Trans2D(v, :) = mri.contrue_retest_remAutoRoi_Up_trans(extractor);
    contrueRetestRemAuto_Low_Trans2D(v, :) = mri.contrue_retest_remAutoRoi_Low_trans(extractor);

end

% Create average reliability of raw and AR(1) corrected time courses 
% Also provide confidence intervals for the sample mean and sample standard deviation
naiveRaw = average_rel(RelTrans, RelTransLow, RelTransUp, NumSimulation, 1);
naiveAR1 = average_rel(RelRemAutoTrans, RelRemAutoTransLow, RelRemAutoTransUp, NumSimulation, 1);
cleanRaw = average_rel(RelTrans, RelTransLow, RelTransUp, NumSimulation, 2);
cleanAR1 = average_rel(RelRemAutoTrans, RelRemAutoTransLow, RelRemAutoTransUp, NumSimulation, 2);

% Prepare the test-retest reliability table
rawTab = [struct2table(naiveRaw) struct2table(cleanRaw)];
AR1Tab = [struct2table(naiveAR1) struct2table(cleanAR1)];
RelAllTab = rows2vars(vertcat(rawTab, AR1Tab));

RelAllTab.Properties.VariableNames{2} = 'reliability';
RelAllTab.Properties.VariableNames{3} = 'reliabilityAR(1)';
% Write the mean reliability table (S1 table: Grand mean reliability)
writetable(RelAllTab, fullfile(path, 'RelMain.csv'));

% Prepare reliability per region table
relroitab = cell2table(fwf_region_name);
relroitab(:, 2) = table(tanh(mean(RelTrans, 2))); % Mean reliability per region
relroitab(:, 3) = table(tanh(std(RelTrans'))'); % Standard deviation of reliability per region
relroitab(:, 4) = table(tanh(mean(RelRemAutoTrans, 2))); % Mean reliability with AR(1) correction
relroitab(:, 5) = table(tanh(std(RelRemAutoTrans'))'); % Standard deviation with AR(1) correction

relroitab.Properties.VariableNames{2} = 'rel';
relroitab.Properties.VariableNames{3} = 'sample_std';
relroitab.Properties.VariableNames{4} = 'rel_Rem(AR1)';
relroitab.Properties.VariableNames{5} = 'sample_std_Rem(AR1)';
% Write reliability per region table (S2 Table: Reliability per brain region)
writetable(relroitab, fullfile(path, 'RelRoi.csv'));

% Create the correlation between brain reliability and behavior
% Reliability per brain region (not included in the paper)
RelBehBrain = cell2table(fwf_region_name);
RelBehBrain(:, 2) = table(corr((RelTrans)', atanh(relReact(1:nSamp)'))); % Correlation of raw reliability with behavior
RelBehBrain(:, 3) = table(corr((RelRemAutoTrans)', atanh(relReact(1:nSamp)'))); % Correlation of AR(1) corrected reliability with behavior
RelBehBrain.Properties.VariableNames{2} = 'rel';
RelBehBrain.Properties.VariableNames{3} = 'rel_Rem(AR1)';

% Set negative reliability values of response behavior to zero
relReact(relReact < 0) = 0;
mRel = mean(RelTrans); % Calculate the mean reliability across subjects
mRel(mRel < 0) = 0; % Set negative mean reliability values to zero
mRelRemA = mean(RelRemAutoTrans); % Calculate the mean reliability with AR(1) correction
mRelRemA(mRelRemA < 0) = 0; % Set negative mean reliability values with AR(1) correction to zero

% Create the correlation between mean brain reliability and behavior
% Reliability for censored reliability
RelBrainBeh = corr(mRel', atanh(relReact(1:nSamp)')); % Correlation of mean reliability with behavior
RelBrainAR_1Beh = corr(mRelRemA', atanh(relReact(1:nSamp)')); % Correlation of AR(1) corrected mean reliability with behavior
CorrBrainBeh = table(RelBrainBeh, RelBrainAR_1Beh); % Create a table for the correlations
% Write brain-behavior correlations to a CSV file
writetable(CorrBrainBeh, fullfile(path, 'CorrBrainBeh.csv'));

% Estimate the mean detectable connectivity as well as the confidence intervals of the sample mean and
% sample standard deviation for time series that were corrected for serial correlations (AR(1))
% and raw time series for a test and retest run
ConTest = struct2table(average_con(NumSimulation, connectionTestRaw_Trans2D, connectionTestRaw_Low_Trans2D, connectionTestRaw_Up_Trans2D));
ConRetest = struct2table(average_con(NumSimulation, connectionRetestRaw_Trans2D, connectionRetestRaw_Low_Trans2D, connectionRetestRaw_Up_Trans2D));
ConTestAR1 = struct2table(average_con(NumSimulation, connectionTestRemAuto_Trans2D, connectionTestRemAuto_Low_Trans2D, connectionTestRemAuto_Up_Trans2D));
ConRetestAR1 = struct2table(average_con(NumSimulation, connectionRetestRemAuto_Trans2D, connectionRetestRemAuto_Low_Trans2D, connectionRetestRemAuto_Up_Trans2D));
ConTrueTest = struct2table(average_con(NumSimulation, contrueTestRaw_Trans2D, contrueTestRaw_Low_Trans2D, contrueTestRaw_Up_Trans2D));
ConTrueRetest = struct2table(average_con(NumSimulation, contrueRetestRaw_Trans2D, contrueRetestRaw_Low_Trans2D, contrueRetestRaw_Up_Trans2D));
ConTrueTestAR1 = struct2table(average_con(NumSimulation, contrueTestRemAuto_Trans2D, contrueTestRemAuto_Low_Trans2D, contrueTestRemAuto_Up_Trans2D));
ConTrueRetestAR1 = struct2table(average_con(NumSimulation, contrueRetestRemAuto_Trans2D, contrueRetestRemAuto_Low_Trans2D, contrueRetestRemAuto_Up_Trans2D));

% Combine all connectivity tables into one
ConALLTab = vertcat(ConTest, ConRetest, ConTestAR1, ConRetestAR1, ConTrueTest, ConTrueRetest, ConTrueTestAR1, ConTrueRetestAR1);
% Create row names for the combined connectivity table
rowNames = {'ConnectivityTest', 'ConnectivityRetest', 'ConnectivityTestAR1', 'ConnectivityRetestAR1', 'DetectableConnectivityTest', 'DetectableConnectivityRetest', 'DetectableConnectivityTestAR1', 'DetectableConnectivityRetestAR1'};

% Add the row names to the combined connectivity table
ConALLTab.Properties.RowNames = rowNames;
% Write the main statistics table (Table 1 and S3 Table: Grand mean connectivity) to a CSV file
writetable(ConALLTab, fullfile(path, 'ConnectivityMain.csv'));

% Perform 2D t-tests to evaluate the connectivity paths of the connectomes against zero
% This will help estimate the Dice overlap in subsequent steps
[hTconRaw, pTconRaw, ciTconRaw, statsTconRaw] = ttest(connectionTestRaw_Trans2D);
[hRconRaw, pRconRaw, ciRconRaw, statsRconRaw] = ttest(connectionRetestRaw_Trans2D);

[hTcontRaw, pTcontRaw, ciTcontRaw, statsTcontRaw] = ttest(contrueTestRaw_Trans2D);
[hRcontRaw, pRcontRaw, ciRcontRaw, statsRcontRaw] = ttest(contrueRetestRaw_Trans2D);

[hTconRemAuto, pTconRemAuto, ciTconRemAuto, statsTconRemAuto] = ttest(connectionTestRemAuto_Trans2D);
[hRconRemAuto, pRconRemAuto, ciRconRemAuto, statsRconRemAuto] = ttest(connectionRetestRemAuto_Trans2D);

[hTcontRemAuto, pTcontRemAuto, ciTcontRemAuto, statsTcontRemAuto] = ttest(contrueTestRemAuto_Trans2D);
[hRcontRemAuto, pRcontRemAuto, ciRcontRemAuto, statsRcontRemAuto] = ttest(contrueRetestRemAuto_Trans2D);

% Perform conjunction/Dice overlap analysis for plotting purposes
[RawCon, RawDice] = Dice(pTconRaw, pRconRaw, l, 2);
[RawTCon, RawTDice] = Dice(pTcontRaw, pRcontRaw, l, 2);
[RemCon, RemDice] = Dice(pTconRemAuto, pRconRemAuto, l, 2);
[RemTCon, RemTDice] = Dice(pTcontRemAuto, pRcontRemAuto, l, 2);

% Conduct conjunction analysis over the t-statistics
minstatCraw = min(statsTconRaw.tstat, statsRconRaw.tstat); % Minimum t-statistics for raw connectivity
minstatCTraw = min(statsTcontRaw.tstat, statsRcontRaw.tstat); % Minimum t-statistics for detectable connectivity
minstatCremA = min(statsTconRemAuto.tstat, statsRconRemAuto.tstat); % Minimum t-statistics for AR(1) corrected connectivity
minstatCTremA = min(statsTcontRemAuto.tstat, statsRcontRemAuto.tstat); % Minimum t-statistics for AR(1) corrected detectable connectivity

% Create a matrix for the median t-statistics
Tdata = [minstatCraw; minstatCTraw; minstatCremA; minstatCTremA];

% Create a table for median t-statistics
rowNames = {'Connectivity', 'DetectableConnectivity', 'ConnectivityAR1', 'DetectableConnectivityAR1'};
medianTstat = table(median(Tdata, 2));
medianTstat.Properties.RowNames = rowNames;
% Write the median t-statistics table to a CSV file
writetable(medianTstat, fullfile(path, 'medianTstat.csv'));

% Estimate classic reliability across subjects
CorRaw = tanh(mean(atanh(diag(corr(connectionTestRaw_Trans2D, connectionRetestRaw_Trans2D)))));
CorRawTrue = tanh(mean(atanh(diag(corr(contrueTestRaw_Trans2D, contrueRetestRaw_Trans2D)))));
CorRemAuto = tanh(mean(atanh(diag(corr(connectionTestRemAuto_Trans2D, connectionRetestRemAuto_Trans2D)))));
CorRemAutoTrue = tanh(mean(atanh(diag(corr(contrueTestRemAuto_Trans2D, contrueRetestRemAuto_Trans2D)))));

% Create a table for between-subject reliability
BetweenSubRel = table(CorRaw, CorRemAuto);
writetable(BetweenSubRel, fullfile(path, 'BetweenSubRel.csv'));

% Perform conjunction analysis over the t-statistics
minRaw = min(connectionTestRaw_Trans2D, connectionRetestRaw_Trans2D);
minRawTrue = min(contrueTestRaw_Trans2D, contrueRetestRaw_Trans2D);
minRemAuto = min(connectionTestRemAuto_Trans2D, connectionRetestRemAuto_Trans2D);
minRemAutoTrue = min(contrueTestRemAuto_Trans2D, contrueRetestRemAuto_Trans2D);

% Compile minimum values into a data matrix
Cdata = [minRaw; minRawTrue; minRemAuto; minRemAutoTrue];


% Estimate the relationship between sample size, test-retest reliability,
% d(etectable) connectivity, and the amount of reproducibility that can be 
% obtained at the group level.
BO = bootRelSamplesize(connectionTestRaw_Trans2D, connectionRetestRaw_Trans2D, contrueTestRaw_Trans2D, contrueRetestRaw_Trans2D, connectionTestRemAuto_Trans2D,...
    connectionRetestRemAuto_Trans2D, contrueTestRemAuto_Trans2D, contrueRetestRemAuto_Trans2D, RelRemAutoTrans, RelTrans,NumSimulation,l);




