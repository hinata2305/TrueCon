% Warranty this code has not undergone peer code review yet. We do not take responsibility for its correctness.
% This code may not be used for medical applications. ©Koten Schüppen
% 25/01/2024

% This function creates first-level statistics for fMRI reliability and connectivity studies and analyzes how reproducible group studies are.
% This function estimates TSNR, autocorrelation of time course, connectivity, and true (detectable) connectivity.

% add the path where the 
current = 
addpath current

load fMRI_clean.mat
load relResponseTime.mat
load fwf_region_name.mat
load namelist_n50smaller3mm_tasksmaller3wrong.mat
load StroopCorrectAwnser.mat
load MemoryCorrectAwnser.mat



% size of the connectivity matrix
s = [34 34];
% matrix used to extract lower triangle of connectivity data
extractor = tril(true(s(2)),-1);

% threshold span for dice overlap
l = 250;


numboot=10000;
NumSimulation=10000;

%the subject loop
parfor v=1:50

    % Name of subject
    subject=namelist{v}
    % read in timcourses from a test and retest run per subject
    timesec1=Clean_SPM_p_1(:,:,v);
    timesec2=Clean_SPM_p_2(:,:,v);
    % estaimte test retest reliabilty with C.I. from raw timecourses
    % and timecoutses that were correct for serial correaltions
    rmi=remove_AutoBlockBoot_short(timesec1,timesec2,numboot);

    RelTrans(:,v)=rmi.relRoi_raw_trans;
    RelTransUp(:,v)=rmi.relRoi_raw_Up_trans;
    RelTransLow(:,v)=rmi.relRoi_raw_Low_trans;

    RelRemAutoTrans(:,v)=rmi.relRoi_remAuto_trans;
    RelRemAutoTransUp(:,v)=rmi.relRoi_remAuto_Up_trans;
    RelRemAutoTransLow(:,v)=rmi.relRoi_remAuto_Low_trans;

    tic
    mri = remove_AutoBlockBoot(timesec1,timesec2,numboot);
    toc

    % estimate (true=detectable) connectivity with C.I. from raw timecourses
    % and timecourses that were corrected for serial correlations. For
    % convience some matrice are tranformed into 1D and gathered into 2D

    connection_test_raw_temp=mri.connection_test_raw_trans;
    connectionTestRawTrans(v,:,:)=connection_test_raw_temp;
    connectionTestRaw_Trans2D(v,:)=connection_test_raw_temp(extractor);
    connectionTestRaw_Up_Trans2D(v,:)=mri.connection_test_raw_Up_trans(extractor)
    connectionTestRaw_Low_Trans2D(v,:)=mri.connection_test_raw_Low_trans(extractor)

    connection_retest_raw_temp=mri.connection_retest_raw_trans;
    connectionRetestRawTrans(v,:,:)=connection_retest_raw_temp;
    connectionRetestRaw_Trans2D(v,:)=connection_retest_raw_temp(extractor);
    connectionRetestRaw_Up_Trans2D(v,:)=mri.connection_retest_raw_Up_trans(extractor)
    connectionRetestRaw_Low_Trans2D(v,:)=mri.connection_retest_raw_Low_trans(extractor)

    contrue_test_raw_temp=mri.contrue_test_raw_trans;
    contrueTestRawTrans(v,:,:)=contrue_test_raw_temp;
    contrueTestRaw_Trans2D(v,:)=contrue_test_raw_temp(extractor);
    contrueTestRaw_Up_Trans2D(v,:)=mri.contrue_test_raw_Up_trans(extractor);
    contrueTestRaw_Low_Trans2D(v,:)=mri.contrue_test_raw_Low_trans(extractor);

    contrue_retest_raw_temp=mri.contrue_retest_raw_trans;
    contrueRetestRawTrans(v,:,:)=contrue_retest_raw_temp;
    contrueRetestRaw_Trans2D(v,:)=contrue_retest_raw_temp(extractor);
    contrueRetestRaw_Up_Trans2D(v,:)=mri.contrue_retest_raw_Up_trans(extractor);
    contrueRetestRaw_Low_Trans2D(v,:)=mri.contrue_retest_raw_Low_trans(extractor);

    connection_test_remAutoRoi_temp=mri.connection_test_remAutoRoi_trans;
    connectionTestRemAutoTrans(v,:,:)=connection_test_remAutoRoi_temp;
    connectionTestRemAuto_Trans2D(v,:)=connection_test_remAutoRoi_temp(extractor);
    connectionTestRemAuto_Up_Trans2D(v,:)=mri.connection_test_remAutoRoi_Up_trans(extractor);
    connectionTestRemAuto_Low_Trans2D(v,:)=mri.connection_test_remAutoRoi_Low_trans(extractor);

    connection_retest_remAutoRoi_temp=mri.connection_retest_remAutoRoi_trans;
    connectionRetestRemAutoTrans(v,:,:)=connection_retest_remAutoRoi_temp;
    connectionRetestRemAuto_Trans2D(v,:)=connection_retest_remAutoRoi_temp(extractor);
    connectionRetestRemAuto_Up_Trans2D(v,:)=mri.connection_retest_remAutoRoi_Up_trans(extractor);
    connectionRetestRemAuto_Low_Trans2D(v,:)=mri.connection_retest_remAutoRoi_Low_trans(extractor);

    contrue_test_remAutoRoi_temp=mri.contrue_test_remAutoRoi_trans;
    contrueTestRemAutoTrans(v,:,:)=contrue_test_remAutoRoi_temp;
    contrueTestRemAuto_Trans2D(v,:)=contrue_test_remAutoRoi_temp(extractor);
    contrueTestRemAuto_Up_Trans2D(v,:)=mri.contrue_test_remAutoRoi_Up_trans(extractor);
    contrueTestRemAuto_Low_Trans2D(v,:)=mri.contrue_test_remAutoRoi_Low_trans(extractor);

    contrue_retest_remAutoRoi_temp=mri.contrue_retest_remAutoRoi_trans;
    contrueRetestRemAutoTrans(v,:,:)=contrue_retest_remAutoRoi_temp;
    contrueRetestRemAuto_Trans2D(v,:)=contrue_retest_remAutoRoi_temp(extractor);
    contrueRetestRemAuto_Up_Trans2D(v,:)=mri.contrue_retest_remAutoRoi_Up_trans(extractor);
    contrueRetestRemAuto_Low_Trans2D(v,:)=mri.contrue_retest_remAutoRoi_Low_trans(extractor);

end

sampN=50;

%create average rel
naiveRaw=average_rel(RelTrans,RelTransLow,RelTransUp,NumSimulation,1);
naiveAR1=average_rel(RelRemAutoTrans,RelRemAutoTransLow,RelRemAutoTransUp,NumSimulation,1);
cleanRaw=average_rel(RelTrans,RelTransLow,RelTransUp,NumSimulation,2);
cleanAR1=average_rel(RelRemAutoTrans,RelRemAutoTransLow,RelRemAutoTransUp,NumSimulation,2);

rawTab=([struct2table(naiveRaw) struct2table(cleanRaw)]);
AR1Tab=([struct2table(naiveAR1) struct2table(cleanAR1)]);
RelAllTab=rows2vars(vertcat(rawTab, AR1Tab));

RelAllTab.Properties.VariableNames{2} = 'reliability';
RelAllTab.Properties.VariableNames{3} = 'reliabilityAR(1)';

writetable(RelAllTab,'/data/backup/Graz/FWF/TRUE_CON/Stat_Result/RelMain.csv')

relroitab=cell2table(fwf_region_name)
relroitab(:,2)=table(tanh(mean(RelTrans,2)))
relroitab(:,3)=table(tanh(std(RelTrans'))')
relroitab(:,4)=table(tanh(mean(RelRemAutoTrans,2)))
relroitab(:,5)=table(tanh(std(RelRemAutoTrans'))')

relroitab.Properties.VariableNames{2} = 'rel';
relroitab.Properties.VariableNames{3} = 'sample_std';
relroitab.Properties.VariableNames{4} = 'rel_Rem(AR1)';
relroitab.Properties.VariableNames{5} = 'sample_std_Rem(AR1)';

writetable(relroitab,'/data/backup/Graz/FWF/TRUE_CON/Stat_Result/RelRoi.csv')

%here we create the correaltion between brain reliability and behavior
%reliabilty

relReact=(diag(corr(MixRespTest',MixRespRetest','rows','pairwise')));
RelBehBrain=cell2table(fwf_region_name);
RelBehBrain(:,2) = table(corr((RelTrans)',atanh(relReact)));
RelBehBrain(:,3) = table(corr((RelRemAutoTrans)',atanh(relReact)));
RelBehBrain.Properties.VariableNames{2} = 'rel';
RelBehBrain.Properties.VariableNames{3} = 'rel_Rem(AR1)';

relReact(relReact<0)=0;
mRel=mean(RelTrans);
mRel(mRel<0)=0;
mRelRemA=mean(RelRemAutoTrans);
mRelRemA(mRelRemA<0)=0;
RelBrainBeh=corr(mRel',atanh(relReact));
RelBrainAR_1Beh=corr(mRelRemA',atanh(relReact));
CorrBrainBeh=table(RelBrainBeh,RelBrainAR_1Beh);

writetable(CorrBrainBeh,'/data/backup/Graz/FWF/TRUE_CON/Stat_Result/CorrBrainBeh.csv')

ConTest=struct2table(average_con(NumSimulation,connectionTestRaw_Trans2D,connectionTestRaw_Low_Trans2D,connectionTestRaw_Up_Trans2D));
ConRetest=struct2table(average_con(NumSimulation,connectionRetestRaw_Trans2D,connectionRetestRaw_Low_Trans2D,connectionRetestRaw_Up_Trans2D));
ConTestAR1=struct2table(average_con(NumSimulation,connectionTestRemAuto_Trans2D,connectionTestRemAuto_Low_Trans2D,connectionTestRemAuto_Up_Trans2D));
ConRetestAR1=struct2table(average_con(NumSimulation,connectionRetestRemAuto_Trans2D,connectionRetestRemAuto_Low_Trans2D,connectionRetestRemAuto_Up_Trans2D));
ConTrueTest=struct2table(average_con(NumSimulation,contrueTestRaw_Trans2D,contrueTestRaw_Low_Trans2D,contrueTestRaw_Up_Trans2D));
ConTrueRetest=struct2table(average_con(NumSimulation,contrueRetestRaw_Trans2D,contrueRetestRaw_Low_Trans2D,contrueRetestRaw_Up_Trans2D));
ConTrueTestAR1=struct2table(average_con(NumSimulation,contrueTestRemAuto_Trans2D,contrueTestRemAuto_Low_Trans2D,contrueTestRemAuto_Up_Trans2D));
ConTrueRetestAR1=struct2table(average_con(NumSimulation,contrueRetestRemAuto_Trans2D,contrueRetestRemAuto_Low_Trans2D,contrueRetestRemAuto_Up_Trans2D));

ConALLTab= vertcat(ConTest,ConRetest,ConTestAR1,ConRetestAR1,ConTrueTest,ConTrueRetest,ConTrueTestAR1,ConTrueRetestAR1);
% Create the row names
rowNames = {'ConnectivityTest', 'ConnectivityRetest', 'ConnectivityTestAR1', 'ConnectivityRetestAR1', 'DetectableConnectivityTest', 'DetectableConnectivityRetest', 'DetectableConnectivityTestAR1', 'DetectableConnectivityRetestAR1'};

% Add the row names
ConALLTab.Properties.RowNames = rowNames;

writetable(ConALLTab,'/data/backup/Graz/FWF/TRUE_CON/Stat_Result/ConnectivityMain.csv')

%perform 2D ttest
[hTconRaw,pTconRaw,ciTconRaw,statsTconRaw] = ttest(connectionTestRaw_Trans2D);
[hRconRaw,pRconRaw,ciRconRaw,statsRconRaw]= ttest(connectionRetestRaw_Trans2D);

[hTcontRaw,pTcontRaw,ciTcontRaw,statsTcontRaw]  = ttest(contrueTestRaw_Trans2D);
[hRcontRaw,pRcontRaw,ciRcontRaw,statsRcontRaw] = ttest(contrueRetestRaw_Trans2D);

[hTconRemAuto,pTconRemAuto,ciTconRemAuto,statsTconRemAuto] = ttest(connectionTestRemAuto_Trans2D);
[hRconRemAuto,pRconRemAuto,ciRconRemAuto,statsRconRemAuto] = ttest(connectionRetestRemAuto_Trans2D);

[hTcontRemAuto,pTcontRemAuto,ciTcontRemAuto,statsTcontRemAuto] = ttest(contrueTestRemAuto_Trans2D);
[hRcontRemAuto,pRcontRemAuto,ciRcontRemAuto,statsRcontRemAuto] = ttest(contrueRetestRemAuto_Trans2D);

%here we perform the conjunction/dice overlap for plotting purposes
[RawCon,RawDice]=Dice(pTconRaw,pRconRaw,l,2);
[RawTCon,RawTDice]=Dice(pTcontRaw,pRcontRaw,l,2);
[RemCon,RemDice]=Dice(pTconRemAuto,pRconRemAuto,l,2);
[RemTCon,RemTDice]=Dice(pTcontRemAuto,pRcontRemAuto,l,2);

%here we perform conjunction analysis over the t-statistics
minstatCraw=min(statsTconRaw.tstat,statsRconRaw.tstat);
minstatCTraw=min(statsTcontRaw.tstat,statsRcontRaw.tstat);
minstatCremA=min(statsTconRemAuto.tstat,statsRconRemAuto.tstat);
minstatCTremA=min(statsTcontRemAuto.tstat,statsRcontRemAuto.tstat);

Tdata=[minstatCraw; minstatCTraw; minstatCremA; minstatCTremA];

rowNames = {'Connectivity','DetectableConnectivity','ConnectivityAR1',  'DetectableConnectivityAR1'};
medianTstat=table(median(Tdata,2));
medianTstat.Properties.RowNames=rowNames
writetable(medianTstat,'/data/backup/Graz/FWF/TRUE_CON/Stat_Result/medianTstat.csv')

%here we estimate classic rel over subjects
CorRaw=tanh(mean(atanh(diag(corr(connectionTestRaw_Trans2D,connectionRetestRaw_Trans2D)))));
CorRawTrue=tanh(mean(atanh(diag(corr(contrueTestRaw_Trans2D,contrueRetestRaw_Trans2D)))));
CorRemAuto=tanh(mean(atanh(diag(corr(connectionTestRemAuto_Trans2D,connectionRetestRemAuto_Trans2D)))));
CorRemAutoTrue=tanh(mean(atanh(diag(corr(contrueTestRemAuto_Trans2D,contrueRetestRemAuto_Trans2D)))));

BetweenSubRel = table(CorRaw, CorRemAuto);
writetable(BetweenSubRel,'/data/backup/Graz/FWF/TRUE_CON/Stat_Result/BetweenSubRel.csv')

minRaw=min( (connectionTestRaw_Trans2D),  (connectionRetestRaw_Trans2D));
minRawTrue=min( (contrueTestRaw_Trans2D),  (contrueRetestRaw_Trans2D));
minRemAuto=min( (connectionTestRemAuto_Trans2D),  (connectionRetestRemAuto_Trans2D));
minRemAutoTrue=min( (contrueTestRemAuto_Trans2D),  (contrueTestRemAuto_Trans2D));

Cdata=[minRaw; minRawTrue; minRemAuto; minRemAutoTrue];




[fList,pList] = matlab.codetools.requiredFilesAndProducts('Estimate_ResAuto_TrueConBlockBoot.m');


newFolder = '/data/backup/Graz/FWF/TRUE_CON/Code/';



for i = 1:numel(fList)

    copyfile( fList{i}, newFolder);

end

