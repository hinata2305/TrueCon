
clear all



%define functional data folder

addpath /data/backup/Graz/FWF/TRUE_CON/Code/
addpath /home/edvz/koten/opt/matlab/toolboxes/surfstat/
%addpath /bif/home/jan/linux
startdir=dir('/data/backup/Graz/FWF/SG_FILTER/DATA_MATLAB_REPRODUCTION/Logfiles/');
folName = '/data/backup/Graz/FWF/SG_FILTER/DATA_MATLAB_REPRODUCTION/FUNC/'
load namelist_n50smaller3mm_tasksmaller3wrong.mat
load ('mask.mat')


startdirAna='/data/backup/Graz/FWF/TRUE_CON/Code/';

%the size of the patch on the surface
diametre=8;
%length of timecourse
hemiS=487;
%number of timecourse in connectome
numTimecourses=34;
%number of bootstraps ot be
numBootstrap=10000;
%number of bootstraps
numSubject=length(namelist);
%size of connectome
numDatapoints=561;

%read surface
surf = SurfStatReadSurf({[startdirAna 'lh.sphere.reg']});
%read fmri mask 1=no data. find data that are empty
indexL=find((sum(maskRetestL,2)+sum(maskTestL,2))>1);
indexR=find((sum(maskRetestR,2)+sum(maskTestR,2))>1);

   s=size(maskTestL)
   ls=s(1);
   
    maskL_orig=ones(ls,1);
    maskR_orig=ones(ls,1);

    %The borders of the empty mask are dilated with roi size
 
    for mm=1:length(indexL)
        tempindex=find((SurfStatROI(indexL(mm),diametre.*2,surf)==1));
        maskL_orig(tempindex)=0;
   
    end
 
    for mm=1:length(indexR)
        tempindex=find((SurfStatROI(indexR(mm),diametre.*2,surf)==1));
        maskR_orig(tempindex)=0;
    end
 
 %here we define the regions for random selection. We concat the left and right hemisphere
 %as we want to conider them simultianiously
 MASKall=find(([maskL_orig ;maskR_orig])==1);
 SM=size(MASKall);
 %here we create the 3D random matrix
 RandIndices = randi(SM(1),numTimecourses,numBootstrap);

 clear maskR_orig maskL_orig
 %here we extract the relevant vertex numbers thatare realted to time courses     
 for i=1:numBootstrap
 RanCon(:,i)=MASKall(RandIndices(:,i));
 end


       %here we predefine the matrices for the output
      TempConTN=zeros(numBootstrap,numDatapoints,numSubject);
      TempConR=zeros(numBootstrap,numDatapoints,numSubject);
      TempTrueT=zeros(numBootstrap,numDatapoints,numSubject);
      TempTrueR=zeros(numBootstrap,numDatapoints,numSubject);
      TempConAR1T=zeros(numBootstrap,numDatapoints,numSubject);
      TempConAR1R=zeros(numBootstrap,numDatapoints,numSubject);
      TempTrueAR1T=zeros(numBootstrap,numDatapoints,numSubject);
      TempTrueAR1R=zeros(numBootstrap,numDatapoints,numSubject);
      MeanRelRoi=zeros(numBootstrap,numSubject);
      MeanRelAutoRoi=zeros(numBootstrap,numSubject);

%here is the subject lopp
for k = 1:numSubject

    disp k
    %read time course data of subject

    subject=namelist{k}
    %read in the data for the subject
    tic
    datafile=load([folName subject '/fulltimecourse.mat']);
    toc
    % extract reaction time data
    [testa,testb,encode1a,encode2a,encode1b,encode2b,reactST,reactSR,reactMT,reactMR,respST,respSR,respMT,respMR] = extract_rt_rel_response_filter(['/data/backup/Graz/FWF/SG_FILTER/DATA_MATLAB_REPRODUCTION/Logfiles/' subject]);

    % timecourse
    s1=size(datafile.time1);
    s2=size(datafile.time2);

    const1=s1(2)/2+1;
    const2=s2(2)/2+1;

    %here we extract the relevant sections of the timecourse sec1 test sec2
    %note that the timecourse of the left and right hemi are concatenated

    sec1=[encode2a(1) encode2a(1)+hemiS encode2a(1)+const1  encode2a(1)+hemiS+const1];
    sec2=[encode2b(1) encode2b(1)+hemiS encode2b(1)+const2  encode2b(1)+hemiS+const2];

    % gather the func data
    Test=[datafile.time1(:,sec1(1):sec1(2)); datafile.time1(:,sec1(3):sec1(4))];
    Retest=[datafile.time2(:,sec1(1):sec1(2)); datafile.time2(:,sec1(3):sec1(4))];

    % gather noise
    NoiseT=[datafile.noise_pca1(sec1(1):sec1(2),:) datafile.movpca1(sec1(1):sec1(2),1:2)];
    NoiseR=[datafile.noise_pca2(sec1(1):sec1(2),:) datafile.movpca2(sec1(1):sec1(2),1:2)];

    tl=size(Test);
  
tic
    parfor n=1:numBootstrap
       
      %select 34 vertex of interest
      SelectedVertex=RanCon(:,n);
      %extract timecourse and denoise it
      roiTest=preproROI(SelectedVertex,Test,NoiseT,diametre,surf);
      roiRetest=preproROI(SelectedVertex,Retest,NoiseR,diametre,surf);
      %estimate timecourse reliabilty and perform AR(1)
      TempRel=remove_auto1(roiTest, roiRetest);

      Trel=[TempRel.reltrans];
      TrelAuto=[TempRel.reltransAutoRem];
      %set negative reliabilty at zero
      Trel( Trel < 0) = 0;
      TrelAuto(TrelAuto< 0) = 0;
      MeanRelRoi(n,k)  =tanh(mean(Trel));
      MeanRelAutoRoi(n,k)  =tanh(mean(TrelAuto));
      
      % estimate all form of connectivity and perform AR(1)
      temp=remove_autoRoi(roiTest, roiRetest);

      TempConT(n,:,k)=temp.connection_test_raw_trans;
      TempConR(n,:,k)=temp.connection_retest_raw_trans;
      TempTrueT(n,:,k)=temp.contrue_test_raw_trans;
      TempTrueR(n,:,k)=temp.contrue_retest_raw_trans;
      TempConAR1T(n,:,k)=temp.connection_test_remAutoRoi_trans;
      TempConAR1R(n,:,k)=temp.connection_retest_remAutoRoi_trans;
      TempTrueAR1T(n,:,k)=temp.contrue_test_remAutoRoi_trans;
      TempTrueAR1R(n,:,k)=temp.contrue_retest_remAutoRoi_trans;


    end
toc

end

ConOverlapStat=DiceRoi(TempConT,TempConR);
ConAR1OverlapStat=DiceRoi(TempConAR1T,TempConAR1R);

TrueOverlapStat=DiceRoi(TempTrueT,TempTrueR);
TrueAR1OverlapStat=DiceRoi(TempTrueAR1T,TempTrueAR1R);

% Convert structure variables to tables
ConOverlapStatTable = struct2table(ConOverlapStat);
ConAR1OverlapStatTable = struct2table(ConAR1OverlapStat);
TrueOverlapStatTable = struct2table(TrueOverlapStat);
TrueAR1OverlapStatTable = struct2table(TrueAR1OverlapStat);

combinedTable = vertcat(ConOverlapStatTable, ConAR1OverlapStatTable, TrueOverlapStatTable, TrueAR1OverlapStatTable);
 
% Define the new row names
newRowNames = {'Observed', 'ObservedAR1', 'Detectable', 'DetectableAR1'};

% Add the new row names to the table
combinedTable.Properties.RowNames = newRowNames;

% Display the updated table
disp(combinedTable);

writetable(combinedTable,'/data/backup/Graz/FWF/TRUE_CON/Stat_Result/BootstrapConnectome.csv')  


% Displaying the table
