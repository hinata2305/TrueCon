
%Copyright © 2024 Koten and Schüppen All rights reserved
%Important Notice: This code is not intended for medical applications 
%and does not have legal approval for such use. We strongly recommend 
%using FDA-approved software for any medical purposes. 


% This function calculates the average reliability estimates (rel) and also averages the confidence intervals of these estimates (RelLow RelUp). Additionally, it provides the percentiles of the distributions (Li Ui).
% In version 1, the mean is calculated considering both negative and positive reliability values. However, in version 2, the negative reliability is either set at zero or excluded from the matrix altogether.

function [rel] = average_rel(Rel,RelLow,RelUp,boot,version)

% get the size of the input matrices 
s=size(Rel); 
total=s(1)*s(2);

if version==1


    % reshape the matrices into column vectors
    Rel=reshape(Rel,[],1);
    RelLow=reshape(RelLow,[],1);
    RelUp=reshape(RelUp,[],1);

    % calculate the mean reliability estimates
    rel.Mean=tanh(mean(Rel));
    % calculate the mean confidence intervals for the time course
    rel.TimecourseLow=tanh(mean(RelLow));
    rel.TimecourseUp=tanh(mean(RelUp));

    % sort the reliability estimates in ascending order
    Srel=sort(Rel);
    n=length(Srel);
    % calculate the percentiles of the distributions
    Li=round(2.5./100.*(n+1));
    Ui=round(97.5./100.*(n+1));

    % transform the percentiles of the sample distributions
    rel.SampleLow=tanh(Srel(Li));
    rel.SampleUp=tanh(Srel(Ui));
    % transform the standard deviation of the sample distributions
    rel.SampleSTD=tanh(std(Srel));

    % perform bootstrapping to estimate the confidence intervals of the mean
    conf=sampboot(reshape(Rel,[],1),boot);
    rel.SampleMeanLow=tanh(conf.meanLB);
    rel.SampleMeanUp=tanh(conf.meanUB);

    % calculate the percentage of positive, negative, and mixed c.i. estimates in the time course
    indexPos=RelLow>0&RelUp>0;
    indexNeg=RelLow<0&RelUp<0;
    indexMix=RelLow<0&RelUp>0;

    rel.PosTimecourseSum=((sum(indexPos))/total)*100;
    % calculate the mean confidence intervals for the positive reliability estimates
    rel.PosTimecourseLow=tanh(mean(RelLow(indexPos)));
    rel.PosTimecourseUp=tanh(mean(RelUp(indexPos)));
    rel.NegTimecourseSum=((sum(indexNeg))/total)*100;
    % calculate the mean confidence intervals for the negative reliability estimates
    rel.NegTimecourseLow=tanh(mean(RelLow(indexNeg)));
    rel.NegTimecourseUp=tanh(mean(RelUp(indexNeg)));
    rel.MixTimecourseSum=((sum(indexMix))/total)*100;
    % calculate the mean confidence intervals for the mixed reliability estimates
    rel.MixTimecourseLow=tanh(mean(RelLow(indexMix)));
    rel.MixTimecourseUp=tanh(mean(RelUp(indexMix)));

else


    % reshape the matrices into column vectors
    Rel=reshape(Rel,[],1);
    indexR=Rel<0;

    % set the negative reliability estimates to zero
    Rel(indexR)=0;
    % calculate the mean reliability estimates
    rel.ZeroMean=tanh(mean(Rel));

    % sort the non-negative reliability estimates in ascending order
    Srel=sort(Rel);
    n=length(Srel);
    % calculate the percentiles of the non-negative sample distributions
    Li=round(2.5./100.*(n+1));
    Ui=round(97.5./100.*(n+1));
    % calculate the percentiles of the non-negative sample distributions
    rel.ZeroSampleLow=tanh(Srel(Li));
    rel.ZeroSampleUp=tanh(Srel(Ui));
    % calculate the standard deviation of the non-negative sample distributions
    rel.ZeroSampleSTD=tanh(std(Srel));
    % perform bootstrapping to estimate the confidence intervals of the mean
    conf=sampboot(Rel,boot);
    rel.ZeroSampleMeanLow=tanh(conf.meanLB);
    rel.ZeroSampleMeanUp=tanh(conf.meanUB);

    % remove the negative reliability estimates from the matrices
    Rel(indexR)=[];
    RelLow(indexR)=[];
    RelUp(indexR)=[];

    % calculate the percentage of corrupted samples
    rel.corrupt=((sum(indexR>0))/total).*100;
    % calculate the mean reliability estimates for the reliability estimates
    rel.NoZeroMean=tanh(mean(mean(Rel)));
    % calculate the mean confidence intervals for the time course
    rel.NoZeroTimecourseLow=tanh(mean(mean(RelLow)));
    rel.NoZeroTimecourseUp=tanh(mean(mean(RelUp)));
    % sort the  reliability estimates in ascending order
    Srel=sort(Rel);
    n=length(Srel);
    % calculate the percentiles of the  sample distributions
    Li=round(2.5./100.*(n+1));
    Ui=round(97.5./100.*(n+1));
    % calculate the percentiles of the sample distributions
    rel.NoZeroSampleLow=tanh(Srel(Li));
    rel.NoZeroSampleUp=tanh(Srel(Ui));
    % calculate the standard deviation of the sample distributions
    rel.NoZeroSampSTD=tanh(std(Srel));
    % perform bootstrapping to estimate the confidence intervals of th mean
    conf=sampboot(Rel,boot);
    rel.NoZeroSampleMeanLow=tanh(conf.meanLB);
    rel.NoZeroSampleMeanUp=tanh(conf.meanUB);

end

end
