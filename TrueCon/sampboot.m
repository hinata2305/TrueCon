%Copyright © 2024 Koten and Schüppen All rights reserved
%Important Notice: This code is not intended for medical applications 
%and does not have legal approval for such use. We strongly recommend 
%using FDA-approved software for any medical purposes. 

function [conf] = sampboot(samp,numboot) 
% this function performs a bootstrap on a distribution 
% and estimates the C.I. for the mean and standard deviation

s=length(samp); % predefine matrix 
SampMean=zeros(numboot,1); 
SampSTD=zeros(numboot,1); 

% walk through loop for the number of bootstraps 

for bs = 1 : numboot 
% create random sample 
index=randi([1,s],s,1); 
SampMean(bs)=mean(samp(index)); 
SampSTD(bs)=std(samp(index)); 
end

 

% estimate confidence intervals
Rel=sort(SampMean);
n=length(Rel);
Li=round(2.5/100*(n+1)); 
Ui=round(97.5/100*(n+1)); 
conf.meanLB=Rel(Li);
conf.meanUB=Rel(Ui);

STD=sort(SampSTD);
conf.stdLB=STD(Li);
conf.stdUB=STD(Ui);