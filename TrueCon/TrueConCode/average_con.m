
%Copyright © 2024 Koten and Schüppen All rights reserved
%Important Notice: This code is not intended for medical applications 
%and does not have legal approval for such use. We strongly recommend 
%using FDA-approved software for any medical purposes. 


function [mcon] = average_con(nboot,Con,ConLow,ConUp)

% This function calculates the average of multiple confidence intervals and reverses the signs of both intervals
% if the upper bound is negative.

% Reshape mean and c.i. matrix
Con = reshape(Con,[],1);
ConLow = reshape(ConLow,[],1);
ConUp = reshape(ConUp,[],1);

% Average connectivity
mcon.meanRaw = tanh(nanmean(Con));

% Invert c.i. if upperbound is negative
index = ConUp < 0;
mcon.PercNeg = (sum(index) / length(Con)) * 100;
ConLow(index) = ConLow(index) * -1;
ConUp(index) = ConUp(index) * -1;
Con(index) = Con(index) * -1;

% Estimate mean of the connectivity and c.i.
mcon.mean = tanh(nanmean(Con));
mcon.Low = tanh(nanmean(ConLow));
mcon.Up = tanh(nanmean(ConUp));

% Estimate percentiles of the distribution
Scon = sort(Con);
n = length(Scon);
Li = round(2.5 / 100 * (n + 1));
Ui = round(97.5 / 100 * (n + 1));
mcon.SampLow = tanh(Scon(Li));
mcon.SampUp = tanh(Scon(Ui));

conf = sampboot(Con, nboot);
mcon.SampMeanLow = tanh(conf.meanLB);
mcon.SampMeanUp = tanh(conf.meanUB);

% Estimate std of distribution
mcon.STDSamp = tanh(std(Con));
mcon.STDSampLow = tanh(conf.stdLB);
mcon.STDSampUp = tanh(conf.stdUB);

end
