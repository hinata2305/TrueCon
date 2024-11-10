%Copyright © 2024 Koten and Schüppen All rights reserved
%Important Notice: This code is not intended for medical applications 
%and does not have legal approval for such use. We strongly recommend 
%using FDA-approved software for any medical purposes. 

function [remAuto1,remAuto2,zeta,residual]= autocorfilter(ROIseed,ROItarget)

% Function to remove autocorrelation from a time series

% Step 1: Extract residual time course
[~,~,residual]=regress(ROIseed,[ROItarget ones(length(ROItarget),1)]);

% Step 2: Estimate autocorrelation of residual
zeta=regress(residual(1:end-1),[residual(2:end) ones(length(residual)-1,1)] ); zeta= zeta(1);

% Step 3: Remove residual autocorrelation from time course
remAuto1 = [ROIseed(1:end-1)-(zeta.*ROIseed(2:end));ROIseed(end)];
remAuto2 = [ROItarget(1:end-1)-(zeta.*ROItarget(2:end));ROItarget(end)];

end
