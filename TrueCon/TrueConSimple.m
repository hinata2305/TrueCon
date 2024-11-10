%Copyright © 2024 Koten and Schüppen All rights reserved
%Important Notice: This code is not intended for medical applications 
%and does not have legal approval for such use. We strongly recommend 
%using FDA-approved software for any medical purposes. 

function [reldifTest,contrueTest,reldifRetest,contrueRetest,relmax]=TrueConSimple(conTest,conRetest,relROIseed,relROItarget)

% Compute the sign of the original connectivity values 
signieTest=sign(conTest);
signieRetest=sign(conRetest);
% Take the absolute value of the connectivity values
absconTest=abs(conTest); absconRetest=abs(conRetest);

% Check if both reliability coefficients are positive
if (relROIseed>0 && relROItarget>0)
    % Estimate the upper bound of the detectable connectivity
    relmax=(sqrt(relROIseed*relROItarget));
    % Calculate the absolute overestimation of the detectable connectivity
    reldifTest=relmax-absconTest;
    % Replace the absolute connectivity values with the detectable connectivity,
    % using the sign of the original connectivity data
    if absconTest > relmax
        contrueTest=relmax.*signieTest;
    else
        contrueTest= conTest;
    end

else
    % Set the absolute overestimation, detectable connectivity, and upper bound to zero
    reldifTest=0-absconTest;
    contrueTest=0;
    relmax=0;
end

% Check if both reliability coefficients are positive
if (relROIseed>0 && relROItarget>0)
    % Calculate the absolute overestimation of the detectable connectivity
    reldifRetest=relmax-absconRetest;


    % Replace the absolute connectivity values with the detectable connectivity,
    % using the sign of the original connectivity data
    if absconRetest > relmax
        contrueRetest=relmax.*signieRetest;
    else
        contrueRetest= conRetest;
    end

else
    % Set the absolute overestimation, detectable connectivity, and upper bound to zero
    reldifRetest=0-absconTest;
    contrueRetest=0;
end

end




