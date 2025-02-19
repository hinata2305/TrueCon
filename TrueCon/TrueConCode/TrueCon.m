
%Copyright © 2024 Koten and Schüppen All rights reserved
%Important Notice: This code is not intended for medical applications 
%and does not have legal approval for such use. We strongly recommend 
%using FDA-approved software for any medical purposes. 


function [reldifTest,contrueTest,reldifRetest,contrueRetest,relmax]=TrueCon(conTest,conRetest,relROIseed,relROItarget)

% Compute the sign of the original connectivity values for each element
signieTest=sign(conTest);
signieRetest=sign(conRetest);

% Take the absolute value of the connectivity values for each element
absconTest=abs(conTest);
absconRetest=abs(conRetest);

% Initialize output vectors
reldifTest = zeros(size(conTest));
contrueTest = zeros(size(conTest));
reldifRetest = zeros(size(conRetest));
contrueRetest = zeros(size(conRetest));
relmax = zeros(size(conTest));

% Loop through each element in the input vectors

for i = 1:size(conTest,1)

    % Check if both reliability coefficients are positive for the current element
    if (relROIseed(i)>0 & relROItarget(i)>0)
        % Estimate the upper bound of the detectable connectivity for the current element
        relmax(i) = (sqrt(relROIseed(i)*relROItarget(i)));
        % Calculate the absolute overestimation of the detectable connectivity for the current element
        reldifTest(i) = relmax(i)-absconTest(i);


        % Replace the absolute connectivity value with the detectable connectivity,
        % using the sign of the original connectivity data for the current element
        if absconTest(i) > relmax(i)
            contrueTest(i) = relmax(i).*signieTest(i);
        else
            contrueTest(i) = conTest(i);
        end

    else

        % Set the absolute overestimation, detectable connectivity, and upper bound to zero for the current element
        reldifTest(i) = 0-absconTest(i);
        contrueTest(i) = 0;
        relmax(i) = 0;
    end

    % Check if both reliability coefficients are positive for the current element
    if (relROIseed(i)>0 & relROItarget(i)>0)
        % Calculate the absolute overestimation of the detectable connectivity for the current element
        reldifRetest(i) = relmax(i)-absconRetest(i);


        % Replace the absolute connectivity value with the detectable connectivity,
        % using the sign of the original connectivity data for the current element
        if absconRetest(i) > relmax(i)
            contrueRetest(i) = relmax(i).*signieRetest(i);
        else
            contrueRetest(i) = conRetest(i);
        end

    else % Set the absolute overestimation, detectable connectivity, and upper bound to zero for the current element
        reldifRetest(i) = 0-absconRetest(i);
        contrueRetest(i) = 0;
    end

end
end