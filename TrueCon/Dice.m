
%Copyright © 2024 Koten and Schüppen All rights reserved
%Important Notice: This code is not intended for medical applications 
%and does not have legal approval for such use. We strongly recommend 
%using FDA-approved software for any medical purposes. 


function [Con,DiceOverlap] = Dice(ImageTest,ImageRetest,l,version)

% This function estimates the conjunction (Con) and dice overlap from a test and retest image for a larger range of thresholds given in l. %
% Version 1 uses a linear threshold that looks for values larger than the threshold, while version 2 uses a power-based threshold that looks for values smaller than the threshold. %


% Check which version is being used %

if version == 1

    % Initialize Con and DiceOverlap with NaN values
    Con=nan(l,1);
    DiceOverlap=nan(l,1);

    % Iterate through the thresholds
    for z = 1:l

        % Calculate the threshold value
        t=z*0.2;

        % Calculate the Con and DiceOverlap values based on the threshold
        Con(z) = sum(ImageTest>t & ImageRetest>t);
        DiceOverlap(z) = (2.*(sum(ImageTest>t & ImageRetest>t))) / (sum(ImageTest>t) + sum(ImageRetest>t));

    end

else


    % Initialize Con and DiceOverlap with NaN values
    Con=nan(l,1);
    DiceOverlap=nan(l,1);

    % Iterate through the thresholds
    for z = 1:l

        % Calculate the h and t values for the power-based threshold
        h=z*.1;
        t= power(0.05,h);

        % Calculate the Con and DiceOverlap values based on the threshold
        Con(z) = sum(ImageTest<t & ImageRetest<t);
        DiceOverlap(z) = (2.*(sum(ImageTest<t & ImageRetest<t))) / (sum(ImageTest<t) + sum(ImageRetest<t));

    end

end