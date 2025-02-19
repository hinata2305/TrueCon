
%Copyright © 2024 Koten and Schüppen All rights reserved
%Important Notice: This code is not intended for medical applications 
%and does not have legal approval for such use. We strongly recommend 
%using FDA-approved software for any medical purposes. 




function [confidence] = ConfInterval(con)
% ConfInterval: Calculates the confidence interval for a given set of data
% con: Input data (should be a numerical array)
% confidence: Output structure with lower bound (LB) and upper bound (UB) of the confidence interval


% Remove any zero values from the data
index = con == 0;
con(index) = [];

% Sort the data in ascending order
ConSort = sort(con);

% Determine the number of data points
n = length(ConSort);

% Calculate the indices of the lower and upper bounds of the confidence interval
Li = round(2.5/100 * (n+1));
Ui = round(97.5/100 * (n+1));

% Try to set the lower and upper bounds of the confidence interval
try

    confidence.LB = ConSort(Li);
    confidence.UB = ConSort(Ui);

    % On an error (e.g., empty data), set lower and upper bounds to NaN
catch
    confidence.LB = nan;
    confidence.UB = nan;
end

end