
%Copyright © 2024 Koten and Schüppen All rights reserved
%Important Notice: This code is not intended for medical applications 
%and does not have legal approval for such use. We strongly recommend 
%using FDA-approved software for any medical purposes. 


% This routine aims to estimate the impact of sample size on the reproducibility of fMRI group studies and how the reproducibility of group results is influenced by the reproducibility of the time course.
% The routine requires the following input, which is obtained from the Estimate_ResAuto_TrueCon4Boot function:

%connectionTestRaw_Trans2D: The connectivity of test data for the initial dataset
%connectionRetestRaw_Trans2D: The connectivity of retest data for the initial dataset
%contrueTestRaw_Trans2D: The detectable connectivity for test data constrained by the test-retest reliability (RelTrans)
%contrueRetestRaw_Trans2D: The detectable connectivity for retest data
%connectionTestRemAuto_Trans2D: The connectivity with serial auto correlations removed for the test dataset
%connectionRetestRemAuto_Trans2D: The connectivity with serial auto correlations removed for the retest dataset
%contrueTestRemAuto_Trans2D: The detectable connectivity test constrained by the test-retest reliability (RelRemAutoTransTrans), when serial auto correlations are removed
%contrueRetestRemAuto_Trans2D: The detectable connectivity retest after removing serial auto correlations
%RelRemAutoTrans: The test-retest reliability when serial auto correlations are removed
%RelTrans: The test-retest reliability

% The routine iterates through a loop to increase the sample size and stores relevant information about the sample size in the variable "z". It then runs an inner loop to perform N bootstraps specified in the variable "NumSimulation".

NumSimulation = 10000;

for z = 1:41


    % Increase sample size
    samplesize = (9) + 1 * z;

    % Initialize matrices to store results
    ConjRaw = nan(NumSimulation, l);
    DiceRaw = nan(NumSimulation, l);
    ConjRawTrue = nan(NumSimulation, l);
    DiceRawTrue = nan(NumSimulation, l);
    ConjRauto = nan(NumSimulation, l);
    DiceRauto = nan(NumSimulation, l);
    ConjRautoTrue = nan(NumSimulation, l);
    DiceRautoTrue = nan(NumSimulation, l);
    MeanConnectivityTrans = nan(NumSimulation, 4);
    MeanContrueTrans = nan(NumSimulation, 4);
    meanRelAutoTrans = nan(NumSimulation, 1);
    meanRelTrans = nan(NumSimulation, 1);

    tic

    parfor ss = 1:NumSimulation

        % Randomly select samples
        v = randperm(50, samplesize);

        % Subset the data based on the selected samples
        connectionTestRawTrans_smp_2D = connectionTestRaw_Trans2D(v, :);
        connectionRetestRaw_smp_Trans2D = connectionRetestRaw_Trans2D(v, :);
        contrueTestRaw_smp_Trans2D = contrueTestRaw_Trans2D(v, :);
        contrueRetestRaw_smp_Trans2D = contrueRetestRaw_Trans2D(v, :);
        connectionTestRemAuto_smp_Trans2D = connectionTestRemAuto_Trans2D(v, :);
        connectionRetestRemAuto_smp_Trans2D = connectionRetestRemAuto_Trans2D(v, :);
        contrueTestRemAuto_smp_Trans2D = contrueTestRemAuto_Trans2D(v, :);
        contrueRetestRemAuto_smp_Trans2D = contrueRetestRemAuto_Trans2D(v, :);

        % Calculate mean of RelRemAutoTrans and RelTrans
        meanRelAuto_smp_Trans(ss) = mean(mean(RelRemAutoTrans(:, v)));
        meanRel_smp_Trans(ss) = mean(mean(RelTrans(:, v)));

        % Perform t-tests and calculate Dice coefficient
        [~, p_test] = ttest(connectionTestRawTrans_smp_2D);
        [~, p_retest] = ttest(connectionRetestRaw_smp_Trans2D);
        [ConjRaw(ss, :), DiceRaw(ss, :)] = Dice(p_test, p_retest, l, 2);

        [~, p_test] = ttest(contrueTestRaw_smp_Trans2D);
        [~, p_retest] = ttest(contrueRetestRaw_smp_Trans2D);
        [ConjRawTrue(ss, :), DiceRawTrue(ss, :)] = Dice(p_test, p_retest, l, 2);

        [~, p_test] = ttest(connectionTestRemAuto_smp_Trans2D);
        [~, p_retest] = ttest(connectionRetestRemAuto_smp_Trans2D);
        [ConjRauto(ss, :), DiceRauto(ss, :)] = Dice(p_test, p_retest, l, 2);

        [~, p_test] = ttest(contrueTestRemAuto_smp_Trans2D);
        [~, p_retest] = ttest(contrueRetestRemAuto_smp_Trans2D);
        [ConjRautoTrue(ss, :), DiceRautoTrue(ss, :)] = Dice(p_test, p_retest, l, 2);

        % Estimate connectivity and true connectivity, with and without serial autocorrelations
        MeanConnectivityTrans(ss, :) = [(mean(mean(connectionTestRawTrans_smp_2D))) (mean(mean(connectionRetestRaw_smp_Trans2D))) (mean(mean(connectionTestRemAuto_smp_Trans2D))) (mean(mean(connectionRetestRemAuto_smp_Trans2D)))];
        MeanContrueTrans(ss, :) = [(mean(mean(contrueTestRaw_smp_Trans2D))) (mean(mean(contrueRetestRaw_smp_Trans2D))) (mean(mean(contrueTestRemAuto_smp_Trans2D))) (mean(mean(contrueRetestRemAuto_smp_Trans2D)))];

    end

    toc

    % Calculate 2.5th and 97.5th percentiles
    Left5 = round((NumSimulation ./ 100) * 2.5);
    Right5 = round(NumSimulation - (NumSimulation ./ 100) * 2.5);

    % Sort and store connectivity and Dice coefficients
    ConjAll = reshape([ConjRaw ConjRawTrue ConjRauto ConjRautoTrue], [NumSimulation, l, 4]);
    MConjALL(z, :, :) = squeeze(mean(ConjAll));
    STDConjALL(z, :, :) = squeeze(std(ConjAll));
    ConjAll = sort(ConjAll, 2, 'descend');
    Left5_ConjALL(z, :, :) = squeeze(ConjAll(Left5, :, :));
    Right5_ConjALL(z, :, :) = squeeze(ConjAll(Right5, :, :));

    DiceAll = reshape([DiceRaw DiceRawTrue DiceRauto DiceRautoTrue], [NumSimulation, l, 4]);
    MDiceALL(z, :, :) = squeeze(nanmean(DiceAll));
    DiceNannie(z, :, :) = squeeze(sum(isnan(DiceAll)));
    STDDiceALL(z, :, :) = squeeze(nanstd(DiceAll));
    DiceAll = sort(DiceAll, 2, 'descend');
    Left5_DiceAll(z, :, :) = squeeze(DiceAll(Left5, :, :));
    Right5_DiceAll(z, :, :) = squeeze(DiceAll(Right5, :, :));

    % Sort and store mean connectivity and true connectivity
    MeanConnectivityTrans = sort(MeanConnectivityTrans, 1);
    MConnectivity(z, :) = tanh(mean(MeanConnectivityTrans));
    STDConnectivity(z, :) = tanh(std(MeanConnectivityTrans));
    Left5_Connectivity(z, :) = tanh(MeanConnectivityTrans(Left5, :));
    Right5_Connectivity(z, :) = tanh(MeanConnectivityTrans(Right5, :));

    MeanContrueTrans = sort(MeanContrueTrans, 1);
    MContrue(z, :) = tanh(mean(MeanContrueTrans));
    STDContrue(z, :) = tanh(std(MeanContrueTrans));
    Left5_Contrue(z, :) = tanh(MeanContrueTrans(Left5, :));
    Right5_Contrue(z, :) = tanh(MeanContrueTrans(Right5, :));

    meanRel_smp_Trans = sort(meanRel_smp_Trans, 1);
    meanRel(z) = tanh(mean(meanRel_smp_Trans));
    stdRel(z) = tanh(std(meanRel_smp_Trans));
    Left5Rel(z) = tanh(meanRel_smp_Trans(Left5));
    Right5Rel(z) = tanh(meanRel_smp_Trans(Right5));

    meanRelAuto_smp_Trans = sort(meanRelAuto_smp_Trans, 1);
    meanRelAuto(z) = tanh(mean(meanRelAuto_smp_Trans));
    stdRelAuto(z) = tanh(std(meanRelAuto_smp_Trans));
    Left5RelAuto(z) = tanh(meanRelAuto_smp_Trans(Left5));
    Right5RelAuto(z) = tanh(meanRelAuto_smp_Trans(Right5));

end
