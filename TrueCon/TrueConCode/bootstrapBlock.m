
%Copyright © 2024 Koten and Schüppen All rights reserved
%Important Notice: This code is not intended for medical applications 
%and does not have legal approval for such use. We strongly recommend 
%using FDA-approved software for any medical purposes. 




% This function is designed to execute a stationary block bootstrap for estimating the confidence intervals of reliability (Rel_trans) using three input arguments.
% Additionally, it estimates connectivity (ConTestBoot_trans for testing and ConRetestBoot_trans for retesting), as well as detectable connectivity (contrueTestBoot_trans, contrueTestBoot_trans) using five input arguments.
% If the user intends to include all the mentioned measures, it is required to provide nine input arguments.

function [conf] = bootstrapBlock(numboot,RoiSeedTest,RoiSeedRetest,RoiTargetTest,RoiTargetRetest,RoiSeedRelTest,RoiSeedRelRetest,RoiTargetRelTest,RoiTargetRelRetest)

if nargin==9

    % estaimte optimal block length
    block= [opt_block(RoiSeedTest),opt_block(RoiSeedRetest),opt_block(RoiTargetTest),opt_block(RoiTargetRetest),opt_block(RoiSeedRelTest),opt_block(RoiSeedRelRetest),...
        opt_block(RoiTargetRelTest),opt_block(RoiTargetRelRetest)];

    % estimate the mean block length that serves all
    W = ceil(mean(block(1,:),2));
    % obtain indices for block bootstrap
    [RoiSeedT, indices]=stationary_bootstrap(RoiSeedTest,numboot,W);
    % apply indices
    RoiSeedR=RoiSeedRetest(indices);
    RoiTargetT=RoiTargetTest(indices);
    RoiTargetR=RoiTargetRetest(indices);

    RoiSeedRelT=RoiSeedRelTest(indices);
    RoiSeedRelR=RoiSeedRelRetest(indices);
    RoiTargetRelT=RoiTargetRelTest(indices);
    RoiTargetRelR=RoiTargetRelRetest(indices);

    % Obtain connectivity for test and retest
    ConTestBoot_trans = atanh(diag(corr(RoiSeedT,RoiTargetT)));
    ConRetestBoot_trans = atanh(diag(corr(RoiSeedR,RoiTargetR)));
    % Obtain reliabilty for seed and target
    RelSeed_trans=atanh(diag(corr(RoiSeedRelT,RoiSeedRelR)));
    RelTarget_trans=atanh(diag(corr(RoiTargetRelT,RoiTargetRelR)));

    % Estimate detectable connectivity
    [~,contrueTestBoot_trans,~,contrueRetestBoot_trans,relmaxBoot_trans]=TrueCon(ConTestBoot_trans, ConRetestBoot_trans,RelSeed_trans,RelTarget_trans);

    % Gather confidence intervals for ConTestBoot_trans
    confidence=ConfInterval(ConTestBoot_trans);
    conf.ConTestBoot_trans_LB=confidence.LB;
    conf.ConTestBoot_trans_UB=confidence.UB;

    % Gather confidence intervals for ConRetestBoot_trans
    confidence=ConfInterval(ConRetestBoot_trans);
    conf.ConRetestBoot_trans_LB=confidence.LB;
    conf.ConRetestBoot_trans_UB=confidence.UB;

    % Gather confidence intervals for contrueTestBoot_trans
    confidence=ConfInterval(contrueTestBoot_trans);
    conf.contrueTestBoot_Trans_LB=confidence.LB;
    conf.contrueTestBoot_Trans_UB=confidence.UB;

    % Gather confidence intervals for contrueRetestBoot_trans
    confidence=ConfInterval(contrueRetestBoot_trans);
    conf.contrueRetestBoot_Trans_LB=confidence.LB;
    conf.contrueRetestBoot_Trans_UB=confidence.UB;

    % Gather confidence intervals for relmaxBoot_trans
    confidence=ConfInterval(relmaxBoot_trans);
    conf.relmaxBoot_Trans_LB=confidence.LB;
    conf.relmaxBoot_Trans_UB=confidence.UB;

elseif nargin==5


    % Estimate optimal blocksize using opt_block
    block = [opt_block(RoiSeedTest),opt_block(RoiSeedRetest),opt_block(RoiTargetTest),opt_block(RoiTargetRetest)];
    % Estimate the mean block length that serves all
    W = ceil(mean(block(1,:),2));
    % Obtain indices for block bootstrap
    [RoiSeedT, indices]=stationary_bootstrap(RoiSeedTest,numboot,W);
    % Apply indices
    RoiSeedR=RoiSeedRetest(indices);
    RoiTargetT=RoiTargetTest(indices);
    RoiTargetR=RoiTargetRetest(indices);

    % Obtain connectivity for test and retest
    ConTestBoot_trans = atanh(diag(corr(RoiSeedT,RoiTargetT)));
    ConRetestBoot_trans = atanh(diag(corr(RoiSeedR,RoiTargetR)));
    % Obtain reliabilty for seed and target
    RelSeed_trans=atanh(diag(corr(RoiSeedT,RoiSeedR)));
    RelTarget_trans=atanh(diag(corr(RoiTargetT,RoiTargetR)));

    % Estimate detectable connectivity
    [~,contrueTestBoot_trans,~,contrueRetestBoot_trans,relmaxBoot_trans]=TrueCon(ConTestBoot_trans, ConRetestBoot_trans,RelSeed_trans,RelTarget_trans);

    % Gather confidence intervals for ConTestBoot_trans
    confidence=ConfInterval(ConTestBoot_trans);
    conf.ConTestBoot_trans_LB=confidence.LB;
    conf.ConTestBoot_trans_UB=confidence.UB;

    % Gather confidence intervals for ConRetestBoot_trans
    confidence=ConfInterval(ConRetestBoot_trans);
    conf.ConRetestBoot_trans_LB=confidence.LB;
    conf.ConRetestBoot_trans_UB=confidence.UB;

    % Gather confidence intervals for contrueTestBoot_trans
    confidence=ConfInterval(contrueTestBoot_trans);
    conf.contrueTestBoot_Trans_LB=confidence.LB;
    conf.contrueTestBoot_Trans_UB=confidence.UB;

    % Gather confidence intervals for contrueRetestBoot_trans
    confidence=ConfInterval(contrueRetestBoot_trans);
    conf.contrueRetestBoot_Trans_LB=confidence.LB;
    conf.contrueRetestBoot_Trans_UB=confidence.UB;

    % Gather confidence intervals for relmaxBoot_trans
    confidence=ConfInterval(relmaxBoot_trans);
    conf.relmaxBoot_Trans_LB=confidence.LB;
    conf.relmaxBoot_Trans_UB=confidence.UB;

else


    % Estimate optimal blocksize using opt_block
    block = [opt_block(RoiSeedTest),opt_block(RoiSeedRetest)];
    W = ceil(mean(block(1,:),2));
    % Obtain indices for block bootstrap
    [RoiSeedT, indices]=stationary_bootstrap(RoiSeedTest,numboot,W);
    % Apply indices
    RoiSeedR=RoiSeedRetest(indices);
    % Obtain reliabilty for seed
    Rel_trans=atanh(diag(corr(RoiSeedT,RoiSeedR)));

    Rel=sort(Rel_trans);
    indi=Rel==0;
    Rel(indi)=[];
    n=length(Rel);
    Li=round(2.5./100.*(n+1));
    Ui=round(97.5./100.*(n+1));
    conf.Rel_Trans_LB=Rel(Li);
    conf.Rel_Trans_UB=Rel(Ui);

end

end