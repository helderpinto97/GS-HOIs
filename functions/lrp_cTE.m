%% conditional Granger causality using LRP toolbox
% (single regression, solution of Yule Walker equations for correlations up to lag q)
% computation of GC from Xi to Xj given Xk
% input: VAR parameters Am, Su (with idMVAR)

function [TE,SigmaY,SigmaY_YZ,SigmaY_XYZ] = lrp_cTE(Am,Su,i,j,k,q)

[SigmaY,SigmaY_YZ] = lrp_LinReg(Am,Su,q,j,[j k]);
[~,SigmaY_XYZ] = lrp_LinReg(Am,Su,q,j,[i j k]);

TE=0.5*log(SigmaY_YZ/SigmaY_XYZ);

end