function iir=iir_lin(Am,Su,q,ii,jj,kk)

% Am - Coeficients Matrix 
% Su - Covariance Matrix
% q - Truncation for Correlations
% ii - Target index (y)
% jj, kk - Indexes of the remaining processes

% 1st Step - Individuals MIRs
Iyx = MIR_lin(Am,Su,q,ii,jj);
Iyz=MIR_lin(Am,Su,q,ii,kk);

% 2nd Step - Joint MIR
Iy_xz=MIR_lin(Am,Su,q,ii,[jj kk]);

% 3rd Step - IIR
iir= Iyx+Iyz-Iy_xz;
end