function Ixyz_W=cIIR_lin(Am,Su,q,ii,jj)

% Am - Coefficients Matrix
% Su - Covariance Matrix
% q - Truncation lag Correlations
% ii - Indexes Triplet for IIR (so vector with 3 entries, initial triplet)
% jj - Indexes for Conditioning (can be vector)

ret_12_z=cMIR_lin(Am,Su,q,ii(1),ii(2),jj);
ret_13_z=cMIR_lin(Am,Su,q,ii(1),ii(3),jj);
ret_1_23_z=cMIR_lin(Am,Su,q,ii(1),[ii(2) ii(3)],jj);

Ixyz_W=ret_12_z+ret_13_z-ret_1_23_z;

end