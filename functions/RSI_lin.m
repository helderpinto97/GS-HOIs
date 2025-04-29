function rsi=RSI_lin(Am,Su,q,ii)

% Am -- Coeficients Matrix
% Su -- Covariance Matrix
% q -- Truncation lag for correlation
% ii -- indexes of Processes (the last element of ii must be the element that we are test to be included in X_{k-1})

% Save the X_{k-1}
proc_rem=setdiff(ii,ii(end));

Mir_Ind=0;
for k=1:length(proc_rem)
    Ixk=MIR_lin(Am,Su,q,ii(end),proc_rem(k));
    Mir_Ind=Mir_Ind+Ixk;
end

mir_all=MIR_lin(Am,Su,q,ii(end),proc_rem);

rsi=Mir_Ind-mir_all;

end