function  OIR_delta=OIR_lin(Am,Su,q,ii)

% Am -- Coeficients Matrix
% Su -- Covariance Matrix
% q -- Truncation lag for correlation
% ii -- indexes of Processes (the last element of ii must be the element that we are test to be included in X_{k-1})

% Save the X_{k-1}
x_k_minus_one=setdiff(ii,ii(end));

% Value of k
k=length(ii)-2;

% Initiate the value of the sum
Mir_Ind=0;
for i=1:length(x_k_minus_one)
    x_temp=setdiff(x_k_minus_one,x_k_minus_one(i));
    Ixk=MIR_lin(Am,Su,q,ii(end),x_temp);
    Mir_Ind=Mir_Ind+Ixk;
end

% Joint MIR X_j with X_{k-1}
mir_all=MIR_lin(Am,Su,q,ii(end),x_k_minus_one);

OIR_delta=Mir_Ind-k*mir_all;

end