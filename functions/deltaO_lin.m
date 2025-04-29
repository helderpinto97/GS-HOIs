function delta=deltaO_lin(Am,Su,q,ij,j)

% Faes, L., Mijatovic, G., Antonacci, Y., Pernice, R., Bara, C., Sparacino, L., Sammartino, M., Porta, A., Marinazzo, D., & Stramaglia, S. (2022). 
% A New Framework for the Time-and Frequency-Domain Assessment of High-Order Interactions in Networks of Random Processes. IEEE Transactions on Signal Processing, 70, 5766–5777. 
% https://doi.org/10.1109/TSP.2022.3221892

% Am - Coefficents Matrix (N x p*N)
% Su - Covariance Matrix  (N x N)
% q - Truncantion lag for correlation
% ij - complete vector of the processesa
% j - target process

% Number of processes
N=size(Am,1);

if N<=2
    error('Mininum of 3 processes to estimate this measure');
end

% First Term of the formula - Eq.(8)
I_t=MIR_lin(Am,Su,q,j,ij);

% Sum Terms - Eq. (8)
I_sum=0;
for k=1:N-1
    j_tmp=setdiff(ij,k); 
    tmp=MIR_lin(Am,Su,q,N,j_tmp);
    I_sum=I_sum+tmp;
end

delta=(2-N)*I_t+I_sum;

end