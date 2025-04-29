function [Amt,Sut,ki]=VARreorder(Am,Su,ii,jj)

% extract VAR parameters from Am
M=size(Am,1);
p=size(Am,2)/M;
for i=1:p
    Ak(:,:,i)=Am(:,M*(i-1)+1:M*i);
end

%reorganize to have ii as 1st process and jj as 2nd process
kk=setdiff(1:M,[ii jj]);
ki=[ii jj kk]; %new order (1st ii, 2nd jj)

%permutation according to causal order ki
P=zeros(M,M); 
for i=1:M
    P(i,ki(i))=1;
end
Sut=P*Su*P'; %Su tilda (reordered)
for i=1:p
    Akt(:,:,i)=P*Ak(:,:,i)*P';
end
Amt=[]; % group coefficient blocks Ak in a single matrix Am
for i=1:p
    Amt=[Amt Akt(:,:,i)];
end
end

