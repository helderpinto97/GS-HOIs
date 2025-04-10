%% OIR - computation of MIR and of its time-domain decomposition and frequency-domain expansion for two blocks of processes

function out=oir_mir(A,C,K,Su,Mv,i_1,i_2,Fs,nfft,do_td)

narginchk(7,10)
if nargin<10, do_td='n'; end
if nargin<9, nfft=512; end
if nargin<8, Fs=1; end %default non compute time domain measures through submodels

% indexes of the two blocks to analyze inside the Q time series
[i1,i2]=oir_subindexes(Mv,i_1,i_2);

% reduced model with the two blocks to analyze - ired=[i1 i2]
[CR,KR,VR]=oir_submodel(A,C,K,Su,[i1 i2]); % Eqs. (27-28)

% indexes of the two blocks inside the reduced process
Mr1=length(i1); Mr2=length(i2); Mr=Mr1+Mr2;
j1=1:Mr1; j2=Mr1+1:Mr;

% FREQUENCY DOMAIN GC
[f12,f1_2,f2_1,f1o2,S,freq] = oir_fdGC(A,CR,KR,VR,j1,j2,Fs,nfft);

if do_td=='y'
    % TIME DOMAIN GC
    [~,~,VR1]=oir_submodel(A,C,K,Su,i1); % reduced model with the 1st block only - ired=i1
    F2_1 = log(det(VR1)) - log(det(VR(j1,j1)));

    [~,~,VR2]=oir_submodel(A,C,K,Su,i2); % reduced model with the 2nd block only - ired=i2
    F1_2 = log(det(VR2)) - log(det(VR(j2,j2)));

    F1o2=log(det(VR(j1,j1))) + log(det(VR(j2,j2))) - log(det(VR));
    F12=F1_2+F2_1+F1o2;
else
    F2_1=mean(f2_1);
    F1_2=mean(f1_2); %GC from i1 to i2
    F1o2=mean(f1o2); %GC from i2 to i1
    F12=mean(f12);
end


out.I12=F12/2;
out.T1_2=F1_2/2; %TE from i1 to i2
out.T2_1=F2_1/2; %TE from i2 to i1
out.I1o2=F1o2/2;
out.f12=f12;
out.f1_2=f1_2;
out.f2_1=f2_1;
out.f1o2=f1o2;
out.S=S;
out.freq=freq;


end