clear all; close all; clc;


addpath('SEED-G-toolbox-main\auxiliary functions');
addpath('SEED-G-toolbox-main\dependencies\BCT');
addpath('SEED-G-toolbox-main\dependencies\Fieldtrip');
addpath('SEED-G-toolbox-main\dependencies\MVGC');
addpath('SEED-G-toolbox-main\dependencies\NYH');
addpath('SEED-G-toolbox-main\dependencies\PDC_AsympSt\routines');
addpath('SEED-G-toolbox-main\main');

% savedir = 'G:\Shared drives\APR_HP_TeseM_CMUP\Mathematics_RAITSA\Code_v3\Random_Network\Data_gen\Generated_Data';

Nodes=13;

Trials=1;
ARn=0;  %
den=0.1;
Edges=round((Nodes*(Nodes-1))*den);

Conn_Range =[-0.9 0.9];

popt=3;
DelayRange=[1 popt];
Sw=eye(Nodes);% residual covariance matrix - White Noise (input Barnett Toolbox)
% modifica per lag zero
% Sw(2,4)=0.8;
% Sw(4,2)=0.8;

mtrunc=     0;       %default:0
decayfac=   100;     %default:100
Singtr=     1;
SigLim=     80;      %microvolt lim. for the amplitude of the generated signal
MinDelta=  0.1;
SNR=[inf];
DataLength=[10000];

iter=1;  % number of surrogate datasets

% name_1 ='MixSub_EEG';
% EEG1=loadname(fullfile(datadir,name_1));
% EEG1.samp(:,20,:)=[];

max_att=Trials*10;
% folder=sprintf('Nodes_%s',num2str(Nodes));
% mkdir(fullfile(savedir,folder));

for it=1:iter
    disp(it)
    flag=1;
    while flag==1
        
        %creating ground-truth network
        [CIJ inPos] = makerandCIJ_dir_withPureAR_onMainDiag(Nodes,Edges,ARn);
        
        [Model DelayMatrix]=get_ConnectivityModel(CIJ,Conn_Range, MinDelta, DelayRange);
        %creating corresponding delay matrix
        DelayMatrix = CIJ;
        DelayValues = min(DelayRange):1:max(DelayRange);
        DelayValues_pos = ceil((rand(1,Edges))*length(DelayValues));
        DelayMatrix(find(DelayMatrix)) = DelayValues(DelayValues_pos);
        
        %mixing model and its corresponding delay matrix
        [ModelDel] = rearrangeModel_forMVARbasedDataGen(DelayRange,Model,DelayMatrix);
        
        %generating autoregressive components
        for r=1:ARn
            
            for tt=1:size(EEG1.samp,3)
                EEG_ch = EEG1.samp(:,r,tt);
                orgEEG{tt}=EEG_ch;
            end %cycle on real trial
            clear tt
            
            ar_coef(:,r) = arfit_v2(orgEEG, popt, popt,'AIK','zero');
            
        end
        clear r
        
        %adding AR components to the ground truth
        for ii=1:ARn
            ModelDel(inPos(ii),inPos(ii),:)=ar_coef(:,ii);
        end
        clear ii
        
        for dl=1:length(DataLength)
            
            [Y,E,mtrunc]=var_to_tsdata(ModelDel,Sw,DataLength(dl),Singtr,mtrunc,decayfac);
            % check for stability
            [G,~] = var_to_autocov(ModelDel,Sw,1000);
            % disp(info.error)
            [Y1,E,mtrunc]=var_to_tsdata(ModelDel,eye(size(Model,1)),DataLength(dl),Singtr,mtrunc,decayfac);
%             clear G
%             [A,SIG] = tsdata_to_var(Y,3,'OLS');
%             [G,info] = var_to_autocov(ModelDel,Sw,100000);
%             [F] = autocov_to_pwcgc(G);
            
            

            
            Ctr=find(abs(Y)>SigLim);
%             Ctr1=find(abs(Y_nolag0)>SigLim);
            
            if ~isempty(Ctr) %|| ~isempty(Ctr1)
                flag=1;
                clear Data
                break
            end
            
            for s=1:length(SNR)
                for ch=1:size(Y,1)
                    nc=randn(1,size(Y,2));
                    Ynorm=norm(Y(ch,:));
                    ncnorm=norm(nc);
                    noise=(Ynorm/(ncnorm*sqrt(SNR(s))))*nc;
                    Ynoise(ch,:)=Y(ch,:)+noise;
                end
                samp{s,dl}=Ynoise';
           
            end
            clear Y Ynoise
            if dl==length(DataLength)
                flag=0;
            end
            
        end
        
    end
    
    % Data.Y_lag0=Y;
    Data.Y_nolag0=Y1;
    Data.Model=Model;
    Data.SNR=SNR;
    Data.Sw_lag0=Sw;
    Data.Sw_nolag0=eye(size(Model,1));
    Data.DataLength=DataLength;
    Data.ConnRange=Conn_Range;
    Data.popt=popt;
    Data.ModelDel=ModelDel;
    % name=sprintf('Sim_%s',num2str(it));
    % disp(name)
    % savedir_1=fullfile(savedir,folder);
    % save(fullfile(savedir_1,name),'Data');
    % clear Data
    save Model_Struct_13_Nodes_den_0.05.mat Data
end