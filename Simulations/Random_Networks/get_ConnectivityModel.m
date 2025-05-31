function [Model DelayMatrix MaxDelay] = get_ConnectivityModel(ModelStructure, ValuesRange, MinDelta, DelayRange)

%%%% input:
% ModelStructure: Connectivity Model Structure. It is a NxN (N: number of nodes)
%                 matrix in which entries are equal to 1 in correspondence to a link
% ValuesRange:    range of possible connections values
% MinDelta:       minimum difference between two different values of connections
% DelayRange:     range of possible delay values

%%%% output:
% Model:          NxN connectivity model with connections values belonging to the
%                 selected range (ValuesRange)
% DelayMatrix:    NxN matrix in which the entry ij (i,j=1,...N) different from zero
%                 corresponds to the delay associated with the connection
%                 ij of Model. Delay values belong to the selected range
%                 DelayRange.
% MaxDelay:       maximum possible delay

% Author: Manuela Petti
% Date: 05/06/2014
% version: 1.0.0

% modified 29/05/2015: to allow negative connection values in ValuesRange

Structure = ModelStructure;
N_conn = length(find(ModelStructure));
nNodes = size(Structure,1);

%%% assignment of values to the connections
ModelValues =min(ValuesRange):MinDelta:max(ValuesRange);

if any(ModelValues==0)
    x=find(ModelValues==0);
    ModelValues(x)=[];
end

ModelValues_pos =ceil((rand(1,N_conn))*length(ModelValues));
ModelStructure(find(ModelStructure)) = ModelValues(ModelValues_pos);
Model = ModelStructure;

%%% assignment of delay values to the connections
DelayValues = min(DelayRange):1:max(DelayRange);
DelayValues_pos = ceil((rand(1,N_conn))*length(DelayValues));
Structure(find(Structure)) = DelayValues(DelayValues_pos);
DelayMatrix = Structure;

%%% computation of maximum delay
diagSecondarie = zeros([nNodes-1,1]);
for rr=1:nNodes
    for cc=(rr+1):nNodes
        diagSecondarie(cc-rr) = diagSecondarie(cc-rr)+DelayMatrix(rr,cc);
    end
end
MaxDelay = max(diagSecondarie);
